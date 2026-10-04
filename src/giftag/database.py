"""Load a database directory written by `giftag build`."""

import csv
import json
import re
from collections import defaultdict
from pathlib import Path

import pyhmmer

from giftag import DB_FORMAT, GiftagError
from giftag import search as search_
from giftag.build import cluster_family, dbcan_cluster_key

_NOT_AMINO = re.compile(r"[^ACDEFGHIKLMNPQRSTVWYBJZOUX*]")


class Database:
    def __init__(self, path):
        self.path = Path(path)
        manifest = self.path / "giftag.json"
        if not manifest.exists():
            raise GiftagError(f"{self.path}: no giftag database here; run `giftag build` first")
        with open(manifest) as handle:
            self.manifest = json.load(handle)
        if self.manifest.get("format") != DB_FORMAT:
            raise GiftagError(
                f"{self.path}: database format {self.manifest.get('format')} is not the "
                f"format {DB_FORMAT} this giftag reads; rebuild it with `giftag build --force`")

        self.rules = {}
        self.labels = defaultdict(list)
        with open(self.path / "profiles.tsv", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = (row["source"], row["profile"])
                self.rules[key] = (row["rule"], float(row["threshold"]))
                self.labels[key].append((row["namespace"], row["accession"]))
        with open(self.path / "markers.tsv", newline="") as handle:
            self.markers = list(csv.DictReader(handle, delimiter="\t"))

        self.alphabet = pyhmmer.easel.Alphabet.amino()
        self._hmms = {}

    @property
    def unsearchable(self):
        return [m for m in self.markers if m["status"] != "searchable"]

    def _load(self, name):
        if name not in self._hmms:
            path = self.path / f"{name}.hmm"
            hmms = []
            if path.exists() and path.stat().st_size:
                with pyhmmer.plan7.HMMFile(path) as handle:
                    hmms = list(handle)
            self._hmms[name] = hmms
        return self._hmms[name]

    def _sub_by_family(self):
        if "dbcan_sub_grouped" not in self._hmms:
            grouped = defaultdict(list)
            for hmm in self._load("dbcan_sub"):
                grouped[cluster_family(dbcan_cluster_key(hmm.name))].append(hmm)
            self._hmms["dbcan_sub_grouped"] = dict(grouped)
        return self._hmms["dbcan_sub_grouped"]

    def digitize(self, proteins):
        sequences = []
        for name, seq in proteins:
            seq = _NOT_AMINO.sub("X", seq.upper().rstrip("*"))
            if seq:
                text = pyhmmer.easel.TextSequence(name=name, sequence=seq)
                sequences.append(text.digitize(self.alphabet))
        return pyhmmer.easel.DigitalSequenceBlock(self.alphabet, sequences)

    def search(self, proteins, cpus=0, gate=False, stage=None):
        """Return every accepted `Call` for these `(id, sequence)` proteins.

        `stage`, if given, is called with the name of each source as its search
        starts."""
        stage = stage or (lambda name: None)
        block = self.digitize(proteins)
        if not len(block):
            return []
        calls = []
        if self._load("kofam"):
            stage("KOfam")
            calls += search_.search_kofam(self._load("kofam"), block, self.rules, cpus)
        if self._load("ncbifam"):
            stage("NCBIfam")
            calls += search_.search_cutoff("ncbifam", self._load("ncbifam"), block,
                                           self.rules, cpus, "trusted")
        if self._load("pfam"):
            stage("Pfam")
            calls += search_.search_cutoff("pfam", self._load("pfam"), block,
                                           self.rules, cpus, "gathering")
        if self._load("dbcan"):
            stage("dbCAN")
            dbcan = self.manifest["sources"]["dbcan"]
            calls += search_.search_dbcan(self._load("dbcan"), self._sub_by_family(), block,
                                          cpus, dbcan["z_family"], dbcan.get("z_sub", 1), gate=gate)
        return calls

    def load_all(self):
        """Read every profile now, so the first genome is not charged for it."""
        for name in ("kofam", "ncbifam", "pfam", "dbcan"):
            self._load(name)
        self._sub_by_family()
