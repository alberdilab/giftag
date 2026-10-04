"""Read genome or protein FASTA and call genes with pyrodigal."""

import gzip
import re
from pathlib import Path

import pyrodigal

from giftag import GiftagError

# Prodigal's own advice: train on the genome when there is enough sequence to
# learn its gene model, otherwise use the pre-trained metagenomic models.
SINGLE_MODE_MIN_BP = 100_000

_FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fas", ".faa", ".ffn", ".fsa", ".seq")
_NUCLEOTIDE = set("ACGTUN")


def genome_id(path):
    name = Path(path).name
    if name.endswith(".gz"):
        name = name[:-3]
    for suffix in _FASTA_SUFFIXES:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return name


def expand_inputs(paths):
    """Replace each directory by the FASTA files directly inside it."""
    files = []
    for path in map(Path, paths):
        if path.is_dir():
            found = sorted(
                child for child in path.iterdir() if child.is_file()
                and (child.name[:-3] if child.name.endswith(".gz") else child.name)
                .lower().endswith(_FASTA_SUFFIXES))
            if not found:
                raise GiftagError(f"{path}: no FASTA files in this directory")
            files.extend(found)
        elif not path.exists():
            raise GiftagError(f"{path}: no such file")
        else:
            files.append(path)
    return files


def read_fasta(path):
    """Return `[(id, sequence)]`, reading gzip transparently."""
    with open(path, "rb") as probe:
        gzipped = probe.read(2) == b"\x1f\x8b"
    opener = gzip.open if gzipped else open
    records = []
    name = None
    chunks = []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                header = line[1:].strip()
                if not header:
                    raise GiftagError(f"{path}: a FASTA record has no identifier")
                name = header.split()[0]
                chunks = []
            elif name is not None:
                chunks.append(re.sub(r"\s+", "", line))
            elif line.strip():
                raise GiftagError(f"{path}: not a FASTA file")
    if name is not None:
        records.append((name, "".join(chunks)))
    if not records:
        raise GiftagError(f"{path}: no sequences")
    names = [record[0] for record in records]
    if len(set(names)) != len(names):
        raise GiftagError(f"{path}: duplicate sequence identifiers")
    return records


def looks_nucleotide(records):
    sample = "".join(seq for _, seq in records)[:200_000].upper()
    if not sample:
        return False
    return sum(char in _NUCLEOTIDE for char in sample) / len(sample) >= 0.9


def call_genes(records, mode="auto"):
    """Predict proteins; return `(proteins, mode_used)` with Prodigal-style IDs."""
    total = sum(len(seq) for _, seq in records)
    if mode == "auto":
        mode = "single" if total >= SINGLE_MODE_MIN_BP else "meta"
    if mode == "single":
        if total < 20_000:
            raise GiftagError("single mode needs at least 20,000 bp to train; use --mode meta")
        finder = pyrodigal.GeneFinder(meta=False)
        finder.train(*(seq.encode() for _, seq in records))
    else:
        finder = pyrodigal.GeneFinder(meta=True)
    proteins = []
    for contig, seq in records:
        for index, gene in enumerate(finder.find_genes(seq.encode()), start=1):
            proteins.append((f"{contig}_{index}", gene.translate(include_stop=False)))
    return proteins, mode


def write_fasta(path, records):
    with open(path, "w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for start in range(0, len(seq), 60):
                handle.write(seq[start:start + 60] + "\n")
