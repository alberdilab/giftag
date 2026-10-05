"""Summarize the profile counts and disk footprint of full and giftag libraries."""

import argparse
import csv
import json
from pathlib import Path

from compare import write_tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    project = args.project
    db = project / "data/giftag-db"
    upstream = project / "data/upstream"
    manifest = json.loads((db / "giftag.json").read_text())
    with (db / "profiles.tsv").open(newline="") as handle:
        profiles = list(csv.DictReader(handle, delimiter="\t"))
    counts = {source: len({r["profile"] for r in profiles if r["source"] == source})
              for source in ("kofam", "ncbifam", "pfam", "dbcan", "dbcan_sub")}
    source_meta = manifest["sources"]
    full_kofam = sum(line.rstrip().endswith(".hmm") for line in
                     (project / "provenance/kofam-archive-members.txt").open())
    rows = []

    def add(collection, scope, count, targeted, path, artifact_format="HMM library"):
        rows.append(dict(collection=collection, scope=scope, profiles_searched=count,
                         gifter_target_profiles=targeted,
                         stored_bytes=path.stat().st_size,
                         artifact_format=artifact_format, file=str(path)))

    add("KOfam", "full_archive", full_kofam, counts["kofam"],
        upstream / "kofam/profiles.tar.gz", "compressed tarball")
    add("KOfam", "giftag_subset", counts["kofam"], counts["kofam"],
        db / "kofam.hmm", "concatenated HMM")
    add("dbCAN_family", "full", source_meta["dbcan"]["z_family"], counts["dbcan"],
        upstream / "dbcan/dbCAN.hmm")
    add("dbCAN_family", "giftag", source_meta["dbcan"]["z_family"], counts["dbcan"],
        db / "dbcan.hmm")
    add("dbCAN_family", "run_dbcan_subset", source_meta["dbcan"]["z_family"],
        counts["dbcan"], upstream / "dbcan-subset/dbCAN.hmm")
    add("dbCAN_sub", "full", source_meta["dbcan"]["z_sub"], counts["dbcan_sub"],
        upstream / "dbcan/dbCAN_sub.hmm")
    add("dbCAN_sub", "giftag_subset", source_meta["dbcan"]["sub_profiles_kept"],
        counts["dbcan_sub"], db / "dbcan_sub.hmm")
    add("dbCAN_sub", "run_dbcan_subset", source_meta["dbcan"]["sub_profiles_kept"],
        counts["dbcan_sub"], upstream / "dbcan-subset/dbCAN_sub.hmm")
    add("NCBIFAM", "giftag_subset", counts["ncbifam"], counts["ncbifam"],
        db / "ncbifam.hmm")
    add("Pfam", "giftag_subset", counts["pfam"], counts["pfam"], db / "pfam.hmm")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output, ("collection", "scope", "profiles_searched",
                            "gifter_target_profiles", "stored_bytes",
                            "artifact_format", "file"), rows)


if __name__ == "__main__":
    main()
