"""Record whether each target NCBIFAM model exists in InterPro member releases."""

import argparse
import re
from pathlib import Path

from compare import read_tsv, write_tsv


def base(accession):
    return re.sub(r"\.\d+$", "", accession)


def accessions(path, wanted):
    found = set()
    with open(path, "rb") as handle:
        for line in handle:
            if line.startswith(b"ACC "):
                parts = line.split()
                if len(parts) >= 2:
                    accession = base(parts[1].decode())
                    if accession in wanted:
                        found.add(accession)
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--markers", type=Path, required=True)
    parser.add_argument("--giftag-hmm", type=Path, required=True)
    parser.add_argument("--interpro5-hmm", type=Path, required=True)
    parser.add_argument("--interpro6-hmm", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    markers = sorted({(r["namespace"], r["accession"]) for r in read_tsv(args.markers)
                      if r["status"] == "searchable" and
                      r["namespace"] in {"NCBIFAM", "TIGRFAM"}})
    wanted = {base(accession) for _, accession in markers}
    libraries = {
        "giftag20": accessions(args.giftag_hmm, wanted),
        "interpro5_17": accessions(args.interpro5_hmm, wanted),
        "interpro6_18": accessions(args.interpro6_hmm, wanted),
    }
    rows = [dict(namespace=namespace, accession=accession,
                 giftag20=int(base(accession) in libraries["giftag20"]),
                 interpro5_17=int(base(accession) in libraries["interpro5_17"]),
                 interpro6_18=int(base(accession) in libraries["interpro6_18"]))
            for namespace, accession in markers]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output, ("namespace", "accession", "giftag20",
                            "interpro5_17", "interpro6_18"), rows)
    print({name: sum(row[name] for row in rows) for name in libraries})


if __name__ == "__main__":
    main()
