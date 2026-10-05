"""Normalize reference HMM calls and compare them with giftag marker calls.

The unit of agreement is (genome, protein, marker), not a whole-genome count.
The reference is a comparator, not ground truth. Run on identical protein FASTA
files and matching source releases before interpreting any discordance.
"""

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


FIELDS = ("genome_id", "gene_id", "namespace", "accession")
SOURCE_NAMESPACES = {
    "kofam": {"KO"}, "dbcan_family": {"CAZY"},
    "dbcan_sub": {"CAZY"}, "pfam": {"PFAM"},
    "ncbifam": {"NCBIFAM", "TIGRFAM"},
}


def read_tsv(path):
    with open(path, newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path, columns, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def searchable_markers(path):
    return {(r["namespace"], r["accession"]) for r in read_tsv(path)
            if r["status"] == "searchable"}


def _clean_accession(value):
    return re.sub(r"\.\d+$", "", value.strip())


def normalize(path, fmt, genome, markers):
    """Return projected gifter markers from one upstream output file.

    InterProScan's NCBIFAM accessions can be versioned; use the marker ledger
    to resolve the versionless TIGRFAM alias to the correct gifter namespace.
    dbCAN family submodels imply their parent family, matching giftag's output.
    """
    calls = set()
    if fmt == "kofam-mapper":
        with open(path, newline="") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0] and re.fullmatch(r"K\d{5}", parts[1]):
                    calls.add((genome, parts[0], "KO", parts[1]))
    elif fmt in {"hmmsearch-pfam", "hmmsearch-ncbifam"}:
        with open(path) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split()
                if len(fields) < 18:
                    raise ValueError(f"{path}: malformed hmmsearch domtblout line")
                gene, accession = fields[0], _clean_accession(fields[4])
                namespaces = ({"PFAM"} if fmt == "hmmsearch-pfam" else
                              {"NCBIFAM", "TIGRFAM"} if accession.startswith("TIGR") else
                              {"NCBIFAM"})
                for marker in markers:
                    if marker[0] in namespaces and _clean_accession(marker[1]) == accession:
                        calls.add((genome, gene, *marker))
    elif fmt == "interproscan":
        with open(path, newline="") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) < 5 or not row[0]:
                    continue
                analysis, accession = row[3].upper(), _clean_accession(row[4])
                if analysis == "PFAM":
                    namespaces = {"PFAM"}
                elif analysis in {"NCBIFAM", "TIGRFAM"}:
                    namespaces = {"NCBIFAM", "TIGRFAM"} if accession.startswith("TIGR") else {"NCBIFAM"}
                else:
                    continue
                for marker in markers:
                    if marker[0] in namespaces and _clean_accession(marker[1]) == accession:
                        calls.add((genome, row[0], *marker))
    elif fmt in {"dbcan-family", "dbcan-sub"}:
        # run_dbcan V5 writes a named TSV. Match column names rather than
        # positions because the substrate columns differ between releases.
        for row in read_tsv(path):
            lower = {k.lower().strip().replace(" ", "_"): v for k, v in row.items() if k}
            gene = next((lower[k] for k in ("gene_id", "geneid", "gene", "target_name")
                         if lower.get(k)), None)
            profile = next((lower[k] for k in ("hmm_name", "subfam_name", "hmm_profile",
                                                 "hmmprofile", "hmm", "dbcan_sub",
                                                 "dbcan_subfam", "subfamily", "family")
                            if lower.get(k)), None)
            if not gene or not profile:
                raise ValueError(f"{path}: unrecognized run_dbcan columns: {list(row)}")
            profile = profile.split("|")[0].removesuffix(".hmm")
            if fmt == "dbcan-sub":
                match = re.search(r"(?:AA|CBM|CE|GH|GT|PL)\d+_e\d+", profile)
            else:
                match = re.search(r"(?:AA|CBM|CE|GH|GT|PL)\d+(?:_\d+)?", profile)
            if not match:
                continue
            accession = match.group()
            if fmt == "dbcan-family" and "_" in accession:
                if ("CAZY", accession) in markers:
                    calls.add((genome, gene, "CAZY", accession))
                accession = accession.split("_")[0]
            calls.add((genome, gene, "CAZY", accession))
    else:
        raise ValueError(fmt)
    return {call for call in calls if call[2:] in markers}


def compare(giftag, reference, markers, source):
    allowed = SOURCE_NAMESPACES[source]
    if source == "dbcan_family":
        selected = {m for m in markers if m[0] == "CAZY" and "_e" not in m[1]}
    elif source == "dbcan_sub":
        selected = {m for m in markers if m[0] == "CAZY" and "_e" in m[1]}
    else:
        selected = {m for m in markers if m[0] in allowed}
    a = {tuple(r[k] for k in FIELDS) for r in read_tsv(giftag)
         if (r["namespace"], r["accession"]) in selected}
    b = {tuple(r[k] for k in FIELDS) for r in read_tsv(reference)
         if (r["namespace"], r["accession"]) in selected}
    shared = a & b
    only_a, only_b = a - b, b - a
    summary = {
        "source": source, "searchable_markers": len(selected),
        "giftag_pairs": len(a), "reference_pairs": len(b),
        "shared_pairs": len(shared), "giftag_only_pairs": len(only_a),
        "reference_only_pairs": len(only_b),
        "jaccard": len(shared) / len(a | b) if a | b else None,
        "reference_recall": len(shared) / len(b) if b else None,
        "reference_precision": len(shared) / len(a) if a else None,
    }
    differences = [dict(zip(FIELDS, call), status=status)
                   for status, calls in (("giftag_only", only_a), ("reference_only", only_b))
                   for call in calls]
    differences.sort(key=lambda r: tuple(r[k] for k in FIELDS))
    by_genome = defaultdict(Counter)
    by_marker = defaultdict(Counter)
    for marker in selected:
        by_marker[marker]
    for status, calls in (("shared", shared), ("giftag_only", only_a), ("reference_only", only_b)):
        for call in calls:
            by_genome[call[0]][status] += 1
            by_marker[call[2:]][status] += 1
    def rows(counts, names):
        return [dict(zip(names, key), **{k: counter[k] for k in
                ("shared", "giftag_only", "reference_only")})
                for key, counter in sorted(counts.items())]
    genomes = rows({(k,): v for k, v in by_genome.items()}, ("genome_id",))
    marker_rows = rows(by_marker, ("namespace", "accession"))
    return summary, differences, genomes, marker_rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    norm = sub.add_parser("normalize", help="convert one upstream output to gifter marker pairs")
    norm.add_argument("--format", required=True,
                      choices=("kofam-mapper", "interproscan", "dbcan-family", "dbcan-sub",
                               "hmmsearch-pfam", "hmmsearch-ncbifam"))
    norm.add_argument("--input", type=Path, required=True)
    norm.add_argument("--genome", required=True)
    norm.add_argument("--markers", type=Path, required=True)
    norm.add_argument("--output", type=Path, required=True)
    comp = sub.add_parser("compare", help="write agreement tables for one source")
    comp.add_argument("--giftag", type=Path, required=True)
    comp.add_argument("--reference", type=Path, required=True)
    comp.add_argument("--markers", type=Path, required=True)
    comp.add_argument("--source", choices=tuple(SOURCE_NAMESPACES), required=True)
    comp.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    markers = searchable_markers(args.markers)
    if args.command == "normalize":
        calls = sorted(normalize(args.input, args.format, args.genome, markers))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        write_tsv(args.output, FIELDS, [dict(zip(FIELDS, call)) for call in calls])
    else:
        args.outdir.mkdir(parents=True, exist_ok=True)
        summary, differences, genomes, marker_rows = compare(
            args.giftag, args.reference, markers, args.source)
        (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        write_tsv(args.outdir / "discordant.tsv", FIELDS + ("status",), differences)
        write_tsv(args.outdir / "by_genome.tsv",
                  ("genome_id", "shared", "giftag_only", "reference_only"), genomes)
        write_tsv(args.outdir / "by_marker.tsv",
                  ("namespace", "accession", "shared", "giftag_only", "reference_only"), marker_rows)
        print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
