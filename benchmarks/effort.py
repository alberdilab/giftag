"""Compare giftag with serial sums of equivalent reference search stages.

The reference stages were measured separately. Their sum is an estimate of a
serial workflow, not a measurement of a jointly scheduled pipeline.
"""

import argparse
import csv
import statistics
from pathlib import Path

from compare import write_tsv


STACKS = {
    "targeted_reference": (
        ("KofamScan", "subset", "all"),
        ("run_dbcan", "subset", "all"),
        ("HMMER3", "subset", "pfam"),
        ("HMMER3", "subset", "ncbifam"),
    ),
    "full_reference_interpro5": (
        ("KofamScan", "full", "all"),
        ("run_dbcan", "full", "all"),
        ("InterProScan5", "full", "all"),
    ),
    "full_reference_interpro6": (
        ("KofamScan", "full", "all"),
        ("run_dbcan", "full", "all"),
        ("InterProScan6", "full", "all"),
    ),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args()
    with (args.results / "resources.tsv").open(newline="") as handle:
        source = list(csv.DictReader(handle, delimiter="\t"))
    indexed = {}
    genomes = set()
    for row in source:
        key = (row["genome_id"], row["tool"], row["scope"], row["source"])
        if key in indexed:
            raise ValueError(f"duplicate resource measurement: {key}")
        indexed[key] = row
        genomes.add(row["genome_id"])
    details = []
    for genome in sorted(genomes):
        giftag = indexed[(genome, "giftag", "subset", "all")]
        for stack, components in STACKS.items():
            reference = [indexed[(genome, *component)] for component in components]
            giftag_wall = float(giftag["wall_seconds"])
            giftag_cpu = float(giftag["cpu_seconds"])
            reference_wall = sum(float(row["wall_seconds"]) for row in reference)
            reference_cpu = sum(float(row["cpu_seconds"]) for row in reference)
            details.append(dict(
                genome_id=genome, proteins=int(giftag["proteins"]), stack=stack,
                giftag_wall_seconds=giftag_wall,
                reference_serial_wall_seconds=reference_wall,
                wall_ratio_reference_to_giftag=reference_wall / giftag_wall,
                giftag_cpu_seconds=giftag_cpu, reference_cpu_seconds=reference_cpu,
                cpu_ratio_reference_to_giftag=(reference_cpu / giftag_cpu
                                               if giftag_cpu else "")))
    detail_cols = ("genome_id", "proteins", "stack", "giftag_wall_seconds",
                   "reference_serial_wall_seconds", "wall_ratio_reference_to_giftag",
                   "giftag_cpu_seconds", "reference_cpu_seconds",
                   "cpu_ratio_reference_to_giftag")
    write_tsv(args.results / "pipeline-effort.tsv", detail_cols, details)
    summaries = []
    for stack in STACKS:
        rows = [row for row in details if row["stack"] == stack]
        giftag_wall = sum(row["giftag_wall_seconds"] for row in rows)
        reference_wall = sum(row["reference_serial_wall_seconds"] for row in rows)
        giftag_cpu = sum(row["giftag_cpu_seconds"] for row in rows)
        reference_cpu = sum(row["reference_cpu_seconds"] for row in rows)
        summaries.append(dict(
            stack=stack, genomes=len(rows), proteins=sum(row["proteins"] for row in rows),
            giftag_total_wall_seconds=giftag_wall,
            reference_total_serial_wall_seconds=reference_wall,
            wall_ratio_total=reference_wall / giftag_wall,
            median_per_genome_wall_ratio=statistics.median(
                row["wall_ratio_reference_to_giftag"] for row in rows),
            giftag_total_cpu_seconds=giftag_cpu,
            reference_total_cpu_seconds=reference_cpu,
            cpu_ratio_total=(reference_cpu / giftag_cpu if giftag_cpu else "")))
    write_tsv(args.results / "pipeline-effort-summary.tsv",
              ("stack", "genomes", "proteins", "giftag_total_wall_seconds",
               "reference_total_serial_wall_seconds", "wall_ratio_total",
               "median_per_genome_wall_ratio", "giftag_total_cpu_seconds",
               "reference_total_cpu_seconds", "cpu_ratio_total"), summaries)


if __name__ == "__main__":
    main()
