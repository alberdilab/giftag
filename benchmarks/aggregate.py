"""Collect a complete Mjolnir panel into publication-ready TSVs.

Fails if any expected run is missing. This prevents a partial panel from being
reported as a complete benchmark. Use --genomes for an explicit pilot subset.
"""

import argparse
import csv
import json
import random
import statistics
from collections import defaultdict
from pathlib import Path

from compare import FIELDS, compare, read_tsv, searchable_markers, write_tsv


BASE_COMPARISONS = (
    ("KofamScan", "full", "kofam", "kofam/full", "normalized.tsv"),
    ("KofamScan", "subset", "kofam", "kofam/subset", "normalized.tsv"),
    ("run_dbcan", "full", "dbcan_family", "dbcan/full", "family-normalized.tsv"),
    ("run_dbcan", "full", "dbcan_sub", "dbcan/full", "sub-normalized.tsv"),
    ("run_dbcan", "subset", "dbcan_family", "dbcan/subset", "family-normalized.tsv"),
    ("run_dbcan", "subset", "dbcan_sub", "dbcan/subset", "sub-normalized.tsv"),
    ("HMMER3", "subset", "pfam", "hmmer/pfam", "normalized.tsv"),
    ("HMMER3", "subset", "ncbifam", "hmmer/ncbifam", "normalized.tsv"),
)


def comparison_specs(interproscan_versions):
    specs = list(BASE_COMPARISONS)
    for version in interproscan_versions:
        tool = f"InterProScan{version}"
        relative = "interproscan" if version == "6" else "interproscan5"
        specs.extend((tool, "full", source, relative, "normalized.tsv")
                     for source in ("pfam", "ncbifam"))
    return specs


def bootstrap_jaccard(genomes, n=10000):
    """95% percentile interval from genome-block resampling, fixed seed."""
    if len(genomes) < 2:
        return None, None, 0
    rng = random.Random(20261004)
    values = []
    for _ in range(n):
        sampled = rng.choices(genomes, k=len(genomes))
        shared = sum(r["shared"] for r in sampled)
        union = shared + sum(r["giftag_only"] + r["reference_only"] for r in sampled)
        if union:
            values.append(shared / union)
    if not values:
        return None, None, 0
    values.sort()
    return (values[int(0.025 * (len(values) - 1))],
            values[int(0.975 * (len(values) - 1))], len(values))


def collect(project, outdir, selected=None, interproscan_versions=("6",)):
    project, outdir = Path(project), Path(outdir)
    panel = list(csv.DictReader((project / "source/benchmarks/panel.tsv").open(), delimiter="\t"))
    if selected:
        unknown = set(selected) - {r["genome_id"] for r in panel}
        if unknown:
            raise ValueError(f"unknown genomes: {', '.join(sorted(unknown))}")
        panel = [r for r in panel if r["genome_id"] in selected]
    if not panel:
        raise ValueError("empty panel")
    outdir.mkdir(parents=True, exist_ok=True)
    markers = searchable_markers(project / "data/giftag-db/markers.tsv")
    giftag_rows = []
    resources = []
    measurements = []
    proteins = {r["genome_id"]: int(r["proteins"]) for r in
                read_tsv(project / "data/panel/panel_manifest.tsv")}

    def measured(path, source="all"):
        measure = path / "measurement/measure.json"
        if not measure.is_file():
            raise FileNotFoundError(measure)
        meta = json.loads(measure.read_text())
        if meta["exit_code"]:
            raise ValueError(f"failed benchmark: {measure}")
        genome = path.name
        measurements.append(dict(meta, benchmark_genome_id=genome,
                                 benchmark_source=source))
        resources.append(dict(tool=meta["tool"], scope=meta["scope"], source=source,
                              genome_id=genome, proteins=proteins[genome],
                              wall_seconds=meta["wall_seconds"],
                              cpu_seconds=meta["user_cpu_seconds"] + meta["system_cpu_seconds"],
                              peak_child_rss_bytes=meta["peak_child_rss_bytes"],
                              slurm_job_id=meta["slurm_job_id"] or ""))

    for entry in panel:
        genome = entry["genome_id"]
        path = project / "results/runs/giftag" / genome
        marker_path = path / "annotation/giftag_markers.tsv"
        if not marker_path.is_file():
            raise FileNotFoundError(marker_path)
        giftag_rows.extend(read_tsv(marker_path))
        measured(path)
    giftag_path = outdir / "giftag-combined.tsv"
    write_tsv(giftag_path, FIELDS,
              [{k: r[k] for k in FIELDS} for r in giftag_rows])

    agreement = []
    measured_paths = set()
    comparisons = comparison_specs(interproscan_versions)
    for tool, scope, source, relative, filename in comparisons:
        rows = []
        for entry in panel:
            genome = entry["genome_id"]
            path = project / "results/runs" / relative / genome
            reference_path = path / filename
            if not reference_path.is_file():
                raise FileNotFoundError(reference_path)
            rows.extend(read_tsv(reference_path))
            key = str(path)
            if key not in measured_paths:
                measured(path, source if relative.startswith("hmmer/") else "all")
                measured_paths.add(key)
        reference = outdir / f"{tool}-{scope}-{source}-combined.tsv"
        write_tsv(reference, FIELDS, rows)
        summary, discordance, by_genome, by_marker = compare(
            giftag_path, reference, markers, source)
        observed = {r["genome_id"]: r for r in by_genome}
        blocks = [observed.get(entry["genome_id"],
                               {"shared": 0, "giftag_only": 0, "reference_only": 0})
                  for entry in panel]
        low, high, valid = bootstrap_jaccard(blocks)
        summary.update(tool=tool, scope=scope, genomes=len(panel),
                       jaccard_ci_low=low, jaccard_ci_high=high,
                       bootstrap_valid_resamples=valid)
        agreement.append(summary)
        stem = f"{tool}-{scope}-{source}"
        write_tsv(outdir / f"{stem}-discordant.tsv", FIELDS + ("status",), discordance)
        write_tsv(outdir / f"{stem}-by-genome.tsv",
                  ("genome_id", "shared", "giftag_only", "reference_only"), by_genome)
        write_tsv(outdir / f"{stem}-by-marker.tsv",
                  ("namespace", "accession", "shared", "giftag_only", "reference_only"),
                  by_marker)
    write_tsv(outdir / "agreement.tsv",
              ("tool", "scope", "source", "genomes", "searchable_markers", "giftag_pairs",
               "reference_pairs", "shared_pairs", "giftag_only_pairs", "reference_only_pairs",
               "jaccard", "jaccard_ci_low", "jaccard_ci_high",
               "bootstrap_valid_resamples", "reference_recall",
               "reference_precision"), agreement)
    write_tsv(outdir / "resources.tsv",
              ("tool", "scope", "source", "genome_id", "proteins", "wall_seconds",
               "cpu_seconds", "peak_child_rss_bytes", "slurm_job_id"), resources)
    with (outdir / "measurements.jsonl").open("w") as handle:
        for record in measurements:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    grouped = defaultdict(list)
    for row in resources:
        grouped[(row["tool"], row["scope"], row["source"])].append(row)
    resource_summary = []
    for (tool, scope, source), rows in sorted(grouped.items()):
        total_wall = sum(row["wall_seconds"] for row in rows)
        total_proteins = sum(row["proteins"] for row in rows)
        resource_summary.append(dict(
            tool=tool, scope=scope, source=source, genomes=len(rows),
            proteins=total_proteins,
            median_wall_seconds=statistics.median(row["wall_seconds"] for row in rows),
            total_wall_seconds=total_wall,
            total_cpu_seconds=sum(row["cpu_seconds"] for row in rows),
            max_peak_child_rss_bytes=max(row["peak_child_rss_bytes"] for row in rows),
            wall_seconds_per_1000_proteins=1000 * total_wall / total_proteins))
    write_tsv(outdir / "resource-summary.tsv",
              ("tool", "scope", "source", "genomes", "proteins",
               "median_wall_seconds", "total_wall_seconds", "total_cpu_seconds",
               "max_peak_child_rss_bytes", "wall_seconds_per_1000_proteins"),
              resource_summary)
    (outdir / "manifest.json").write_text(json.dumps(
        {"genomes": [r["genome_id"] for r in panel], "comparisons": len(agreement),
         "measurements": len(measurements),
         "interproscan_versions": list(interproscan_versions),
         "bootstrap_seed": 20261004, "bootstrap_resamples": 10000}, indent=2) + "\n")
    print(json.dumps({"genomes": len(panel), "comparisons": len(agreement)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--genomes", nargs="+", help="explicit pilot subset")
    parser.add_argument("--interproscan-version", nargs="+", choices=("5", "6"),
                        default=["6"], help="require these InterProScan versions")
    args = parser.parse_args()
    collect(args.project, args.outdir, args.genomes, args.interproscan_version)


if __name__ == "__main__":
    main()
