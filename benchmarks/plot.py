"""Draw the publication benchmark's concordance and runtime figure."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


COLORS = {
    "giftag subset": "#1f3f59",
    "KofamScan full": "#2f6b9a", "KofamScan subset": "#6c9eb8",
    "run_dbcan full": "#b35d42", "run_dbcan subset": "#d38b6e",
    "HMMER3 pfam": "#93a46b", "HMMER3 ncbifam": "#637a47",
    "InterProScan5 full": "#81639a", "InterProScan6 full": "#a886c0",
}


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def group_key(row):
    if row["tool"] == "HMMER3":
        return f"HMMER3 {row['source']}"
    return f"{row['tool']} {row['scope']}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    agreement = read_tsv(args.results / "agreement.tsv")
    resources = read_tsv(args.results / "resources.tsv")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12.2, 6.3),
                                     gridspec_kw={"width_ratios": [1.2, 1]})
    labels = []
    for y, row in enumerate(agreement):
        key = group_key(row)
        label = f"{row['tool']} ({row['scope']}), {row['source'].replace('_', ' ')}"
        labels.append(label)
        if row["jaccard"]:
            value = float(row["jaccard"])
            low = float(row["jaccard_ci_low"] or value)
            high = float(row["jaccard_ci_high"] or value)
            ax_a.hlines(y, low, high, color=COLORS[key], linewidth=1.5)
            ax_a.vlines((low, high), y - 0.10, y + 0.10,
                        color=COLORS[key], linewidth=1.2)
            ax_a.plot(value, y, "o", color=COLORS[key], markersize=5)
    ax_a.set_yticks(range(len(labels)), labels)
    lower = min(float(row["jaccard_ci_low"] or row["jaccard"])
                for row in agreement if row["jaccard"])
    ax_a.set_xlim(max(0, min(0.85, lower - 0.02)), 1.01)
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Gene–marker Jaccard agreement (95% genome bootstrap CI)")
    ax_a.grid(axis="x", alpha=0.2)
    ax_a.set_title("A  Marker call agreement", loc="left", fontweight="bold")

    groups = defaultdict(list)
    for row in resources:
        key = group_key(row)
        groups[key].append((int(row["proteins"]), float(row["wall_seconds"])))
    for key, points in groups.items():
        ax_b.scatter([p[0] for p in points], [p[1] for p in points],
                     color=COLORS[key], label=key, s=24, alpha=0.85)
    ax_b.set_xlabel("Proteins in input genome")
    ax_b.set_ylabel("Wall time per genome (s, log scale)")
    ax_b.set_yscale("log")
    ax_b.grid(alpha=0.2)
    ax_b.set_title("B  Search time at eight CPUs", loc="left", fontweight="bold")
    ax_b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),
                ncol=3, fontsize=7, frameon=False)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
