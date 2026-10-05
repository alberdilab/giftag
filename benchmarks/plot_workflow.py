"""Draw the documentation figures: concordance and whole-workflow runtime.

Each chart is written to its own file. Unlike plot.py, the runtime chart
compares giftag with the summed reference stages needed to obtain the same
gifter markers, from pipeline-effort.tsv.
"""

import argparse
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot import read_tsv


GIFTAG = "#0d3d38"
REFERENCE = "#b35d42"
MUTED = "#687672"
TICKS = (30, 100, 300, 1000, 3000)
STACKS = {
    "giftag": "giftag\none run, all four collections",
    "targeted_reference": "Original tools, giftag's profiles\n"
                          "KofamScan + run_dbcan + HMMER",
    "full_reference_interpro6": "Original tools, full libraries\n"
                                "KofamScan + run_dbcan + InterProScan 6",
    "full_reference_interpro5": "Original tools, full libraries\n"
                                "KofamScan + run_dbcan + InterProScan 5",
}


def seconds(value):
    return f"{value:,.0f} s"


def save(fig, ax, path):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", transparent=True)
    plt.close(fig)


def plot_agreement(results, path):
    agreement = read_tsv(results / "agreement.tsv")
    fig, ax_a = plt.subplots(figsize=(8.6, 4.8))
    labels = []
    for y, row in enumerate(agreement):
        labels.append(f"{row['tool']} ({row['scope']}), "
                      f"{row['source'].replace('_', ' ')}")
        value = float(row["jaccard"])
        low = float(row["jaccard_ci_low"] or value)
        high = float(row["jaccard_ci_high"] or value)
        ax_a.hlines(y, low, high, color=GIFTAG, linewidth=1.5)
        ax_a.vlines((low, high), y - 0.10, y + 0.10, color=GIFTAG, linewidth=1.2)
        ax_a.plot(value, y, "o", color=GIFTAG, markersize=5)
    ax_a.set_yticks(range(len(labels)), labels)
    ax_a.set_xlim(0.85, 1.01)
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Gene–marker Jaccard agreement (95% genome bootstrap CI)")
    ax_a.grid(axis="x", alpha=0.2)
    save(fig, ax_a, path)


def plot_workflow(results, path):
    effort = read_tsv(results / "pipeline-effort.tsv")
    summary = {row["stack"]: row
               for row in read_tsv(results / "pipeline-effort-summary.tsv")}
    fig, ax_b = plt.subplots(figsize=(8.6, 4.4))
    walls = defaultdict(dict)
    for row in effort:
        walls["giftag"][row["genome_id"]] = float(row["giftag_wall_seconds"])
        walls[row["stack"]][row["genome_id"]] = float(
            row["reference_serial_wall_seconds"])
    for y, stack in enumerate(STACKS):
        values = list(walls[stack].values())
        color = GIFTAG if stack == "giftag" else REFERENCE
        median = statistics.median(values)
        ax_b.scatter(values, [y] * len(values), color=color, s=34, alpha=0.6,
                     linewidths=0, zorder=2)
        ax_b.vlines(median, y - 0.17, y + 0.17, color=color, linewidth=2.4,
                    zorder=3)
        note = f"median {seconds(median)}"
        if stack != "giftag":
            ratio = float(summary[stack]["wall_ratio_total"])
            note += f"  ·  {ratio:.1f}× giftag's total"
        left = median < 1000
        ax_b.annotate(note, (min(values) if left else max(values), y),
                      xytext=(-4 if left else 4, 20),
                      textcoords="offset points",
                      ha="left" if left else "right", fontsize=9, color=MUTED)
    ax_b.set_yticks(range(len(STACKS)), list(STACKS.values()))
    ax_b.set_ylim(len(STACKS) - 0.5, -0.7)
    ax_b.set_xscale("log")
    ax_b.set_xlim(20, 6000)
    ax_b.set_xticks(TICKS, [f"{tick:,}" for tick in TICKS])
    ax_b.minorticks_off()
    ax_b.set_xlabel("Wall time per genome at eight CPUs (s, log scale)")
    ax_b.grid(axis="x", alpha=0.2)
    save(fig, ax_b, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--format", default="svg")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_agreement(args.results,
                   args.output_dir / f"benchmark-agreement.{args.format}")
    plot_workflow(args.results,
                  args.output_dir / f"benchmark-workflow.{args.format}")


if __name__ == "__main__":
    main()
