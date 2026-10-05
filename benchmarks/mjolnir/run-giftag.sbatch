#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-annotate
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=08:00:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
source_dir="$project/source"
panel="$project/data/panel"
db="$project/data/giftag-db"
python="$project/environments/giftag/bin/python"
[[ -s "$project/provenance/build-complete.txt" ]]

while IFS=$'\t' read -r genome assembly organism reason; do
    [[ "$genome" == genome_id ]] && continue
    [[ -z "${BENCH_GENOME:-}" || "$genome" == "$BENCH_GENOME" ]] || continue
    fasta="$panel/$genome.faa"
    out="$project/results/runs/giftag/$genome"
    if [[ -s "$out/measurement/measure.json" &&
          -s "$out/annotation/giftag_run.json" ]]; then
        echo "already complete: giftag $genome"
        continue
    fi
    mkdir -p "$out"
    "$python" "$source_dir/benchmarks/measure.py" \
        --outdir "$out/measurement" --tool giftag --scope subset --threads 8 \
        --version 0.2.0 --input "$fasta" --database "$db/giftag.json" -- \
        "$python" -m giftag annotate --input "$fasta" --outdir "$out/annotation" \
        --db "$db" --input-type protein --threads 8
done < "$source_dir/benchmarks/panel.tsv"
