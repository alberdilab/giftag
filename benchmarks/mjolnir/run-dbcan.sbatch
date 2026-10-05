#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-dbcan
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
source_dir="$project/source"
panel="$project/data/panel"
db="$project/data/upstream/dbcan"
python="$project/environments/giftag/bin/python"
run_dbcan="$project/environments/dbcan/bin/run_dbcan"
[[ -s "$project/provenance/build-complete.txt" ]]
[[ -s "$project/provenance/reference-envs-complete.txt" ]]
[[ -e "$db/dbCAN-sub.hmm" ]] || ln -s dbCAN_sub.hmm "$db/dbCAN-sub.hmm"

while IFS=$'\t' read -r genome assembly organism reason; do
    [[ "$genome" == genome_id ]] && continue
    [[ -z "${BENCH_GENOME:-}" || "$genome" == "$BENCH_GENOME" ]] || continue
    fasta="$panel/$genome.faa"
    out="$project/results/runs/dbcan/full/$genome"
    if [[ -s "$out/measurement/measure.json" &&
          -s "$out/family-normalized.tsv" &&
          -s "$out/sub-normalized.tsv" ]]; then
        echo "already complete: run_dbcan $genome"
        continue
    fi
    mkdir -p "$out"
    "$python" "$source_dir/benchmarks/measure.py" \
        --outdir "$out/measurement" --tool run_dbcan --scope full --threads 8 \
        --version 5.2.9 --input "$fasta" -- \
        "$run_dbcan" CAZyme_annotation --mode protein \
        --input_raw_data "$fasta" --output_dir "$out/annotation" \
        --db_dir "$db" --methods hmm,dbCANsub --threads 8
    "$python" "$source_dir/benchmarks/compare.py" normalize \
        --format dbcan-family --input "$out/annotation/dbCAN_hmm_results.tsv" \
        --genome "$genome" --markers "$project/data/giftag-db/markers.tsv" \
        --output "$out/family-normalized.tsv"
    "$python" "$source_dir/benchmarks/compare.py" normalize \
        --format dbcan-sub --input "$out/annotation/dbCANsub_hmm_results.tsv" \
        --genome "$genome" --markers "$project/data/giftag-db/markers.tsv" \
        --output "$out/sub-normalized.tsv"
done < "$source_dir/benchmarks/panel.tsv"
