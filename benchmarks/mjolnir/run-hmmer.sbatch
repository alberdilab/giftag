#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-hmmer
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
source_dir="$project/source"
panel="$project/data/panel"
db="$project/data/giftag-db"
python="$project/environments/giftag/bin/python"
hmmsearch="$project/environments/kofam/bin/hmmsearch"
[[ -s "$project/provenance/build-complete.txt" ]]
[[ -s "$project/provenance/reference-envs-complete.txt" ]]

while IFS=$'\t' read -r genome assembly organism reason; do
    [[ "$genome" == genome_id ]] && continue
    [[ -z "${BENCH_GENOME:-}" || "$genome" == "$BENCH_GENOME" ]] || continue
    fasta="$panel/$genome.faa"
    for source in pfam ncbifam; do
        out="$project/results/runs/hmmer/$source/$genome"
        if [[ -s "$out/measurement/measure.json" &&
              -s "$out/normalized.tsv" ]]; then
            echo "already complete: HMMER3 $source $genome"
            continue
        fi
        mkdir -p "$out"
        if [[ "$source" == pfam ]]; then cutoff=--cut_ga; else cutoff=--cut_tc; fi
        "$python" "$source_dir/benchmarks/measure.py" \
            --outdir "$out/measurement" --tool HMMER3 --scope subset \
            --threads 8 --version 3.4 --input "$fasta" -- \
            "$hmmsearch" --cpu 8 "$cutoff" --domtblout "$out/hits.domtblout" \
            -o /dev/null "$db/$source.hmm" "$fasta"
        "$python" "$source_dir/benchmarks/compare.py" normalize \
            --format "hmmsearch-$source" --input "$out/hits.domtblout" \
            --genome "$genome" --markers "$db/markers.tsv" \
            --output "$out/normalized.tsv"
    done
done < "$source_dir/benchmarks/panel.tsv"
