#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-build
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
source_dir="$project/source"
data_dir="$project/data"
env_dir="$project/environments/giftag"
[[ -s "$project/provenance/setup-complete.txt" ]]
[[ -s "$project/provenance/sources-complete.txt" ]]
mkdir -p "$project/results/build" "$project/provenance"

"$env_dir/bin/python" "$source_dir/benchmarks/measure.py" \
    --outdir "$project/results/build/giftag-subset" \
    --tool giftag-build --scope subset --threads 8 --version 0.2.0 \
    --input "$project/data/gifter-markers.tsv" -- \
    "$env_dir/bin/python" -m giftag build \
    --db "$data_dir/giftag-db" \
    --gifter-db "$project/data/gifter-markers.tsv" \
    --kofam-release 2026-08-02 --kofam-dir "$data_dir/upstream/kofam" \
    --ncbifam-release 20.0 --dbcan-release db_v5-2-9_5-5-2026 \
    --dbcan-dir "$data_dir/upstream/dbcan" --threads 8

sha256sum "$data_dir/giftag-db"/{giftag.json,markers.tsv,profiles.tsv,kofam.hmm,ncbifam.hmm,pfam.hmm,dbcan.hmm,dbcan_sub.hmm} \
    > "$project/provenance/giftag-db-sha256.txt"
du -b "$data_dir/giftag-db"/* > "$project/provenance/giftag-db-bytes.tsv"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$project/provenance/build-complete.txt"
