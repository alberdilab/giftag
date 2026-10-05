#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-dbcan-inspect
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
genome=${BENCH_GENOME:?}
gene=${BENCH_GENE:?}
cluster=${BENCH_CLUSTER:?}
variant=${BENCH_VARIANT:-giftag}
out="$project/results/adjudication"
mkdir -p "$out"
if [[ "$variant" == dbcan ]]; then
    python="$project/environments/dbcan/bin/python"
    export PYTHONPATH="$project/source/src"
else
    python="$project/environments/giftag/bin/python"
fi
"$python" \
    "$project/source/benchmarks/inspect_dbcan_hit.py" \
    --db "$project/data/giftag-db" --fasta "$project/data/panel/$genome.faa" \
    --gene "$gene" --cluster "$cluster" > "$out/$genome-$gene-$cluster-$variant.json"
