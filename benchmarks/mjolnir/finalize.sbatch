#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-finalize
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
source_dir="$project/source/benchmarks"
final="$project/results/final"
[[ -s "$project/provenance/gifter-env-complete.txt" ]]
[[ -s "$project/provenance/build-complete.txt" ]]
printf '%s  %s\n' \
    752ff2b0cdce56e33bb12d57ce875fb4ddb459a7fd86e88c356ac1714186a87e \
    "$project/data/gifter.sqlite" | sha256sum --check

IFS=, read -ra ips_versions <<< "${BENCH_INTERPROSCAN_VERSIONS:-6}"
"$project/environments/giftag/bin/python" "$source_dir/aggregate.py" \
    --project "$project" --outdir "$final" \
    --interproscan-version "${ips_versions[@]}"
"$project/environments/giftag/bin/python" "$source_dir/scope.py" \
    --project "$project" --output "$final/database-scope.tsv"
"$project/environments/giftag/bin/python" "$source_dir/effort.py" \
    --results "$final"
"$project/environments/gifter-r/bin/Rscript" "$source_dir/gift_impact.R" \
    "$final" "$project/data/gifter.sqlite"
"$project/environments/dbcan/bin/python" "$source_dir/plot.py" \
    --results "$final" --output "$final/benchmark.pdf"
"$project/environments/dbcan/bin/python" "$source_dir/plot.py" \
    --results "$final" --output "$final/benchmark.svg"
{
    for path in "$project/data/interproscan" "$project/cache/apptainer" \
                "$project/data/upstream/kofam/profiles"; do
        [[ ! -e "$path" ]] || du -sb "$path"
    done
} > "$project/provenance/unpacked-database-bytes.tsv"
job_ids=$(tail -n +2 "$final/resources.tsv" | cut -f9 | sort -u | paste -sd, -)
printf '%s\n' "$job_ids" > "$project/provenance/accounting-job-ids.txt"
sha256sum "$final"/*.tsv "$final"/*.json* "$final"/*.txt \
    "$final"/benchmark.pdf "$final"/benchmark.svg \
    > "$project/provenance/final-sha256.txt"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$project/provenance/final-complete.txt"
