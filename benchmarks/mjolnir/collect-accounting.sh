#!/usr/bin/env bash
# Run on Mjolnir's login node after finalize.sbatch finishes.
set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
resources="$project/results/final/resources.tsv"
output="$project/provenance/slurm-accounting.psv"
[[ -s "$resources" ]]
job_ids=$(tail -n +2 "$resources" | cut -f9 | sort -u | paste -sd, -)
[[ -n "$job_ids" ]]
sacct -j "$job_ids" --parsable2 \
    --format=JobIDRaw,JobName,State,ElapsedRaw,TotalCPU,MaxRSS,AllocCPUS,ReqMem,NodeList \
    > "$output.tmp"
[[ $(wc -l < "$output.tmp") -gt 1 ]]
mv "$output.tmp" "$output"
printf '%s\n' "$job_ids" > "$project/provenance/accounting-job-ids.txt"
printf 'Saved %s accounting rows\n' "$(($(wc -l < "$output") - 1))"
