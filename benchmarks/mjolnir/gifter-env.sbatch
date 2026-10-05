#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-gifter-env
#SBATCH --partition=filetransfer
#SBATCH --qos=filetransfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
env="$project/environments/gifter-r"
source_dir="$project/software/gifter-ce63c152/gifter"
mkdir -p "$project/provenance"
[[ -s "$project/provenance/gifter-source-complete.txt" ]]
if [[ ! -x "$env/bin/R" ]]; then
    /home/jpl786/miniforge3/bin/mamba create -y -p "$env" -c conda-forge \
        'r-base=4.3' 'r-dbi' 'r-rsqlite' 'r-cli' 'r-tibble'
fi
"$env/bin/R" CMD INSTALL "$source_dir"
"$env/bin/Rscript" -e 'stopifnot(as.character(packageVersion("gifter")) == "0.7.3"); sessionInfo()' \
    > "$project/provenance/gifter-r-session-info.txt"
/home/jpl786/miniforge3/bin/conda list -p "$env" --explicit \
    > "$project/provenance/gifter-r-conda-explicit.txt"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$project/provenance/gifter-env-complete.txt"
