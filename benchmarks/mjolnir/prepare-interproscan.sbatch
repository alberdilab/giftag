#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-interpro-prepare
#SBATCH --partition=filetransfer
#SBATCH --qos=filetransfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
software="$project/software"
archive="$software/interproscan6-6.0.2.2.tar.gz"
source_dir="$software/interproscan6-6.0.2.2"
nextflow_jar="$software/nextflow-26.04.6-one.jar"
mkdir -p "$software" "$project/provenance"
if [[ ! -f "$archive" ]]; then
    curl --fail --location --retry 6 \
        --output "$archive.part" \
        https://github.com/ebi-pf-team/interproscan6/archive/refs/tags/6.0.2.2.tar.gz
    mv "$archive.part" "$archive"
fi
sha256sum "$archive" > "$project/provenance/interproscan6-source-sha256.txt"
if [[ ! -f "$nextflow_jar" ]]; then
    curl --fail --location --retry 6 \
        --output "$nextflow_jar.part" \
        https://www.nextflow.io/releases/v26.04.6/nextflow-26.04.6-one.jar
    mv "$nextflow_jar.part" "$nextflow_jar"
fi
sha256sum "$nextflow_jar" > "$project/provenance/nextflow-26.04.6-sha256.txt"
if [[ ! -f "$source_dir/main.nf" ]]; then
    tar -xzf "$archive" -C "$software"
fi
[[ -f "$source_dir/main.nf" ]]
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$project/provenance/interproscan6-source-complete.txt"
