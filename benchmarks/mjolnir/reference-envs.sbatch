#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-ref-envs
#SBATCH --partition=filetransfer
#SBATCH --qos=filetransfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=1-00:00:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
envs="$project/environments"
software="$project/software"
mkdir -p "$envs" "$software" "$project/provenance"

if [[ ! -x "$envs/kofam/bin/ruby" ]]; then
    /home/jpl786/miniforge3/bin/mamba create -y -p "$envs/kofam" \
        -c conda-forge -c bioconda 'ruby=3.3' 'parallel' 'hmmer=3.4'
fi
if [[ ! -x "$software/kofam_scan/exec_annotation" ]]; then
    curl --fail --location --retry 6 \
        --output "$software/kofam_scan-v1.3.0.tar.gz.part" \
        https://github.com/takaram/kofam_scan/archive/refs/tags/v1.3.0.tar.gz
    mv "$software/kofam_scan-v1.3.0.tar.gz.part" "$software/kofam_scan-v1.3.0.tar.gz"
    tar -xzf "$software/kofam_scan-v1.3.0.tar.gz" -C "$software"
    mv "$software/kofam_scan-1.3.0" "$software/kofam_scan"
    chmod +x "$software/kofam_scan/exec_annotation"
fi
sha256sum "$software/kofam_scan-v1.3.0.tar.gz" \
    > "$project/provenance/kofamscan-source-sha256.txt"
/home/jpl786/miniforge3/bin/conda list -p "$envs/kofam" --explicit \
    > "$project/provenance/kofamscan-conda-explicit.txt"
PATH="$envs/kofam/bin:$PATH" "$software/kofam_scan/exec_annotation" --version \
    > "$project/provenance/kofamscan-version.txt"

if [[ ! -x "$envs/dbcan/bin/python" ]]; then
    /home/jpl786/miniforge3/bin/mamba create -y -p "$envs/dbcan" \
        -c conda-forge -c bioconda 'python=3.11' 'pip'
fi
"$envs/dbcan/bin/python" -m pip install 'dbcan==5.2.9'
"$envs/dbcan/bin/python" -m pip freeze > "$project/provenance/dbcan-pip-freeze.txt"
/home/jpl786/miniforge3/bin/conda list -p "$envs/dbcan" --explicit \
    > "$project/provenance/dbcan-conda-explicit.txt"
"$envs/dbcan/bin/run_dbcan" version > "$project/provenance/run-dbcan-version.txt"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$project/provenance/reference-envs-complete.txt"
