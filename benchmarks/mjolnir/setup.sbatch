#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-setup
#SBATCH --partition=filetransfer
#SBATCH --qos=filetransfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=2-00:00:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
source_dir="$project/source"
data_dir="$project/data"
env_dir="$project/environments/giftag"
mkdir -p "$data_dir" "$project/logs" "$project/provenance" "$project/environments"
export PATH="$env_dir/bin:$PATH"

if [[ ! -x "$env_dir/bin/python" ]]; then
    /home/jpl786/miniforge3/bin/mamba create -y -p "$env_dir" \
        -c conda-forge -c bioconda \
        'python=3.11' 'pip' 'pyhmmer>=0.12,<0.14' 'pyrodigal>=3.5,<4' \
        'rich>=13,<15' 'typer>=0.12,<1' 'ncbi-datasets-cli>=18,<20'
fi
"$env_dir/bin/python" -m pip install --no-deps -e "$source_dir"
"$env_dir/bin/python" -m pip freeze > "$project/provenance/giftag-pip-freeze.txt"
/home/jpl786/miniforge3/bin/conda list -p "$env_dir" --explicit \
    > "$project/provenance/giftag-conda-explicit.txt"
"$env_dir/bin/python" -m giftag version > "$project/provenance/giftag-version.txt"
"$env_dir/bin/datasets" --version > "$project/provenance/ncbi-datasets-version.txt"

gifter="$project/data/gifter-markers.tsv"
[[ -s "$project/provenance/gifter-source-complete.txt" ]]
printf '%s  %s\n' \
    552fc134bb989ed37e0223ce05a17e0632a217416a197fbd6e36ce1d536667b9 "$gifter" \
    | sha256sum --check
"$env_dir/bin/python" "$source_dir/benchmarks/fetch_panel.py" \
    --panel "$source_dir/benchmarks/panel.tsv" --outdir "$data_dir/panel"
"$env_dir/bin/python" - \
    "$source_dir/benchmarks/data/panel-manifest.tsv" \
    "$data_dir/panel/panel_manifest.tsv" <<'PY'
import csv
import sys

def load(path):
    with open(path, newline="") as handle:
        return {row["genome_id"]: row for row in csv.DictReader(handle, delimiter="\t")}

expected, observed = map(load, sys.argv[1:])
if expected.keys() != observed.keys():
    raise SystemExit("NCBI panel genome IDs changed")
for genome in expected:
    for field in ("assembly", "proteins", "benchmark_faa_sha256"):
        if expected[genome][field] != observed[genome][field]:
            raise SystemExit(f"NCBI panel changed: {genome} {field}")
print(f"verified {len(expected)} pinned protein FASTA files")
PY

sha256sum "$source_dir/benchmarks/panel.tsv" "$gifter" \
    > "$project/provenance/input-sha256.txt"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$project/provenance/setup-complete.txt"
