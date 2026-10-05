#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-interpro-models
#SBATCH --partition=filetransfer
#SBATCH --qos=filetransfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
five=/opt/shared_software/shared_envmodules/conda/interproscan-5.75-106.0/interproscan-5.75-106.0/data
six="$project/data/interproscan"
"$project/environments/giftag/bin/python" \
    "$project/source/benchmarks/interpro_model_presence.py" \
    --markers "$project/data/giftag-db/markers.tsv" \
    --giftag-hmm "$project/data/giftag-db/ncbifam.hmm" \
    --interpro5-hmm "$five/ncbifam/17.0/ncbifam.hmm" \
    --interpro6-hmm "$six/ncbifam/18.0/ncbifam.hmm" \
    --output "$project/provenance/interpro-ncbifam-model-presence.tsv"
