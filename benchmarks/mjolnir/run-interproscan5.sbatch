#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-interpro5
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=01:00:00

set -euo pipefail
source /etc/profile.d/modules.sh
module use /opt/shared_software/shared_envmodules/modules
module load openjdk/20.0.0 interproscan/5.75-106.0

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
source_dir="$project/source"
panel="$project/data/panel"
python="$project/environments/giftag/bin/python"
[[ -s "$project/provenance/build-complete.txt" ]]
interproscan=$(command -v interproscan.sh)

while IFS=$'\t' read -r genome assembly organism reason; do
    [[ "$genome" == genome_id ]] && continue
    [[ -z "${BENCH_GENOME:-}" || "$genome" == "$BENCH_GENOME" ]] || continue
    fasta="$panel/$genome.faa"
    out="$project/results/runs/interproscan5/$genome"
    if [[ -s "$out/measurement/measure.json" &&
          -s "$out/normalized.tsv" ]]; then
        echo "already complete: InterProScan5 $genome"
        continue
    fi
    mkdir -p "$out/tmp"
    "$python" "$source_dir/benchmarks/measure.py" \
        --outdir "$out/measurement" --tool InterProScan5 --scope full \
        --threads 8 --version 5.75-106.0 --input "$fasta" -- \
        "$interproscan" -i "$fasta" -f TSV -appl Pfam,NCBIFAM \
        -cpu 8 -dp -T "$out/tmp" -o "$out/annotation.tsv"
    "$python" "$source_dir/benchmarks/compare.py" normalize \
        --format interproscan --input "$out/annotation.tsv" \
        --genome "$genome" --markers "$project/data/giftag-db/markers.tsv" \
        --output "$out/normalized.tsv"
done < "$source_dir/benchmarks/panel.tsv"
