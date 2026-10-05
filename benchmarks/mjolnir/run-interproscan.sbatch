#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-interpro
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00

set -euo pipefail
source /etc/profile.d/modules.sh
module use /opt/shared_software/shared_envmodules/modules
module load openjdk/20.0.0 nextflow/25.10.4 apptainer/1.5.3

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
source_dir="$project/source"
panel="$project/data/panel"
python="$project/environments/giftag/bin/python"
[[ -s "$project/provenance/build-complete.txt" ]]
export NXF_HOME="$project/cache/nextflow"
export NXF_VER=26.04.6
export NXF_BIN="$project/software/nextflow-26.04.6-one.jar"
export NXF_SYNTAX_PARSER=v2
export NXF_SINGULARITY_CACHEDIR="$project/cache/apptainer"
export APPTAINER_CACHEDIR="$project/cache/apptainer"
workflow="$project/software/interproscan6-6.0.2.2"
run_root=${BENCH_RESULT_ROOT:-$project/results/runs/interproscan}
[[ -s "$project/provenance/interproscan6-source-complete.txt" ]]
[[ -s "$NXF_BIN" ]]
mkdir -p "$NXF_HOME" "$NXF_SINGULARITY_CACHEDIR" "$project/data/interproscan"

while IFS=$'\t' read -r genome assembly organism reason; do
    [[ "$genome" == genome_id ]] && continue
    [[ -z "${BENCH_GENOME:-}" || "$genome" == "$BENCH_GENOME" ]] || continue
    fasta="$panel/$genome.faa"
    out="$run_root/$genome"
    if [[ -s "$out/measurement/measure.json" &&
          -s "$out/normalized.tsv" ]]; then
        echo "already complete: InterProScan6 $genome"
        continue
    fi
    mkdir -p "$out/annotation"
    "$python" "$source_dir/benchmarks/measure.py" \
        --outdir "$out/measurement" --tool InterProScan6 --scope full \
        --threads 8 --version 6.0.2.2+Nextflow26.04.6 --input "$fasta" -- \
        nextflow run "$workflow" \
        -profile singularity -w "$out/work" \
        --input "$fasta" --datadir "$project/data/interproscan" \
        --interpro 108.0 --applications Pfam,NCBIFAM \
        --formats tsv --no-matches-api --cpus 1 --max-workers 8 \
        --outdir "$out/annotation" --outprefix "$genome"
    "$python" "$source_dir/benchmarks/compare.py" normalize \
        --format interproscan --input "$out/annotation/$genome.tsv" \
        --genome "$genome" --markers "$project/data/giftag-db/markers.tsv" \
        --output "$out/normalized.tsv"
done < "$source_dir/benchmarks/panel.tsv"
