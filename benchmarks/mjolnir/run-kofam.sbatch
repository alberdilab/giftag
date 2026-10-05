#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-kofam
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
source_dir="$project/source"
root="$project/data/upstream/kofam"
panel="$project/data/panel"
env="$project/environments/kofam"
python="$project/environments/giftag/bin/python"
scan="$project/software/kofam_scan/exec_annotation"
[[ -s "$project/provenance/kofam-prepared.txt" ]]
[[ -s "$project/provenance/reference-envs-complete.txt" ]]
gzip -cd "$root/ko_list.gz" > "$root/ko_list"
export PATH="$env/bin:$PATH"
tmpdir=
cleanup() {
    [[ -z "$tmpdir" ]] || rm -rf -- "$tmpdir"
}
trap cleanup EXIT

while IFS=$'\t' read -r genome assembly organism reason; do
    [[ "$genome" == genome_id ]] && continue
    [[ -z "${BENCH_GENOME:-}" || "$genome" == "$BENCH_GENOME" ]] || continue
    fasta="$panel/$genome.faa"
    for scope in full subset; do
        if [[ "$scope" == full ]]; then
            profiles="$root/profiles"
        else
            profiles="$root/gifter-subset.hal"
        fi
        out="$project/results/runs/kofam/$scope/$genome"
        if [[ -s "$out/measurement/measure.json" &&
              -s "$out/normalized.tsv" ]]; then
            echo "already complete: KofamScan $scope $genome"
            continue
        fi
        mkdir -p "$out"
        tmpdir=$(mktemp -d "${SLURM_TMPDIR:-/tmp}/giftag-kofam-${genome}-${scope}-XXXXXX")
        "$python" "$source_dir/benchmarks/measure.py" \
            --outdir "$out/measurement" --tool KofamScan --scope "$scope" \
            --threads 8 --version 1.3.0 --input "$fasta" \
            --database "$root/ko_list.gz" --database "$root/gifter-subset.hal" -- \
            "$scan" -p "$profiles" -k "$root/ko_list" --cpu=8 \
            --tmp-dir "$tmpdir" -f mapper -o "$out/mapper.tsv" "$fasta"
        "$python" "$source_dir/benchmarks/compare.py" normalize \
            --format kofam-mapper --input "$out/mapper.tsv" --genome "$genome" \
            --markers "$project/data/giftag-db/markers.tsv" \
            --output "$out/normalized.tsv"
        cleanup
        tmpdir=
    done
done < "$source_dir/benchmarks/panel.tsv"
