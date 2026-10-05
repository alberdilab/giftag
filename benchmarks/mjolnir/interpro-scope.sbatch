#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-interpro-scope
#SBATCH --partition=filetransfer
#SBATCH --qos=filetransfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
five=/opt/shared_software/shared_envmodules/conda/interproscan-5.75-106.0/interproscan-5.75-106.0/data
six="$project/data/interproscan"
out="$project/provenance/interpro-member-scope.tsv"
printf 'tool\tinterpro_release\tcollection\tmember_release\tprofiles\tbytes\tsha256\tfile\n' > "$out"
add() {
    local tool=$1 release=$2 collection=$3 member=$4 file=$5
    local count size sha
    count=$(grep -c '^NAME ' "$file")
    size=$(stat -c%s "$file")
    sha=$(sha256sum "$file" | cut -d' ' -f1)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$tool" "$release" "$collection" "$member" "$count" "$size" "$sha" "$file" >> "$out"
}
add InterProScan5 106.0 Pfam 37.4 "$five/pfam/37.4/pfam_a.hmm"
add InterProScan5 106.0 NCBIFAM 17.0 "$five/ncbifam/17.0/ncbifam.hmm"
add InterProScan6 108.0 Pfam 38.1 "$six/pfam/38.1/pfam_a.hmm"
add InterProScan6 108.0 NCBIFAM 18.0 "$six/ncbifam/18.0/ncbifam.hmm"
