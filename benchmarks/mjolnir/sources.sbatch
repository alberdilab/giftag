#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-sources
#SBATCH --partition=filetransfer
#SBATCH --qos=filetransfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
root="$project/data/upstream"
mkdir -p "$root/kofam" "$root/dbcan" "$project/provenance"

fetch() {
    local url=$1 destination=$2
    if [[ ! -s "$destination" ]]; then
        curl --fail --location --retry 8 --retry-delay 15 --continue-at - \
            --output "$destination.part" "$url"
        mv "$destination.part" "$destination"
    fi
}

kofam_base=https://www.genome.jp/ftp/db/kofam/archives/2026-08-02
fetch "$kofam_base/ko_list.gz" "$root/kofam/ko_list.gz"
if [[ ! -s "$root/kofam/profiles.tar.gz" ]]; then
    # GenomeNet's FTP mirror is much faster for this large archive.
    fetch ftp://ftp.genome.jp/pub/db/kofam/archives/2026-08-02/profiles.tar.gz \
        "$root/kofam/profiles.tar.gz" || \
        fetch "$kofam_base/profiles.tar.gz" "$root/kofam/profiles.tar.gz"
fi
gzip -t "$root/kofam/ko_list.gz"
tar -tzf "$root/kofam/profiles.tar.gz" > "$project/provenance/kofam-archive-members.txt"

dbcan_base=https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_5-5-2026
fetch "$dbcan_base/dbCAN.hmm" "$root/dbcan/dbCAN.hmm"
fetch "$dbcan_base/dbCAN_sub.hmm" "$root/dbcan/dbCAN_sub.hmm"
fetch "$dbcan_base/fam-substrate-mapping.tsv" "$root/dbcan/fam-substrate-mapping.tsv"

sha256sum "$root/kofam/ko_list.gz" "$root/kofam/profiles.tar.gz" \
          "$root/dbcan/dbCAN.hmm" "$root/dbcan/dbCAN_sub.hmm" \
          "$root/dbcan/fam-substrate-mapping.tsv" \
          > "$project/provenance/upstream-sha256.txt"
du -b "$root/kofam/ko_list.gz" "$root/kofam/profiles.tar.gz" \
      "$root/dbcan/dbCAN.hmm" "$root/dbcan/dbCAN_sub.hmm" \
      > "$project/provenance/upstream-bytes.tsv"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$project/provenance/sources-complete.txt"
