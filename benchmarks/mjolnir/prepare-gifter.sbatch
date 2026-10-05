#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-gifter-source
#SBATCH --partition=filetransfer
#SBATCH --qos=filetransfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
commit=ce63c15200b64a6614985aedd118da4d92013ea1
root="$project/software/gifter-ce63c152/gifter"
archive="$project/data/upstream/gifter/$commit.tar.gz"
mkdir -p "$(dirname "$archive")" "$root" "$project/provenance"
if [[ ! -s "$archive" ]]; then
    curl --fail --location --retry 6 --retry-delay 10 \
        --output "$archive.part" \
        "https://github.com/alberdilab/gifter/archive/$commit.tar.gz"
    mv "$archive.part" "$archive"
fi
printf '%s  %s\n' \
    8891549c49cfa25535bc4b8853e490d7d4e6ea2cc9ba28f5f6e3007ef01fe345 "$archive" \
    | sha256sum --check
sha256sum "$archive" > "$project/provenance/gifter-source-archive-sha256.txt"
if [[ ! -f "$root/DESCRIPTION" ]]; then
    tar -xzf "$archive" --strip-components=1 -C "$root"
fi

markers="$root/inst/extdata/database-source/markers.tsv"
sqlite="$root/inst/extdata/gifter.sqlite"
printf '%s  %s\n' \
    552fc134bb989ed37e0223ce05a17e0632a217416a197fbd6e36ce1d536667b9 "$markers" \
    752ff2b0cdce56e33bb12d57ce875fb4ddb459a7fd86e88c356ac1714186a87e "$sqlite" \
    | sha256sum --check
cp "$markers" "$project/data/gifter-markers.tsv"
cp "$sqlite" "$project/data/gifter.sqlite"
sha256sum "$project/data/gifter-markers.tsv" "$project/data/gifter.sqlite" \
    > "$project/provenance/gifter-input-sha256.txt"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$project/provenance/gifter-source-complete.txt"
