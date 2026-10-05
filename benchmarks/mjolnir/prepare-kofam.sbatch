#!/usr/bin/env bash
#SBATCH --job-name=giftag-bench-kofam-prepare
#SBATCH --partition=cpuqueue
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

project=/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026
root="$project/data/upstream/kofam"
[[ -s "$project/provenance/build-complete.txt" ]]
if [[ ! -d "$root/profiles" ]]; then
    tar -xzf "$root/profiles.tar.gz" -C "$root"
fi
"$project/environments/giftag/bin/python" - "$project/data/giftag-db/profiles.tsv" "$root" <<'PY'
import csv
import sys
from pathlib import Path

table, root = Path(sys.argv[1]), Path(sys.argv[2])
with table.open(newline="") as handle:
    profiles = sorted({r["profile"] for r in csv.DictReader(handle, delimiter="\t")
                       if r["source"] == "kofam"})
if not profiles:
    raise SystemExit("no KOfam profiles in giftag database")
paths = [(root / "profiles" / f"{profile}.hmm").resolve() for profile in profiles]
missing = [str(path) for path in paths if not path.is_file()]
if missing:
    raise SystemExit(f"{len(missing)} KOfam profiles missing; first: {missing[0]}")
(root / "gifter-subset.hal").write_text("".join(f"{path}\n" for path in paths))
print(f"KofamScan subset: {len(paths)} profiles")
PY
sha256sum "$root/gifter-subset.hal" \
    > "$project/provenance/kofamscan-subset-sha256.txt"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$project/provenance/kofam-prepared.txt"
