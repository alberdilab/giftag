"""Verify copied publication results against their Mjolnir SHA-256 manifest."""

import argparse
import hashlib
from pathlib import Path


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=root / "results/final")
    parser.add_argument("--manifest", type=Path,
                        default=root / "provenance/final-sha256.txt")
    args = parser.parse_args()
    seen = set()
    for line in args.manifest.read_text().splitlines():
        expected, remote_path = line.split("  ", 1)
        name = Path(remote_path).name
        if name in seen:
            raise ValueError(f"duplicate result filename: {name}")
        seen.add(name)
        path = args.results / name
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        if digest.hexdigest() != expected:
            raise ValueError(f"SHA-256 mismatch: {path}")
    if not seen:
        raise ValueError("empty result checksum manifest")
    print(f"verified {len(seen)} publication result files")


if __name__ == "__main__":
    main()
