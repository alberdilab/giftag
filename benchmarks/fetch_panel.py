"""Download the fixed NCBI RefSeq protein panel and assign unique protein IDs.

Requires the NCBI `datasets` CLI. Keeps the original package and header map so
all transformed sequences can be traced to the public accession and release.
"""

import argparse
import csv
import hashlib
import json
import subprocess
import zipfile
from pathlib import Path


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def records(data):
    header, sequence = None, []
    for line in data.decode("utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(sequence)
            header, sequence = line[1:], []
        elif line.strip():
            if header is None:
                raise ValueError("sequence before first FASTA header")
            sequence.append(line.strip())
    if header is not None:
        yield header, "".join(sequence)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    panel = list(csv.DictReader(args.panel.open(newline=""), delimiter="\t"))
    summary = []
    for entry in panel:
        gid, accession = entry["genome_id"], entry["assembly"]
        zip_path = args.outdir / f"{accession}.zip"
        if not zip_path.exists():
            part = zip_path.with_suffix(".zip.part")
            part.unlink(missing_ok=True)
            subprocess.run(["datasets", "download", "genome", "accession", accession,
                            "--include", "protein", "--filename", str(part)], check=True)
            with zipfile.ZipFile(part) as archive:
                if archive.testzip() is not None:
                    raise ValueError(f"{accession}: invalid NCBI data package")
            part.rename(zip_path)
        with zipfile.ZipFile(zip_path) as archive:
            path = f"ncbi_dataset/data/{accession}/protein.faa"
            if path not in archive.namelist():
                raise ValueError(f"{accession}: no protein.faa in NCBI data package")
            data = archive.read(path)
            report_path = f"ncbi_dataset/data/assembly_data_report.jsonl"
            if report_path in archive.namelist():
                (args.outdir / f"{accession}.assembly.jsonl").write_bytes(
                    archive.read(report_path))
        proteins = list(records(data))
        if not proteins or any(not sequence for _, sequence in proteins):
            raise ValueError(f"{accession}: missing or empty protein sequences")
        fasta = args.outdir / f"{gid}.faa"
        mapping = args.outdir / f"{gid}.headers.tsv"
        with open(fasta, "w") as faa, open(mapping, "w", newline="") as tsv:
            writer = csv.writer(tsv, delimiter="\t", lineterminator="\n")
            writer.writerow(("gene_id", "original_header", "length_aa"))
            for index, (header, sequence) in enumerate(proteins, 1):
                gene_id = f"{gid}_p{index:06d}"
                faa.write(f">{gene_id}\n{sequence}\n")
                writer.writerow((gene_id, header, len(sequence)))
        summary.append(dict(genome_id=gid, assembly=accession, proteins=len(proteins),
                            source_zip_sha256=sha256(zip_path.read_bytes()),
                            source_protein_sha256=sha256(data),
                            benchmark_faa_sha256=sha256(fasta.read_bytes())))
    with open(args.outdir / "panel_manifest.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, summary[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    print(json.dumps({"genomes": len(summary), "proteins": sum(r["proteins"] for r in summary)}))


if __name__ == "__main__":
    main()
