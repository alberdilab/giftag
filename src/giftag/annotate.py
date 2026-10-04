"""`giftag annotate`: genomes or proteins in, a gifter marker table out."""

import csv
import json
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from giftag import GiftagError, __version__
from giftag.database import Database
from giftag.fetch import log
from giftag.genes import call_genes, genome_id, looks_nucleotide, read_fasta, write_fasta

MARKER_COLUMNS = (
    "genome_id", "gene_id", "namespace", "accession", "source", "profile", "rule",
    "score", "evalue", "threshold", "coverage", "target_from", "target_to",
)
GENOME_COLUMNS = ("genome_id", "input", "input_type", "gene_calling", "sequences",
                  "length_bp", "proteins", "marker_genes", "markers")


def annotate(inputs, outdir, db_dir, threads=0, mode="auto", input_type="auto"):
    started = time.monotonic()
    database = Database(db_dir)
    inputs = [Path(p) for p in inputs]
    if not inputs:
        raise GiftagError("no input files")
    ids = [genome_id(p) for p in inputs]
    duplicated = sorted(gid for gid, n in Counter(ids).items() if n > 1)
    if duplicated:
        raise GiftagError(f"inputs share genome IDs: {', '.join(duplicated)}")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    unsearchable = database.unsearchable
    if unsearchable:
        log(f"{len(unsearchable)} gifter markers cannot be searched with this database and will "
            f"read as absent in gifter; see markers.tsv in {database.path}")

    # Tables are written genome by genome, so a long run that stops part way
    # keeps what it finished. giftag_run.json is written last and marks a run
    # as complete.
    (outdir / "giftag_run.json").unlink(missing_ok=True)
    n_rows = 0
    with _Table(outdir / "giftag_markers.tsv", MARKER_COLUMNS) as markers, \
            _Table(outdir / "giftag_genomes.tsv", GENOME_COLUMNS) as genomes:
        for index, (path, gid) in enumerate(zip(inputs, ids), start=1):
            records = read_fasta(path)
            kind = input_type
            if kind == "auto":
                kind = "nucleotide" if looks_nucleotide(records) else "protein"
            if kind == "nucleotide":
                proteins, used = call_genes(records, mode)
                (outdir / "proteins").mkdir(exist_ok=True)
                write_fasta(outdir / "proteins" / f"{gid}.faa", proteins)
                calling = f"pyrodigal {used}"
            else:
                proteins, calling = records, "none (protein input)"

            rows = _label(gid, database.search(proteins, cpus=threads), database.labels)
            markers.write(rows)
            n_rows += len(rows)
            n_markers = len({(r["namespace"], r["accession"]) for r in rows})
            genomes.write([dict(
                genome_id=gid, input=str(path), input_type=kind, gene_calling=calling,
                sequences=len(records),
                length_bp=sum(len(s) for _, s in records) if kind == "nucleotide" else "",
                proteins=len(proteins), marker_genes=len({r["gene_id"] for r in rows}),
                markers=n_markers,
            )])
            log(f"[{index}/{len(inputs)}] {gid}: {len(proteins)} proteins, {n_markers} distinct markers")

    manifest = database.manifest
    run = {
        "giftag_version": __version__,
        "finished_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "seconds": round(time.monotonic() - started, 1),
        "database": {
            "path": str(database.path.resolve()),
            "built_utc": manifest["built_utc"],
            "gifter": manifest["gifter"],
            "sources": {name: {k: v for k, v in meta.items() if k != "terms"}
                        for name, meta in manifest["sources"].items()},
            "marker_status": manifest["marker_status"],
        },
        "parameters": {"mode": mode, "input_type": input_type, "threads": threads},
        "genomes": len(inputs),
        "rows": n_rows,
    }
    with open(outdir / "giftag_run.json", "w") as handle:
        json.dump(run, handle, indent=2)
        handle.write("\n")
    log(f"wrote {n_rows} marker rows for {len(inputs)} genomes to {outdir}")
    return n_rows


def _label(gid, calls, labels):
    """Translate accepted calls into gifter's namespace/accession, one row per
    gene and marker, keeping the strongest supporting hit."""
    best = {}
    for call in calls:
        for namespace, accession in labels.get((call.source, call.profile), ()):
            key = (call.gene_id, namespace, accession)
            if key not in best or call.score > best[key].score:
                best[key] = call
    rows = []
    for (gene, namespace, accession), call in best.items():
        rows.append(dict(
            genome_id=gid, gene_id=gene, namespace=namespace, accession=accession,
            source=call.source, profile=call.profile, rule=call.rule,
            score=f"{call.score:.1f}", evalue=f"{call.evalue:.3g}",
            threshold=f"{call.threshold:g}",
            coverage="" if call.coverage is None else f"{call.coverage:.3f}",
            target_from=call.target_from or "", target_to=call.target_to or "",
        ))
    rows.sort(key=lambda r: (r["gene_id"], r["namespace"], r["accession"]))
    return rows


class _Table:
    def __init__(self, path, columns):
        self._handle = open(path, "w", newline="")
        self._writer = csv.DictWriter(self._handle, columns, delimiter="\t", lineterminator="\n")
        self._writer.writeheader()

    def write(self, rows):
        self._writer.writerows(rows)
        self._handle.flush()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self._handle.close()
