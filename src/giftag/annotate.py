"""`giftag annotate`: genomes or proteins in, a gifter marker table out."""

import csv
import json
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from rich.markup import escape

from giftag import GiftagError, __version__, ui
from giftag.build import cazy_family
from giftag.database import Database
from giftag.genes import call_genes, genome_id, looks_nucleotide, read_fasta, write_fasta

MARKER_COLUMNS = (
    "genome_id", "gene_id", "namespace", "accession", "source", "profile", "rule",
    "score", "evalue", "threshold", "coverage", "target_from", "target_to",
)
GENOME_COLUMNS = ("genome_id", "input", "input_type", "gene_calling", "sequences",
                  "length_bp", "proteins", "marker_genes", "markers")


def annotate(inputs, outdir, db_dir, threads=0, mode="auto", input_type="auto",
             gate_subfamilies=False):
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
        ui.warning(f"{len(unsearchable)} gifter markers cannot be searched with this database "
                   f"and will read as absent in gifter (see markers.tsv in the database)")
    with ui.Task("loading profiles", done="profiles loaded in {elapsed}"):
        database.load_all()

    # Tables are written genome by genome, so a long run that stops part way
    # keeps what it finished. giftag_run.json is written last and marks a run
    # as complete.
    (outdir / "giftag_run.json").unlink(missing_ok=True)
    n_rows = 0
    summary = []
    width = len(str(len(inputs)))
    with _Table(outdir / "giftag_markers.tsv", MARKER_COLUMNS) as markers, \
            _Table(outdir / "giftag_genomes.tsv", GENOME_COLUMNS) as genomes, \
            ui.Task("annotating", total=len(inputs), steps=False) as task:
        for index, (path, gid) in enumerate(zip(inputs, ids), start=1):
            began = time.monotonic()
            name = escape(gid)

            def stage(step):
                task.update(description=f"{name} · {step}")

            stage("reading")
            records = read_fasta(path)
            kind = input_type
            if kind == "auto":
                kind = "nucleotide" if looks_nucleotide(records) else "protein"
            if kind == "nucleotide":
                stage("calling genes")
                proteins, used = call_genes(records, mode)
                (outdir / "proteins").mkdir(exist_ok=True)
                write_fasta(outdir / "proteins" / f"{gid}.faa", proteins)
                calling = f"pyrodigal {used}"
            else:
                proteins, calling = records, "none (protein input)"

            calls = database.search(proteins, cpus=threads, gate=gate_subfamilies, stage=stage)
            rows = _label(gid, calls, database.labels)
            markers.write(rows)
            n_rows += len(rows)
            row = dict(
                genome_id=gid, input=str(path), input_type=kind, gene_calling=calling,
                sequences=len(records),
                length_bp=sum(len(s) for _, s in records) if kind == "nucleotide" else "",
                proteins=len(proteins), marker_genes=len({r["gene_id"] for r in rows}),
                markers=len({(r["namespace"], r["accession"]) for r in rows}),
            )
            genomes.write([row])
            summary.append(row)
            task.advance()
            ui.info(f"[dim]{index:>{width}}/{len(inputs)}[/dim] [bold]{name}[/bold] · "
                    f"{len(proteins):,} proteins · {row['markers']:,} markers · "
                    f"{ui.format_duration(time.monotonic() - began)}")

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
        "parameters": {"mode": mode, "input_type": input_type, "threads": threads,
                       "gate_subfamilies": gate_subfamilies},
        "genomes": len(inputs),
        "rows": n_rows,
    }
    with open(outdir / "giftag_run.json", "w") as handle:
        json.dump(run, handle, indent=2)
        handle.write("\n")
    return {"outdir": outdir, "rows": n_rows, "genomes": summary,
            "seconds": run["seconds"], "unsearchable": len(unsearchable)}


def _label(gid, calls, labels):
    """Translate accepted calls into gifter's namespace/accession, one row per
    gene and marker, keeping the strongest supporting hit."""
    best = {}
    for call in calls:
        found = list(labels.get((call.source, call.profile), ()))
        if call.source == "dbcan" and cazy_family(call.profile) != call.profile:
            # A hit to an official subfamily model (GH43_18) is a hit to its
            # family: the narrower evidence licenses the broader marker. The
            # row keeps the subfamily model as its `profile`.
            found += labels.get(("dbcan", cazy_family(call.profile)), ())
        for namespace, accession in found:
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
