"""Command line: `giftag build`, `giftag annotate` (or bare `giftag -i ...`), `giftag info`."""

import os
import sys
from collections import Counter, defaultdict
from enum import Enum
from pathlib import Path
from typing import Annotated, List, Optional

import typer
from rich.console import Console
from rich.markup import escape
from rich.panel import Panel
from rich.table import Table

from giftag import GiftagError, __version__, ui
from giftag import sources as S

COMMANDS = ("build", "annotate", "info", "version")
TAGLINE = "Annotate genomes with exactly the markers gifter evaluates"

DATABASE_PANEL = "Database"
SOURCES_PANEL = "Sources"
IO_PANEL = "Input and output"
SEARCH_PANEL = "Gene calling and search"
EXECUTION_PANEL = "Execution"
LOGGING_PANEL = "Logging"

_SOURCE_TITLES = {"kofam": "KOfam", "ncbifam": "NCBIfam", "pfam": "Pfam", "dbcan": "dbCAN",
                  "": "(no source)"}
_STATUS_TEXT = {
    "no_threshold": "the source publishes no threshold",
    "missing_from_source": "not in the source release",
    "unsupported_namespace": "no sequence profile exists for this namespace",
    "source_skipped": "source left out of this build",
}


class Mode(str, Enum):
    auto = "auto"
    single = "single"
    meta = "meta"


class InputType(str, Enum):
    auto = "auto"
    nucleotide = "nucleotide"
    protein = "protein"


def default_db():
    if os.environ.get("GIFTAG_DB"):
        return Path(os.environ["GIFTAG_DB"])
    base = os.environ.get("XDG_DATA_HOME") or os.path.join(os.path.expanduser("~"), ".local", "share")
    return Path(base) / "giftag"


def _option(*names, panel, **kwargs):
    return typer.Option(*names, rich_help_panel=panel, **kwargs)


_HELP_REQUESTED = False

app = typer.Typer(
    add_completion=False,
    help=f"{TAGLINE}.",
    no_args_is_help=False,
    pretty_exceptions_enable=False,
    rich_markup_mode="rich",
)

DbOption = Annotated[Optional[Path], _option(
    "--db", "-d", panel=DATABASE_PANEL, metavar="DIR",
    help="Database directory. Defaults to [bold]$GIFTAG_DB[/bold], else ~/.local/share/giftag.")]


@app.callback()
def _configure(
    context: typer.Context,
    verbose: Annotated[bool, _option("--verbose", "-v", panel=LOGGING_PANEL,
                                     help="Show debug messages.")] = False,
    quiet: Annotated[bool, _option("--quiet", "-q", panel=LOGGING_PANEL,
                                   help="Show warnings and errors only; no progress.")] = False,
    log_file: Annotated[Optional[Path], _option(
        "--log-file", panel=LOGGING_PANEL, metavar="PATH",
        help="Also write a UTC-stamped plain-text log to this file.")] = None,
):
    ui.configure(verbose, quiet, log_file)
    context.obj = {"quiet": quiet}
    if context.invoked_subcommand in ("build", "annotate") and not _HELP_REQUESTED:
        ui.info(f"[bold cyan]giftag[/bold cyan] {__version__} · {context.invoked_subcommand}")


@app.command("build", rich_help_panel="Workflows")
def build_command(
    context: typer.Context,
    db: DbOption = None,
    gifter_db: Annotated[str, _option(
        "--gifter-db", panel=DATABASE_PANEL, metavar="PATH|URL",
        help="gifter SQLite database or marker TSV. Defaults to gifter main on GitHub.")] = S.GIFTER_DB_URL,
    force: Annotated[bool, _option("--force", panel=DATABASE_PANEL,
                                   help="Replace an existing database.")] = False,
    sources: Annotated[str, _option(
        "--sources", panel=SOURCES_PANEL,
        help="Comma-separated sources to include.")] = ",".join(S.SOURCES),
    kofam_release: Annotated[str, _option(
        "--kofam-release", panel=SOURCES_PANEL,
        help="'current', or a GenomeNet archive date such as 2026-08-02.")] = S.KOFAM_DEFAULT_RELEASE,
    kofam_dir: Annotated[Optional[Path], _option(
        "--kofam-dir", panel=SOURCES_PANEL, metavar="DIR",
        help="Local KOfam copy: ko_list and profiles/ or profiles.tar.gz.")] = None,
    ncbifam_release: Annotated[str, _option("--ncbifam-release", panel=SOURCES_PANEL,
                                            help="NCBIfam release.")] = S.NCBIFAM_DEFAULT_RELEASE,
    ncbifam_dir: Annotated[Optional[Path], _option(
        "--ncbifam-dir", panel=SOURCES_PANEL, metavar="DIR",
        help="Local NCBIfam copy: hmm_PGAP.tsv and hmm_PGAP.HMM/.")] = None,
    pfam_dir: Annotated[Optional[Path], _option(
        "--pfam-dir", panel=SOURCES_PANEL, metavar="DIR",
        help="Local Pfam profiles, as PF00001.hmm[.gz] files.")] = None,
    dbcan_release: Annotated[str, _option("--dbcan-release", panel=SOURCES_PANEL,
                                          help="dbCAN release.")] = S.DBCAN_DEFAULT_RELEASE,
    dbcan_dir: Annotated[Optional[Path], _option(
        "--dbcan-dir", panel=SOURCES_PANEL, metavar="DIR",
        help="Local run_dbcan database: dbCAN.hmm and dbCAN_sub.hmm.")] = None,
    threads: Annotated[int, _option("--threads", "-t", panel=EXECUTION_PANEL,
                                    help="Parallel downloads.")] = 8,
):
    """Download the pinned sources and compile the profile database.

    Each source is downloaded [bold]on this machine[/bold] and only the profiles
    gifter's markers need are kept. giftag redistributes nothing.
    """
    from giftag.build import build
    db = db or default_db()
    build(db, gifter_db, [s.strip() for s in sources.split(",") if s.strip()],
          kofam_release=kofam_release, kofam_dir=kofam_dir,
          ncbifam_release=ncbifam_release, ncbifam_dir=ncbifam_dir,
          pfam_dir=pfam_dir, dbcan_release=dbcan_release, dbcan_dir=dbcan_dir,
          force=force, threads=threads)
    if not context.obj["quiet"]:
        _database_report(db, ui.console, built=True)


@app.command("annotate", rich_help_panel="Workflows")
def annotate_command(
    context: typer.Context,
    input: Annotated[List[Path], _option(
        "--input", "-i", panel=IO_PANEL, metavar="FASTA...",
        help="Genome or protein FASTA files, gzipped or not, or directories of them. "
             "Several may follow one [bold]-i[/bold].")],
    outdir: Annotated[Path, _option("--outdir", "-o", panel=IO_PANEL, metavar="DIR",
                                    help="Output directory.")],
    db: DbOption = None,
    input_type: Annotated[InputType, _option(
        "--input-type", panel=IO_PANEL,
        help="Treat inputs as nucleotide or protein instead of detecting each file.")] = InputType.auto,
    mode: Annotated[Mode, _option(
        "--mode", panel=SEARCH_PANEL,
        help="Gene calling: train on each genome ([bold]single[/bold]), use metagenomic models "
             "([bold]meta[/bold]), or single from 100 kb up ([bold]auto[/bold]).")] = Mode.auto,
    gate_subfamilies: Annotated[bool, _option(
        "--gate-subfamilies", panel=SEARCH_PANEL,
        help="Search dbCAN-sub clusters only on proteins with a domain of the cluster's "
             "family: faster, but stricter than run_dbcan.")] = False,
    threads: Annotated[int, _option("--threads", "-t", panel=EXECUTION_PANEL,
                                    help="CPU threads for the search; 0 uses all.")] = 0,
):
    """Annotate genomes or proteins and write the marker table gifter reads.

    [bold]giftag -i ... -o ...[/bold] is shorthand for this command.
    """
    from giftag.annotate import annotate
    from giftag.genes import expand_inputs
    result = annotate(expand_inputs(input), outdir, db or default_db(), threads=threads,
                      mode=mode.value, input_type=input_type.value,
                      gate_subfamilies=gate_subfamilies)
    if not context.obj["quiet"]:
        _annotation_report(result, ui.console)


@app.command("info", rich_help_panel="Information")
def info_command(db: DbOption = None):
    """Describe a giftag database: sources, releases, markers and terms of use."""
    _database_report(db or default_db(), Console(highlight=False), built=False)


@app.command("version", rich_help_panel="Information")
def version_command():
    """Print the installed giftag version."""
    typer.echo(f"giftag {__version__}")


# -- reports ------------------------------------------------------------------

def _table(*columns):
    """A borderless table; a column named `Name>r` is right-aligned."""
    table = Table(show_header=True, header_style="bold", box=None, pad_edge=False,
                  padding=(0, 2, 0, 0))
    for column in columns:
        name, _, right = column.partition(">")
        first = not table.columns
        table.add_column(name, justify="right" if right else "left",
                         style="bold cyan" if first else None, no_wrap=first or bool(right))
    return table


def _directory_size(path):
    return sum(f.stat().st_size for f in Path(path).iterdir() if f.is_file())


def _database_report(db, console, built):
    from giftag.database import Database
    database = Database(db)
    manifest = database.manifest
    gifter = manifest["gifter"]
    version = gifter.get("gifter_db_version") or escape(Path(str(gifter["source"])).name)
    console.print()
    console.print(Panel.fit(
        f"[bold cyan]giftag database[/bold cyan]{' ready' if built else ''}\n"
        f"gifter database {version} · {gifter['markers']:,} markers · "
        f"{ui.format_bytes(_directory_size(db))}\n"
        f"built {manifest['built_utc']} by giftag {manifest['giftag_version']}",
        border_style="cyan"))
    console.print(f"Location: [bold]{escape(str(Path(db).resolve()))}[/bold]", soft_wrap=True)

    by_source = defaultdict(Counter)
    for marker in database.markers:
        by_source[S.NAMESPACE_SOURCE.get(marker["namespace"], "")][marker["status"]] += 1
    table = _table("Source", "Release", "Markers>r", "Searchable>r", "Not searchable>r")
    for source in (*S.SOURCES, ""):
        counts = by_source.get(source)
        if not counts:
            continue
        total = sum(counts.values())
        missing = total - counts["searchable"]
        table.add_row(
            _SOURCE_TITLES[source], str(manifest["sources"].get(source, {}).get("release", "–")),
            f"{total:,}", f"{counts['searchable']:,}",
            f"[yellow]{missing:,}[/yellow]" if missing else "0")
    console.print(table)

    unsearchable = database.unsearchable
    if unsearchable:
        console.print(
            f"\n[yellow]{len(unsearchable):,} markers cannot be searched[/yellow] and will read "
            "as absent in gifter:")
        reasons = _table("Namespace", "Markers>r", "Reason", "Examples")
        grouped = defaultdict(list)
        for marker in unsearchable:
            grouped[(marker["namespace"], marker["status"])].append(marker["accession"])
        for (namespace, status), accessions in sorted(grouped.items()):
            shown = ", ".join(accessions[:4]) + (" …" if len(accessions) > 4 else "")
            reasons.add_row(namespace, f"{len(accessions):,}", _STATUS_TEXT.get(status, status), shown)
        console.print(reasons)
        console.print("[dim]Full list: markers.tsv in the database directory[/dim]")

    if built:
        console.print("\nNext: [bold]giftag -i genomes/*.fa.gz -o annotations/[/bold]")
    else:
        console.print("\n[bold]Terms of use[/bold] [dim](giftag redistributes none of these profiles)[/dim]")
        terms = Table(show_header=False, box=None, pad_edge=False, padding=(0, 2, 0, 0))
        terms.add_column(style="bold cyan", no_wrap=True)
        terms.add_column()
        for source, meta in manifest["sources"].items():
            terms.add_row(_SOURCE_TITLES.get(source, source), meta.get("terms", ""))
        console.print(terms)


def _annotation_report(result, console):
    genomes = result["genomes"]
    console.print()
    if len(genomes) <= 25:
        table = _table("Genome", "Input", "Proteins>r", "Marker genes>r", "Markers>r")
        for row in genomes:
            table.add_row(escape(row["genome_id"]), row["input_type"], f"{row['proteins']:,}",
                          f"{row['marker_genes']:,}", f"{row['markers']:,}")
        console.print(table)
    else:
        markers = sorted(row["markers"] for row in genomes)
        console.print(
            f"[bold]{len(genomes):,} genomes[/bold] · markers per genome: "
            f"min {markers[0]:,}, median {markers[len(markers) // 2]:,}, max {markers[-1]:,}")
    console.print(Panel.fit(
        f"[bold cyan]{result['rows']:,} marker rows[/bold cyan] for {len(genomes):,} "
        f"genome{'s' if len(genomes) != 1 else ''} in {ui.format_duration(result['seconds'])}",
        border_style="cyan"))
    files = Table(show_header=False, box=None, pad_edge=False, padding=(0, 2, 0, 0))
    files.add_column(style="bold cyan", no_wrap=True)
    files.add_column()
    files.add_row("giftag_markers.tsv", "the table gifter reads")
    files.add_row("giftag_genomes.tsv", "per-genome summary")
    files.add_row("giftag_run.json", "database and parameters of this run")
    if any(row["input_type"] == "nucleotide" for row in genomes):
        files.add_row("proteins/", "predicted proteins, one file per genome")
    console.print(f"Written to [bold]{escape(str(result['outdir']))}[/bold]", soft_wrap=True)
    console.print(files)


def show_main_help(console=None):
    console = console or Console(highlight=False)
    console.print(Panel.fit(f"[bold cyan]giftag[/bold cyan] {__version__}\n{TAGLINE}",
                            border_style="cyan"))
    table = _table("Command", "Purpose")
    table.add_row("giftag build [OPTIONS]", "Download the pinned sources and compile the database")
    table.add_row("giftag annotate -i FASTA... -o DIR", "Annotate genomes or proteins")
    table.add_row("giftag -i FASTA... -o DIR", "The same, in short")
    table.add_row("giftag info", "Describe the database")
    table.add_row("giftag version", "Print the installed version")
    console.print(table)
    console.print("\nUse [bold]giftag COMMAND --help[/bold] for a command's options.")
    console.print("Logging options precede the command: [bold]--verbose[/bold], "
                  "[bold]--quiet[/bold], or [bold]--log-file PATH[/bold].")
    return 0


# -- entry point -----------------------------------------------------------------

def normalize(argv):
    """Rewrite the two conveniences Typer does not parse itself: the bare
    `giftag -i ...` shorthand, and several files after one `-i`."""
    argv = list(argv)
    index = 0
    while index < len(argv) and argv[index] in ("-v", "--verbose", "-q", "--quiet", "--log-file"):
        index += 2 if argv[index] == "--log-file" else 1
    if index < len(argv) and argv[index] not in COMMANDS and argv[index] not in ("-h", "--help"):
        argv.insert(index, "annotate")
    expanded = []
    position = 0
    while position < len(argv):
        token = argv[position]
        expanded.append(token)
        position += 1
        if token in ("-i", "--input"):
            first = True
            while position < len(argv) and not (
                    argv[position].startswith("-") and len(argv[position]) > 1):
                if not first:
                    expanded.append(token)
                expanded.append(argv[position])
                first = False
                position += 1
    return expanded


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv in (["--help"], ["-h"], ["help"]):
        return show_main_help()
    if argv == ["--version"]:
        argv = ["version"]
    global _HELP_REQUESTED
    _HELP_REQUESTED = any(token in ("--help", "-h") for token in argv)
    try:
        app(args=normalize(argv), prog_name="giftag")
    except SystemExit as exit_:
        return exit_.code if isinstance(exit_.code, int) else (1 if exit_.code else 0)
    except GiftagError as error:
        ui.console.print(f"[bold red]Error:[/bold red] {escape(str(error))}")
        return 1
    except KeyboardInterrupt:
        ui.console.print("[yellow]Cancelled.[/yellow]")
        return 130
    return 0
