"""Command line: `giftag build`, `giftag annotate` (or bare `giftag -i ...`), `giftag info`."""

import argparse
import os
import sys
from collections import Counter
from pathlib import Path

from giftag import GiftagError, __version__
from giftag import sources as S

COMMANDS = ("build", "annotate", "info")


def default_db():
    if os.environ.get("GIFTAG_DB"):
        return os.environ["GIFTAG_DB"]
    base = os.environ.get("XDG_DATA_HOME") or os.path.join(os.path.expanduser("~"), ".local", "share")
    return os.path.join(base, "giftag")


def parser():
    top = argparse.ArgumentParser(
        prog="giftag",
        description="Annotate genomes with exactly the markers gifter evaluates.",
        epilog="`giftag -i ...` is shorthand for `giftag annotate -i ...`.",
    )
    top.add_argument("--version", action="version", version=f"giftag {__version__}")
    commands = top.add_subparsers(dest="command", metavar="COMMAND")

    build = commands.add_parser(
        "build", help="download the pinned sources and compile the profile database",
        description="Download each pinned source on this machine and keep the profiles "
                    "gifter's markers need. Nothing is redistributed by giftag.")
    build.add_argument("-d", "--db", default=default_db(),
                       help="database directory (default: $GIFTAG_DB or %(default)s)")
    build.add_argument("--gifter-db", default=S.GIFTER_DB_URL, metavar="PATH|URL",
                       help="gifter SQLite database or marker TSV (default: gifter main on GitHub)")
    build.add_argument("--sources", default=",".join(S.SOURCES),
                       help="comma-separated sources to include (default: %(default)s)")
    build.add_argument("--kofam-release", default=S.KOFAM_DEFAULT_RELEASE,
                       help="'current' or a GenomeNet archive date such as 2026-08-02")
    build.add_argument("--kofam-dir", help="use a local KOfam copy (ko_list and profiles/ or profiles.tar.gz)")
    build.add_argument("--ncbifam-release", default=S.NCBIFAM_DEFAULT_RELEASE)
    build.add_argument("--ncbifam-dir", help="use a local NCBIfam copy (hmm_PGAP.tsv and hmm_PGAP.HMM/)")
    build.add_argument("--pfam-dir", help="use local Pfam profiles (PF00001.hmm[.gz] files)")
    build.add_argument("--dbcan-release", default=S.DBCAN_DEFAULT_RELEASE)
    build.add_argument("--dbcan-dir", help="use a local run_dbcan database (dbCAN.hmm, dbCAN_sub.hmm)")
    build.add_argument("-t", "--threads", type=int, default=8, help="parallel downloads")
    build.add_argument("--force", action="store_true", help="replace an existing database")

    annotate = commands.add_parser(
        "annotate", help="annotate genomes or proteins",
        description="Call genes if needed, search the giftag database, and write a marker "
                    "table gifter reads directly.")
    annotate.add_argument("-i", "--input", nargs="+", required=True, metavar="FASTA",
                          help="genome (nucleotide) or protein FASTA files, optionally gzipped")
    annotate.add_argument("-o", "--outdir", required=True, help="output directory")
    annotate.add_argument("-d", "--db", default=default_db(), help="database directory")
    annotate.add_argument("-t", "--threads", type=int, default=0,
                          help="CPU threads for the search (default: all)")
    annotate.add_argument("--mode", choices=("auto", "single", "meta"), default="auto",
                          help="gene calling: train on each genome (single), use metagenomic "
                               "models (meta), or single from 100 kb up (auto)")
    annotate.add_argument("--input-type", choices=("auto", "nucleotide", "protein"), default="auto")

    info = commands.add_parser("info", help="describe a giftag database")
    info.add_argument("-d", "--db", default=default_db(), help="database directory")
    return top


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] not in COMMANDS and argv[0] not in ("-h", "--help", "--version"):
        argv.insert(0, "annotate")
    args = parser().parse_args(argv)
    if args.command is None:
        parser().print_help()
        return 1
    try:
        if args.command == "build":
            from giftag.build import build
            build(args.db, args.gifter_db, [s.strip() for s in args.sources.split(",") if s.strip()],
                  kofam_release=args.kofam_release, kofam_dir=args.kofam_dir,
                  ncbifam_release=args.ncbifam_release, ncbifam_dir=args.ncbifam_dir,
                  pfam_dir=args.pfam_dir, dbcan_release=args.dbcan_release,
                  dbcan_dir=args.dbcan_dir, force=args.force, threads=args.threads)
        elif args.command == "annotate":
            from giftag.annotate import annotate
            annotate(args.input, args.outdir, args.db, threads=args.threads,
                     mode=args.mode, input_type=args.input_type)
        else:
            _info(args.db)
    except GiftagError as error:
        print(f"giftag: error: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        return 130
    return 0


def _info(db):
    from giftag.database import Database
    database = Database(db)
    manifest = database.manifest
    gifter = manifest["gifter"]
    print(f"giftag database {Path(db).resolve()}")
    print(f"  built {manifest['built_utc']} by giftag {manifest['giftag_version']}")
    print(f"  gifter database {gifter.get('gifter_db_version', '?')} ({gifter['markers']} markers)")
    for name, meta in manifest["sources"].items():
        print(f"  {name}: release {meta.get('release')}")
    print("  markers by status:")
    for status, n in manifest["marker_status"].items():
        print(f"    {status}: {n}")
    unsearchable = Counter((m["namespace"], m["status"]) for m in database.unsearchable)
    if unsearchable:
        print("  not searchable, by namespace:")
        for (namespace, status), n in sorted(unsearchable.items()):
            print(f"    {namespace} {status}: {n}")
    print("  terms of use (giftag redistributes none of these profiles):")
    for name, meta in manifest["sources"].items():
        print(f"    {name}: {meta.get('terms', '')}")
