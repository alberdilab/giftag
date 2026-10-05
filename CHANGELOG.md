# giftag changelog

Code, command-line and database-build changes. Newest first. Timestamps are
UTC.

giftag has two versions. The **package version** below follows this file. The
**database format** (`DB_FORMAT` in `giftag/__init__.py`, currently 1) changes
only when a directory written by `giftag build` can no longer be read by
`giftag annotate`; a format change is called out in its entry.

---

## 1.0.0 — 2026-10-05T04:15Z

### First release on PyPI

**What changed.** giftag is published on PyPI, so `pip install giftag` works. A GitHub Actions workflow builds the distributions, runs the tests against the built wheel, and uploads to PyPI through trusted publishing when a GitHub release is published. The source distribution now holds only the package, its offline tests, the README, the licence and this changelog; the benchmarks and the documentation site stay in the repository. The "under development" badge and note are replaced by a PyPI version badge.

**Why.** The quickstart has told users to run `pip install giftag` since the documentation site went up, but the package had never been uploaded. The source distribution was also shipping about 470 KB of benchmark results and documentation sources that an installation never uses.

**Effect.** Packaging and release process only. Annotation, search rules, output files and the database format (still 1) are unchanged from 0.2.0. The version is 1.0.0 because the command line and output files are now treated as stable.

### A dedicated documentation site connects giftag to gifter

**What changed.** A Quarto site now documents the build and annotation commands, output schema, marker-search coverage, source rules, validation limits, and the direct handoff to gifter. It uses gifter's palette, typography and logo, with navigation back to gifter. GitHub Actions renders it on pull requests and deploys it from `main`. The README is now a short entry point to those guides.

**Why.** giftag's CLI and coverage ledger need their own reference, while the interpretation of a GIFT call belongs in gifter's documentation.

**Effect.** Documentation and deployment workflow only. Annotation, search rules and output files are unchanged.

## 0.2.0 — 2026-10-04T05:05Z

### dbCAN-sub is searched on every protein

**What changed.** `annotate` now searches the dbCAN-sub clusters on every
protein. The earlier behaviour, searching a family's clusters only on proteins
with a domain of that family, is the opt-in `--gate-subfamilies`.

**Why.** The gate was meant to keep a subfamily call from standing without its
family. Measured against run_dbcan on *Bacteroides thetaiotaomicron* VPI-5482
it dropped 74 of 157 subfamily calls. Most of that was a bug: the gate looked
for a family hit named exactly `GH43` and ignored hits to dbCAN's official
subfamily models such as `GH43_18`. With that fixed it still dropped 13, and
gifter's subfamily evidence was curated against what run_dbcan emits.

**Effect.** Subfamily calls match run_dbcan on the kept clusters, 157 of 157.
The dbCAN step takes about twice as long; a genome takes about 45 s on 8 cores.
With `--gate-subfamilies` the result is 144 of 157.

### A hit to an official subfamily model counts as its family

**What changed.** When a dbCAN family-library hit is to a subfamily model
(`GH43_18`), giftag also emits the parent family marker (`GH43`). The row keeps
the subfamily model in its `profile` column.

**Why.** run_dbcan reports only the model that wins a region, so a GH43 protein
whose best hit is `GH43_18` carries no `GH43` accession. Fourteen of gifter's 53
CAZy family markers are in families with subfamily models, and those markers
missed every such protein. The narrower hit licenses the broader marker, never
the reverse.

**Effect.** This is the one place giftag's markers differ from run_dbcan's
literal output. On *B. thetaiotaomicron* it raises family marker rows from 104
to 191.

### dbCAN-sub is downloaded by byte range

**What changed.** `build` locates the blocks of the CAZy families gifter uses in
`dbCAN_sub.hmm` with range requests and downloads only those. A table of
profile counts per family ships for the pinned release
(`giftag/data/dbcan_sub_db_v5-2-9_5-5-2026.tsv`); every fetched block is checked
against it, and it supplies run_dbcan's Z. On a mismatch, or for a release
without a table, the whole library is streamed as before.

**Why.** The library is 5.1 GB and gifter needs the clusters of 47 of its 500
families.

**Effect.** 1.25 GB is downloaded instead of 5.1 GB, and the kept profiles are
byte-identical to the streamed result. A full build downloads about 2.9 GB.

### KOfam is fetched over FTP first

**What changed.** KOfam is downloaded from GenomeNet's FTP server, with HTTPS as
the fallback.

**Why.** GenomeNet throttles each HTTPS connection to tens of KB/s, which made
the 1.5 GB archive a day-long download. FTP is about thirty times faster.

**Effect.** KOfam takes about 25 minutes on an ordinary connection.

### Command line on Typer and Rich

**What changed.** The command line is rebuilt on Typer and Rich: grouped option
panels, live download and genome progress bars, elapsed-time log lines, and a
summary table after `build`, `annotate` and `info`. New global options
`--verbose`, `--quiet` and `--log-file PATH` go before the command. `-i`
accepts directories. A `version` command is added.

**Why.** A build downloads for the better part of an hour and an annotation run
may cover thousands of genomes; both need visible progress and a log that can
be kept.

**Effect.** Off a terminal the bars become plain periodic lines. `annotate()`
now returns a summary dictionary instead of a row count. `typer` and `rich` are
new dependencies. `--dbcan-dir` accepts `dbCAN-sub.hmm` as well as
`dbCAN_sub.hmm`. `annotate` writes its tables genome by genome, so an
interrupted run keeps what it finished; `giftag_run.json` marks a complete run.

## 0.1.0 — 2026-10-04T03:13Z

**What changed.** First version. `giftag build` compiles a profile database
from the KOfam, NCBIfam, Pfam and dbCAN releases gifter is curated against,
keeping the profiles gifter's markers need and the dbCAN profiles they compete
with. `giftag annotate` calls genes with pyrodigal, applies each source's own
acceptance rule, and writes the `genome_id`, `gene_id`, `namespace`,
`accession` table gifter reads.

**Why.** gifter's markers span four profile collections, and gifter cannot tell
an unsearched marker from an absent gene.

**Effect.** Database format 1. Every gifter marker gets a recorded status, and
those that cannot be searched are reported.
