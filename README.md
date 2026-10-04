# giftag

giftag annotates genomes with exactly the markers
[gifter](https://github.com/alberdilab/gifter) evaluates, and writes them in
the table gifter reads.

```sh
giftag build                                  # once: compile the profile database
giftag -i genomes/*.fa.gz -o annotations/     # annotate
```

```r
markers <- read.delim("annotations/giftag_markers.tsv")
calls <- gifter::evaluate_gifts_community(markers)
```

## Why a dedicated annotator

gifter's evidence comes from four profile collections: KOfam, dbCAN (CAZy
families and dbCAN-sub clusters), NCBIfam and Pfam. A genome annotated with
only one of them leaves the rest of gifter's markers unobserved, and gifter
cannot distinguish an unobserved marker from an absent gene. giftag reads the
marker list from gifter's own database, searches every marker it can, and says
which markers it cannot search and why.

Because it keeps only the profiles gifter uses, the database is about 1.5 GB
and a genome takes under a minute on 8 cores.

## Install

```sh
pip install giftag        # Python 3.9+, Linux, macOS or Windows
```

The search runs on [pyhmmer](https://github.com/althonos/pyhmmer) and gene
calling on [pyrodigal](https://github.com/althonos/pyrodigal); the command line
uses Typer and Rich. No external binaries are needed.

## `giftag build`

```sh
giftag build [-d DIR] [--gifter-db PATH|URL] [--sources kofam,ncbifam,pfam,dbcan]
```

`build` downloads each source **on your machine** and keeps only what gifter
needs. The database goes to `$GIFTAG_DB`, or `~/.local/share/giftag` by
default, and takes about 1.5 GB. The download is about 2.9 GB:

- **KOfam**, 1.5 GB: the archive is streamed and only gifter's KOs are kept.
  It is fetched from GenomeNet's FTP server first, because GenomeNet throttles
  each HTTPS connection to tens of KB/s. HTTPS is the fallback when FTP is
  blocked.
- **dbCAN-sub**, about 1.2 GB of a 5.1 GB library: giftag locates the blocks
  of the CAZy families gifter uses with byte-range requests and downloads only
  those. Each block is checked against the profile count giftag ships for the
  pinned release; on any mismatch, or for a release giftag has no counts for,
  the whole library is streamed instead.
- **dbCAN families**, 130 MB, and a few hundred small NCBIfam and Pfam files.

| Source | Release (default) | Kept | Acceptance rule |
|---|---|---|---|
| KOfam | GenomeNet current, or `--kofam-release 2026-08-02` | the KOs gifter uses | per-KO adaptive threshold from `ko_list`, on the full or best-domain score (KofamScan's rule) |
| NCBIfam / TIGRFAM | `hmm_PGAP/20.0` (gifter's pin) | the accessions gifter uses | trusted cutoffs (`--cut_tc`) |
| Pfam | InterPro current | the accessions gifter uses | gathering cutoffs (`--cut_ga`) |
| dbCAN | `db_v5-2-9_5-5-2026` (gifter's pin) | the whole family library, plus every dbCAN-sub cluster of each family gifter uses | run_dbcan's rule: i-Evalue < 1e-15, HMM coverage > 0.35, overlap resolution |

If you already have a copy, use `--kofam-dir`, `--ncbifam-dir`, `--pfam-dir`
or `--dbcan-dir` instead of downloading. A run_dbcan database directory works
as is.

The database directory holds:

- `giftag.json`: the gifter database version, each source's release, URL and
  SHA-256 checksum, and its terms of use;
- `markers.tsv`: every gifter marker with its status: `searchable`, or why not
  (`no_threshold`, `missing_from_source`, `unsupported_namespace`,
  `source_skipped`);
- `profiles.tsv`: the rule and threshold for every searchable profile;
- one `.hmm` file per source.

`giftag info` summarises the database.

## `giftag annotate`

```sh
giftag annotate -i FASTA... -o DIR [-d DIR] [-t THREADS] [--mode auto|single|meta] [--gate-subfamilies]
giftag -i FASTA... -o DIR        # the same
```

Inputs are genome (nucleotide) or protein FASTA files, gzipped or not; the type
is detected per file. Several files may follow one `-i`, and a directory stands
for the FASTA files inside it, which avoids shell argument limits for
thousands of genomes. Genes are called with pyrodigal: trained on the genome
from 100 kb up, and with metagenomic models below that or with `--mode meta`.
Each file is one genome, and its file name minus extensions is its `genome_id`.

Outputs:

- `giftag_markers.tsv`: one row per gene and marker, with `genome_id`,
  `gene_id`, `namespace` and `accession`, which is what gifter reads. It also
  carries the source, profile, rule, score, E-value, threshold, HMM coverage and
  position of the hit that supports the row;
- `giftag_genomes.tsv`: per genome, the gene calling used and the counts;
- `giftag_run.json`: the database and parameters that produced the run;
- `proteins/<genome_id>.faa`: predicted proteins, for nucleotide input.

## Progress and logs

On a terminal, downloads show a live bar with size, speed and time remaining,
and annotation shows the genome being processed and its current search stage.
In a job log or a pipe the bars become plain lines. Every run ends with a
summary table.

Logging options go before the command:

```sh
giftag --log-file run.log -i genomes/ -o annotations/   # also write a UTC-stamped log
giftag --quiet -i genomes/ -o annotations/              # warnings and errors only
giftag --verbose build
```

## Design decisions

**Each source keeps its own rule.** KOfam, NCBIfam, Pfam and dbCAN calibrate
their thresholds differently. giftag applies each rule as the source defines it
and reports the score and threshold, instead of inventing a combined
confidence.

**dbCAN competitors are searched too.** run_dbcan keeps only the best-scoring
profile where domains overlap. Searching only the profiles gifter asks for
would let a weaker one win a region a stronger, unrequested one owns. That
would be a false and over-specific subfamily call, so the competing profiles
stay in the database and are searched: the whole family library, and every
dbCAN-sub cluster of each family gifter uses.

**dbCAN E-values are run_dbcan's.** run_dbcan computes E-values with Z set to
the profile count of the full library. giftag records that count at build time
and uses it, so its E-values match even though it keeps fewer profiles.

**A subfamily model counts as its family.** dbCAN's family library holds
models for official CAZy subfamilies (`GH43_18`) beside the family models, and
run_dbcan reports whichever wins a region. A `GH43_18` protein is a GH43
protein, so giftag also emits the `GH43` marker for it and keeps `GH43_18` in
the row's `profile` column. This is the one place giftag's markers differ from
run_dbcan's literal output. Without it, gifter's family markers would miss
every protein whose best hit is a subfamily model.

**Subfamily gating is optional.** By default dbCAN-sub is searched on every
protein, as run_dbcan does. `--gate-subfamilies` searches a family's clusters
only on proteins with a domain of that family. It more than halves the dbCAN
time and is stricter than run_dbcan: it drops clusters that hit a protein with
no family-level hit.

**Unsearchable markers are reported.** For example, EC numbers are activity
evidence in gifter, not sequence profiles, and some KOs have no KOfam
threshold. These markers are listed in `markers.tsv` and counted at the start
of every run, because gifter will read them as absent.

## Validation

Checked on *Escherichia coli* K-12 MG1655 and *Bacteroides thetaiotaomicron*
VPI-5482, gene-called by giftag, against gifter database 2026.35.1.

| Check | Result |
|---|---|
| gifter markers searchable | 1,569 of 1,574 |
| dbCAN family calls vs run_dbcan (*B. thetaiotaomicron*) | 458 of 458 gene-family pairs identical |
| dbCAN-sub calls vs run_dbcan on the kept clusters, same Z | 157 of 157 identical |
| dbCAN-sub with `--gate-subfamilies` | 144 of 157 |
| Speed, 8 cores | about 45 s per genome |

The five markers that cannot be searched are three KOs for which KOfam
publishes no threshold (K19098, K19099, K24159) and two EC numbers.

Not yet checked:

- **KOfam against KofamScan itself.** giftag applies KofamScan's rule with the
  same HMMER engine, but the two have not been run side by side.
- **dbCAN-sub against the full library.** The comparison used the clusters
  giftag keeps. A cluster from a family gifter does not use could still
  out-compete a kept one in run_dbcan's overlap resolution.

KOfam's thresholds are conservative, and giftag reproduces that. For example,
*E. coli*'s dihydroorotase PyrC scores 192.5 against K01465, whose threshold is
320.5, so it is not called. Such misses come from the source, not from giftag.

## Licensing

giftag's code is MIT-0. giftag redistributes no profiles: `giftag build`
downloads them from their publishers, so each source's terms apply to you as
they would for the source's own tool.

- **KOfam / KEGG**: academic use is free. Non-academic use, and academic use
  that provides a service to others, needs a licence from Pathway Solutions
  ([KEGG legal](https://www.kegg.jp/kegg/legal.html)). Without one, build with
  `--sources ncbifam,pfam,dbcan`.
- **dbCAN**: the [Open Data release](https://registry.opendata.aws/run_dbcan/)
  states no restrictions on use. Cite dbCAN3 and CAZy.
- **NCBIfam**: NCBI-built models are public domain; TIGRFAM-origin models are
  CC BY-SA 4.0.
- **Pfam**: CC0.
