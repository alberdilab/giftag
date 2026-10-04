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

Because it keeps only the profiles gifter uses, the database is small and the
search takes seconds to a few minutes per genome.

## Install

```sh
pip install giftag        # Python 3.9+, Linux, macOS or Windows
```

The only dependencies are [pyhmmer](https://github.com/althonos/pyhmmer) and
[pyrodigal](https://github.com/althonos/pyrodigal). No external binaries are
needed.

## `giftag build`

```sh
giftag build [-d DIR] [--gifter-db PATH|URL] [--sources kofam,ncbifam,pfam,dbcan]
```

`build` downloads each source **on your machine**, streams it, and keeps only
what gifter needs. The database goes to `$GIFTAG_DB`, or
`~/.local/share/giftag` by default. The download is about 6.5 GB: the 1.5 GB
KOfam archive and the 4.9 GB dbCAN-sub library. Only the kept profiles are
stored, about 1 GB. Expect an hour or so on an ordinary connection.

KOfam is fetched from GenomeNet's FTP server first. GenomeNet throttles each
HTTPS connection to tens of KB/s, which would make the archive a day-long
download. HTTPS is used only when FTP is blocked.

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
giftag annotate -i FASTA... -o DIR [-d DIR] [-t THREADS] [--mode auto|single|meta]
giftag -i FASTA... -o DIR        # the same
```

Inputs are genome (nucleotide) or protein FASTA files, gzipped or not; the type
is detected per file. Genes are called with pyrodigal: trained on the genome
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

## Design decisions

**Each source keeps its own rule.** KOfam, NCBIfam, Pfam and dbCAN calibrate
their thresholds differently. giftag applies each rule as the source defines it
and reports the score and threshold, instead of inventing a combined
confidence.

**dbCAN competitors are searched too.** run_dbcan keeps only the best-scoring
profile where domains overlap. Searching only the profiles gifter asks for
would let a weaker one win a region a stronger, unrequested one owns. That
would be a false and over-specific subfamily call, so the competing profiles
stay in the database and are searched.

**A dbCAN-sub cluster needs its family.** dbCAN-sub is searched only on
proteins that carry a domain of the cluster's parent family. A subfamily call
narrows a family call, so it should not stand without it. This also keeps the
9,000-profile subfamily search small. It is the one place where giftag is
stricter than run_dbcan.

**dbCAN E-values are run_dbcan's.** run_dbcan computes E-values with Z set to
the profile count of the full library. giftag records that count at build time
and uses it, so its E-values are comparable even though it keeps fewer
profiles.

**Unsearchable markers are reported.** For example, EC numbers are activity
evidence in gifter, not sequence profiles, and some KOs have no KOfam
threshold. These markers are listed in `markers.tsv` and counted at the start
of every run, because gifter will read them as absent.

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
