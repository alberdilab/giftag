# Publication benchmark

See the [results and interpretation](RESULTS.md) and the
[machine-readable final tables](results/final/).
After copying results from Mjolnir, run `python benchmarks/verify_results.py`
from the repository root to check every published table and figure hash.

This directory is the reproducible benchmark for giftag against software that
searches the same profile collections. The eight fixed [RefSeq assemblies](panel.tsv)
include seven bacteria and one archaeon. `fetch_panel.py` downloads NCBI
protein FASTA files and replaces FASTA identifiers with unique, stable
`<genome>_p<index>` IDs. It retains the NCBI packages, original headers,
checksums and protein counts. Every tool receives the **same protein sequences**;
`setup.sbatch` checks the regenerated FASTA hashes against the pinned
[panel manifest](data/panel-manifest.tsv). This study measures annotation, not
gene prediction. The marker set comes from
the public [gifter snapshot](https://github.com/alberdilab/gifter/tree/ce63c15200b64a6614985aedd118da4d92013ea1)
at database version 2026.37.1. `prepare-gifter.sbatch` downloads that exact
commit and verifies the marker TSV (SHA-256
`552fc134bb989ed37e0223ce05a17e0632a217416a197fbd6e36ce1d536667b9`)
and SQLite database. `gifter-env.sbatch` installs gifter 0.7.3 from the pinned
source. The expected hashes are in
[`data/marker-provenance.json`](data/marker-provenance.json). The gifter
source and curated database are fetched into Mjolnir scratch rather than
redistributed by giftag.
The giftag source used on Mjolnir is commit
`ab5769fa216107eef551011ae1a4de6c1c676458`; the
[`source checksum manifest`](provenance/giftag-source-sha256.txt) verifies
the actual files installed there. A later documentation-only repository
commit does not affect these measurements.

## Comparisons

| Collection | Reference method | Full versus subset | Agreement level |
|---|---|---|---|
| KOfam | KofamScan 1.3.0 | Its full 2026-08-02 profile directory and a `.hal` containing exactly giftag's retained KOs | KO per protein |
| dbCAN family and dbCAN-sub | run_dbcan 5.2.9, HMM methods only | Full pinned `db_v5-2-9_5-5-2026` libraries and a subset with giftag's exact HMM bytes; both retain the full family library | CAZy family or cluster per protein |
| Pfam and NCBIFAM | HMMER 3.4 `domtblout` with the same giftag profiles and source cutoffs | giftag subset | Pfam or NCBIFAM accession per protein |
| Pfam and NCBIFAM | InterProScan 6.0.2.2, local computation only | Full Pfam 38.1 and NCBIFAM 18.0 member libraries for InterPro 108.0 | Pfam or NCBIFAM accession per protein |
| Pfam and NCBIFAM | InterProScan 5.75-106.0, Mjolnir module, local computation only | Full Pfam 37.4 and NCBIfam 17.0 member libraries | Pfam or NCBIFAM accession per protein |

The HMMER comparison isolates implementation and cutoff handling because the
profile bytes are identical. InterProScan provides a common complete-pipeline
reference, but its data release can differ from the individual current Pfam
files fetched by giftag; that difference must be recorded when explaining a
discordant call. The installed InterProScan 5 member releases are older than
giftag's NCBIfam 20.0 and Pfam 38.2 models, so it is a release-sensitivity
comparison rather than an implementation control. KofamScan's full and subset
runs use one `ko_list` and one dated KOfam archive. run_dbcan counts profiles
in the supplied library to set its statistical `Z`: its dbCAN-sub subset run
therefore uses 9,025, while giftag preserves the full-library value of
53,411. Interpret subset call differences with this search-space change in
mind. run_dbcan's DIAMOND, gene-cluster and substrate
integration methods are excluded because they ask a different biological
question. [KofamScan documents `.hal` subsets and its adaptive KO thresholds](https://github.com/takaram/kofam_scan),
[run_dbcan documents its HMM-only methods](https://run-dbcan.readthedocs.io/en/latest/user_guide/CAZyme_annotation.html),
and [InterProScan documents analysis selection and local computation](https://interproscan6.readthedocs.io/stable/options/).

`run_dbcan` may let an unrelated dbCAN-sub family displace a retained cluster.
The full-library comparison tests this explicitly. For family markers, an
official dbCAN model such as `GH43_18` also supports its parent `GH43`, as
giftag's documented projection does. All comparisons restrict to markers that
the built giftag database calls `searchable`.

## Metrics and provenance

The agreement unit is a unique `(genome_id, gene_id, namespace, accession)`
tuple. `compare.py` writes intersections, calls unique to each method,
Jaccard agreement and precision/recall **relative to the reference**. These are
concordance measures, not biological sensitivity or specificity. It writes
each discordant tuple and per-genome and per-marker tables. `aggregate.py`
uses 10,000 genome-block bootstrap resamples (seed `20261004`) for a 95%
Jaccard interval. Absence of calls from both tools is excluded from Jaccard;
the searchable-marker ledger gives the denominator for coverage. The aggregate
table records how many bootstrap resamples had at least one call, particularly
important for the single Pfam target marker. These are descriptive intervals
from eight genomes, not confidence limits for all microbes.

`gift_impact.R` replaces one source's giftag calls at a time, then replaces
all five sources with each targeted or full reference stack and reevaluates
the panel with the pinned gifter 0.7.3 source and SQLite database. Its
gained/lost complete GIFT counts show whether call-level discordance changes
the biological summary. No reference is treated as ground truth.

`measure.py` records wall time, user and system CPU time, child-process peak
RSS, exit status, command, input hash, hashes for supplied database files,
host and SLURM job ID. The database and member-library checksum manifests
cover larger libraries separately. Child peak RSS is
the highest single process, so it can understate the combined memory of tools
that launch several workers. SLURM cgroup accounting should be added when
reporting total job memory. `measurements.jsonl` preserves all individual
measurement records and their input hashes. Timed runs exclude one-time environment creation,
profile downloads and database construction; those stages are reported
separately. All timed jobs request eight CPUs. Database size and download
bytes come from the saved provenance files, not estimates. The
[`database-scope.tsv`](data/database-scope.tsv) table labels artifact formats:
KOfam's full archive is compressed, while the giftag subset is an uncompressed
concatenated HMM, so those two byte counts are not a like-for-like storage
ratio. The
[`InterPro member inventory`](provenance/interpro-member-scope.tsv) records
full Pfam and NCBIFAM profile counts, releases, file sizes and SHA-256 hashes.
[`NCBIFAM model presence`](provenance/interpro-ncbifam-model-presence.tsv)
checks each of giftag's 150 target models against both older member releases;
this separates unavailable models from cutoff or model-version disagreements.
`resource-summary.tsv` reports medians, totals and wall seconds per 1,000
input proteins from the eight independent genome measurements; it does not
include database download or build time. `pipeline-effort.tsv` and its summary
add separately measured reference stages as a **serial workflow estimate**:
the targeted stack uses KofamScan subset, run_dbcan subset and direct HMMER;
each full stack uses full KofamScan, full run_dbcan and one InterProScan
release. These sums do not measure a jointly scheduled pipeline and the
InterPro member releases differ from giftag's source models.

## Mjolnir workflow

The scripts in [`mjolnir/`](mjolnir/) use the lab scratch project
`/projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026`, following
that workspace's `AGENTS.md`. Transfer the repository source to
`$project/source`, then submit these stages with `sbatch --dependency=afterok`
for their prerequisites:

1. `prepare-gifter.sbatch` fetches and verifies the pinned public gifter
   snapshot. `setup.sbatch` and `sources.sbatch` create the giftag environment,
   fetch the NCBI panel and pinned full KOfam/dbCAN libraries.
   `reference-envs.sbatch` installs KofamScan and run_dbcan separately.
   `gifter-env.sbatch` installs the pinned gifter source into an isolated R
   environment. Submit `setup.sbatch` and `gifter-env.sbatch` after
   `prepare-gifter.sbatch` succeeds.
2. `build.sbatch` depends on setup and sources. It times the giftag database
   build and saves its ledger and profile checksums.
3. `prepare-kofam.sbatch` depends on the build. It extracts the full KOfam
   archive and makes the exact giftag KO subset `.hal`.
4. `prepare-interproscan.sbatch` pins the InterProScan 6 workflow and Nextflow
   runtime via HTTPS. `run-giftag.sbatch`, `run-kofam.sbatch`,
   `run-dbcan.sbatch`, `run-dbcan-subset.sbatch`, `run-hmmer.sbatch`,
   `run-interproscan.sbatch`, and the
   optional `run-interproscan5.sbatch` each run the same protein panel.
   `BENCH_GENOME=ecoli_mg1655` limits a job to the pilot genome. The first
   InterProScan 6 run should set
   `BENCH_RESULT_ROOT=$project/results/warmup/interproscan`: it fetches member
   databases and containers. Rerun the genome without that override, after
   the cache is populated, for the measured annotation time.
5. Once all eight genomes and comparators have succeeded, run:

   ```sh
   python benchmarks/aggregate.py \
     --project /projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026 \
     --outdir /projects/alberdilab/scratch/jpl786/projects/giftag-benchmark-2026/results/final \
     --interproscan-version 6
   ```

   Use `--interproscan-version 5 6` to require and compare both versions.
   `finalize.sbatch` runs aggregation, GIFT impact and PDF/SVG plotting after
   all required jobs complete; set `BENCH_INTERPROSCAN_VERSIONS=5,6` for both.

`aggregate.py` fails if a required run is missing or failed. After
`finalize.sbatch`, run `bash source/benchmarks/mjolnir/collect-accounting.sh`
on Mjolnir's login node; the Slurm accounting database is unavailable on some
compute nodes. Copy the compact `results/final` and `provenance` directories
back to a durable repository for publication. Keep the large source profiles
in Mjolnir scratch under their upstream terms. Inspect every discordance
against the source model, score,
threshold, release and overlapping competitors before assigning a cause.
The [script checksum manifest](provenance/benchmark-scripts-sha256.txt)
identifies the replay code and panel definition in this repository.

## Sources

- [Aramaki et al. 2020, KofamKOALA and KofamScan](https://pmc.ncbi.nlm.nih.gov/articles/PMC7141845/)
- [Zheng et al. 2023, dbCAN3](https://pmc.ncbi.nlm.nih.gov/articles/PMC10320055/) (method lineage; this benchmark pins run_dbcan 5.2.9 separately)
- [Eddy 2011, HMMER3](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002195)
- [Jones et al. 2014, InterProScan 5](https://pmc.ncbi.nlm.nih.gov/articles/PMC3998142/)
- [Blum et al. 2026, InterProScan 6](https://pmc.ncbi.nlm.nih.gov/articles/PMC13221978/)
- [NCBI Datasets protein package and `--include protein`](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/)
- [KofamScan usage and `.hal` profile subsets](https://github.com/takaram/kofam_scan)
- [run_dbcan HMM methods and output](https://run-dbcan.readthedocs.io/en/latest/user_guide/CAZyme_annotation.html)
- [InterProScan 6 analysis catalogue](https://interproscan6.readthedocs.io/stable/analyses/)
- [InterProScan 6 TSV format](https://interproscan6.readthedocs.io/stable/output/)
