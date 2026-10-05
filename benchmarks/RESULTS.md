# Publication benchmark results

The benchmark uses eight pinned RefSeq assemblies spanning seven bacteria and
one archaeon (35,848 proteins). Every method receives the same protein
sequences and stable protein identifiers. It measures giftag 0.2.0 at source
commit `ab5769fa216107eef551011ae1a4de6c1c676458`. The marker catalogue
is the pinned public gifter database snapshot 2026.37.1; 1,571 of its 1,576
markers have searchable
profiles in the built giftag database. Three markers lack a usable threshold
and two require unsupported evidence. All timed annotations request eight
CPUs. The [workflow and exact inputs](README.md), [assembly hashes](data/panel-manifest.tsv),
[software environments](provenance/), and [machine-readable results](results/final/)
make the comparisons reproducible.

The comparisons measure **gene–marker concordance**, not accuracy against a
biological truth set. A call is a unique genome, protein, namespace and marker
accession tuple. Jaccard is the shared-call count divided by the union of
calls. Calls absent from both methods do not inflate agreement. The intervals
are descriptive 95% percentile intervals from 10,000 genome-block bootstrap
resamples of the eight genomes. They do not establish performance on all
microbial lineages. Reference calls outside the searchable gifter marker set
are excluded; a full-library method's additional annotations are outside the
question this benchmark measures.

## Search scope

| Source | Full reference search | giftag search |
| --- | ---: | ---: |
| KOfam | 27,864 profiles | 900 retained KO profiles |
| dbCAN family | 875 profiles | All 875 profiles; 53 target marker families |
| dbCAN-sub | 53,411 profiles | 9,025 profiles covering 467 target clusters |
| Pfam | 24,736 profiles in InterProScan 5 or 27,481 in InterProScan 6 | 1 target profile |
| NCBIFAM | 17,431 profiles in InterProScan 5 or 18,054 in InterProScan 6 | 150 target profiles |

The exact [database scope table](results/final/database-scope.tsv) contains
stored byte counts and artifact formats. KOfam's full size is a compressed
archive and the subset size is an uncompressed HMM; their raw byte counts
should not be interpreted as a storage reduction factor. The [InterPro member
inventory](provenance/interpro-member-scope.tsv) records the member library
releases and SHA-256 hashes. InterProScan is run with **only** the Pfam and
NCBIFAM analyses, using the complete libraries for those analyses. Both
InterProScan releases use older member models than giftag's NCBIFAM 20.0 and
Pfam 38.2 models. The extracted full KOfam profile directory occupied
7.67 GB on Mjolnir, versus 204 MB for giftag's concatenated KO HMM; those
directory and file layouts still differ. The full dbCAN-sub HMM is 5.13 GB
versus 1.13 GB for giftag's retained HMMs; the 0.13 GB dbCAN family library
remains complete.
The InterProScan 5/6 full Pfam libraries are 1.88/2.07 GB and the full
NCBIFAM libraries are 2.63/2.72 GB, versus 0.14 MB and 23.8 MB for giftag's
one Pfam and 150 NCBIFAM target models.

## Gene–marker agreement

| Reference | Source | Shared | giftag only | Reference only | Jaccard (95% genome bootstrap interval) |
| --- | --- | ---: | ---: | ---: | ---: |
| KofamScan, full and subset | KO | 2,781 | 0 | 0 | 1.000 (1.000–1.000) |
| run_dbcan, full and subset | CAZy family | 483 | 0 | 0 | 1.000 (1.000–1.000) |
| run_dbcan, full and subset | dbCAN-sub cluster | 71 | 0 | 1 | 0.986 (0.959–1.000) |
| Direct HMMER | NCBIFAM | 239 | 0 | 0 | 1.000 (1.000–1.000) |
| Direct HMMER | Pfam | 1 | 0 | 0 | 1.000 (1.000–1.000) |
| InterProScan 5 | NCBIFAM | 229 | 10 | 3 | 0.946 (0.882–0.979) |
| InterProScan 6 | NCBIFAM | 233 | 6 | 0 | 0.975 (0.926–0.997) |
| InterProScan 5 and 6 | Pfam | 1 | 0 | 0 | 1.000 (1.000–1.000) |

The [complete table](results/final/agreement.tsv) has a separate row for all
12 comparisons, including both dbCAN scopes and both InterProScan releases.
The [publication figure](results/final/benchmark.pdf) is also available as
[editable SVG](results/final/benchmark.svg).

**Figure caption.** (A) Jaccard agreement of protein–marker calls within the
searchable gifter marker set. Whiskers show 95% genome-block bootstrap
intervals; the horizontal axis is expanded to 0.85–1.01. (B) Annotation wall
time for each of eight identical protein FASTA inputs, with eight requested
CPUs per run; the vertical axis is logarithmic. One-time setup and downloads
are excluded, and InterProScan 6 had an initial cache warmup. Full refers
to the complete selected reference libraries, while subset refers to the
targeted library specified in the methods. Points are independent genome
runs, not technical repeats.

The [per-genome](results/final/KofamScan-full-kofam-by-genome.tsv),
[per-marker](results/final/InterProScan5-full-ncbifam-by-marker.tsv), and
discordant-call files identify where differences occur. The Pfam interval
is conditional on the 6,510 bootstrap resamples containing the sole observed
Pfam call; it should not be read as evidence of broad Pfam equivalence.

## Computational effort

| Timed annotation | Median seconds per genome | Sum of eight genome runs, seconds | Maximum single-child RSS, GB |
| --- | ---: | ---: | ---: |
| giftag, all target sources | 50.8 | 435 | 1.62 |
| KofamScan, full | 1,841.8 | 14,157 | 0.85 |
| KofamScan, 900-profile subset | 47.3 | 366 | 0.13 |
| run_dbcan, full | 254.5 | 2,208 | 2.92 |
| run_dbcan, reduced dbCAN-sub library | 63.9 | 543 | 0.92 |
| InterProScan 5, full Pfam and NCBIFAM members | 644.2 | 5,286 | 5.42 |
| InterProScan 6, full Pfam and NCBIFAM members | 547.7 | 4,926 | 1.74 |

All runs used eight requested CPUs. The values are annotation time after
environment and database preparation; the InterProScan 6 initial container
and member-database download was separately warmed and recorded. The
[resource table](results/final/resource-summary.tsv) also reports process CPU
time and seconds per 1,000 proteins; [individual records](results/final/measurements.jsonl)
retain commands, hashes, hosts and job IDs. Child RSS is the largest single
process, which particularly understates programs with concurrent workers.
The measured giftag database build took 59.0 seconds from pre-staged KOfam
and dbCAN inputs, and its output occupied 1.49 GB; that build time is
separate from the annotation values above.

Adding independently timed reference stages gives a serial-workflow
estimate for the same eight genomes:

| Reference stack | Sum of stage wall time, seconds | Ratio to giftag's 435 seconds | Ratio of process CPU time |
| --- | ---: | ---: | ---: |
| Targeted KofamScan + run_dbcan + direct HMMER | 937 | 2.15× | 1.26× |
| Full KofamScan + run_dbcan + InterProScan 5 | 21,652 | 49.75× | 17.24× |
| Full KofamScan + run_dbcan + InterProScan 6 | 21,291 | 48.92× | 17.19× |

These are sums of separately measured runs, not elapsed times from a jointly
scheduled pipeline. The [paired genome values](results/final/pipeline-effort.tsv)
and [summary](results/final/pipeline-effort-summary.tsv) expose the
calculation. The full stacks search many models outside gifter's marker set;
that is the purpose of the full-database effort comparison.

## Discordance adjudication

The direct HMMER control searches giftag's exact Pfam and NCBIFAM profile
bytes with the same thresholds. It therefore isolates giftag's parsing and
cutoff handling from model release differences.

The one dbCAN-sub difference is `GH95_e26` on
`btheta_vpi5482_p004502`. Both methods detect competing GH95 models.
run_dbcan's pyhmmer 0.11.0 reports E-values underflowed to zero for both,
so its overlap resolution selects `GH95_e26`; giftag's pyhmmer 0.12.3
retains distinct E-values and selects `GH95_e1`. The raw scores, E-values,
software versions and targeted reproduction are in the
[adjudication record](results/adjudication/README.md). The same single
difference appears with the full and reduced run_dbcan libraries.

InterProScan 5 uses NCBIFAM 17.0, InterProScan 6 uses 18.0, and giftag uses
20.0. Of giftag's 150 target NCBIFAM models, 140 occur in each older member
library. The [model-presence ledger](provenance/interpro-ncbifam-model-presence.tsv)
distinguishes unavailable models from calls that differ despite a model being
present. All six giftag-only NCBIFAM calls relative to InterProScan 6 use
models absent from NCBIFAM 18.0. InterProScan 5 has ten giftag-only calls:
six use unavailable models and four use models present in its older library.
It also has three reference-only calls, all for `TIGR01001.2`. The latter
seven differences are release-dependent observations; this study does not
attribute each one to a specific change in profile scoring or cutoff.

## Consequences for complete GIFT calls

The pinned gifter evaluation finds 555 complete genome–GIFT calls in the
giftag baseline, across 1,304 assessed genome–GIFT pairs. Replacing one
marker source at a time changes none of those complete calls for KofamScan,
run_dbcan, direct HMMER or InterProScan 6. InterProScan 5 NCBIFAM calls cause
one complete call to be lost: `siroheme_to_heme_b` in
*Clostridium acetobutylicum* ATCC 824. On this genome giftag and
InterProScan 6 call `NF040707.3` and `NF040708.3` on separate proteins;
InterProScan 5 calls neither, although both model identifiers are present in
its older member library. Those two NCBIFAM markers distinguish the required
AhbB and AhbA subunits in the curated route. Substituting all sources
together confirms the same result: the targeted and InterProScan 6 reference
stacks preserve all 555 complete calls, while the InterProScan 5 stack loses
this one call. The two rows in the changed-call table represent this same
biological difference under single-source and combined-stack substitutions.
The [GIFT impact summary](results/final/gift-impact-summary.tsv),
[changed-call table](results/final/gift-discordant.tsv), and pinned gifter
source record the full evaluation. The benchmark does not establish which
specific model or cutoff revision made InterProScan 5 miss the two calls.

## Interpretation

The eight-genome panel captures several bacterial phyla and one archaeon,
but is not a random sample of prokaryotic diversity. The gifter catalogue
has only one Pfam target profile and one observed Pfam call in this panel;
its perfect agreement is correspondingly weak evidence for general Pfam
performance. Results are specific to the pinned marker catalogue, database
releases, genomes and software versions. Timed runs exclude one-time
environment setup, downloads and database construction; those stages are
recorded separately. Each genome was measured once, so runtime differences
can include node and shared-filesystem variation. Peak child RSS is the
largest single child process and can understate simultaneous workers' memory.
Process CPU time is also collected through child-process accounting and may
understate work performed by separately launched workflow tasks. The
[Slurm accounting record](provenance/slurm-accounting.psv) is retained alongside
the measurements.
