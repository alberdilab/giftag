# giftag

[![docs](https://github.com/alberdilab/giftag/actions/workflows/docs.yml/badge.svg)](https://github.com/alberdilab/giftag/actions/workflows/docs.yml)
[![Documentation](https://img.shields.io/badge/docs-website-1f425f.svg)](https://alberdilab.github.io/giftag/)
![Status: under development](https://img.shields.io/badge/status-under%20development-orange.svg)
[![License: MIT-0](https://img.shields.io/badge/license-MIT--0-blue.svg)](LICENSE)

giftag annotates genome and protein FASTA files with the markers [gifter](https://alberdilab.github.io/gifter/) evaluates. It writes the gene-by-marker table gifter reads, together with search evidence and a record of which markers could be searched.

**[Read the giftag documentation](https://alberdilab.github.io/giftag/)** for the [quickstart](https://alberdilab.github.io/giftag/get-started.html), [commands](https://alberdilab.github.io/giftag/commands.html), [output files](https://alberdilab.github.io/giftag/output.html), [marker coverage](https://alberdilab.github.io/giftag/coverage.html), [methods and validation](https://alberdilab.github.io/giftag/methods.html), and the [handoff to gifter](https://alberdilab.github.io/giftag/with-gifter.html).

## Quickstart

Python 3.9 or newer is required. The full database build downloads about 2.9 GB from KOfam, NCBIfam, Pfam and dbCAN onto your machine and uses about 1.5 GB of disk space. Read the [source terms](https://alberdilab.github.io/giftag/coverage.html#source-terms) before building.

```sh
pip install giftag
giftag build
giftag annotate -i genomes/ -o annotations/
giftag info
```

Then use the marker table in R:

```r
markers <- read.delim("annotations/giftag_markers.tsv")
calls <- gifter::evaluate_gifts_community(markers)
```

Each FASTA file is one genome. `giftag_markers.tsv` contains `genome_id`, `gene_id`, `namespace` and `accession` for gifter, plus the profile, rule, score and threshold behind each accepted hit. The built database's `markers.tsv` identifies markers giftag could not search. Keep it and `giftag_run.json` with an analysis: an unsearched marker otherwise looks absent to gifter.

The tools make different claims. giftag observes marker evidence; gifter evaluates whether those observations complete a curated genomic capability. Neither tool infers expression, activity or phenotype from a hit.

## Development

```sh
pip install -e '.[test]'
pytest
```

Documentation source is in [`docs/`](docs/). From the repository root, run `quarto render docs` to build the site locally. The GitHub Pages workflow renders it on pull requests and deploys it from `main`.

## Licensing

giftag's code and original documentation are MIT-0. giftag redistributes no profile libraries. `giftag build` downloads them from their publishers, whose terms apply; see [source terms](https://alberdilab.github.io/giftag/coverage.html#source-terms).
