"""Pinned upstream releases, where to fetch them, and the rules each one uses.

giftag redistributes none of these profiles. `giftag build` downloads each
release on the user's machine and keeps only the profiles gifter needs, so the
terms of each source apply to the user exactly as they would for the source's
own tool.

The default releases are the ones gifter's curation was done against
(`inst/extdata/database-source/SOURCES.md` in gifter). KOfam has no gifter pin,
so it defaults to GenomeNet's current release and records the date it got.
"""

GIFTER_DB_URL = (
    "https://raw.githubusercontent.com/alberdilab/gifter/main/inst/extdata/gifter.sqlite"
)

# GenomeNet throttles each HTTPS connection to tens of KB/s, which makes the
# 1.5 GB profile archive a day-long download; its FTP server is about thirty
# times faster. FTP is tried first and HTTPS is the fallback for networks that
# block FTP.
KOFAM_MIRRORS = ("ftp://ftp.genome.jp/pub/db/kofam", "https://www.genome.jp/ftp/db/kofam")
KOFAM_DEFAULT_RELEASE = "current"

NCBIFAM_BASE = "https://ftp.ncbi.nlm.nih.gov/hmm"
NCBIFAM_DEFAULT_RELEASE = "20.0"

DBCAN_BASE = "https://dbcan.s3.us-west-2.amazonaws.com"
DBCAN_DEFAULT_RELEASE = "db_v5-2-9_5-5-2026"

PFAM_ENTRY_URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}?annotation=hmm"
INTERPRO_API_URL = "https://www.ebi.ac.uk/interpro/api/"

# run_dbcan defaults (dbcan/parameter.py, OVERLAP_RATIO_THRESHOLD in
# dbcan/constants/process_utils_constants.py). Both the dbCAN family library
# and dbCAN-sub use the same three values.
DBCAN_EVALUE = 1e-15
DBCAN_COVERAGE = 0.35
DBCAN_OVERLAP = 0.5

SOURCES = ("kofam", "ncbifam", "pfam", "dbcan")

# Which source serves each gifter marker namespace. EC is activity evidence in
# gifter, not a sequence profile, so no source can search for it.
NAMESPACE_SOURCE = {
    "KO": "kofam",
    "NCBIFAM": "ncbifam",
    "TIGRFAM": "ncbifam",
    "PFAM": "pfam",
    "CAZY": "dbcan",
}

TERMS = {
    "kofam": (
        "KOfam, KEGG/GenomeNet. KEGG is not a public database: academic use is "
        "free, non-academic use needs a licence from Pathway Solutions. "
        "Downloaded on this machine and not redistributed by giftag."
    ),
    "ncbifam": (
        "NCBIfam. NCBI-built models are US Government works; TIGRFAM-origin "
        "models are CC BY-SA 4.0. Cite Li et al. 2021, NAR 49:D1020."
    ),
    "pfam": "Pfam via InterPro, CC0. Cite the Pfam and InterPro release papers.",
    "dbcan": (
        "dbCAN (HMMs built from CAZy families). The dbCAN Open Data release "
        "states no restrictions on use. Cite dbCAN3 and CAZy."
    ),
}


def kofam_urls(release):
    """Candidate `(ko_list, profiles)` locations, in the order to try them."""
    bases = [base if release == "current" else f"{base}/archives/{release}" for base in KOFAM_MIRRORS]
    return [f"{b}/ko_list.gz" for b in bases], [f"{b}/profiles.tar.gz" for b in bases]


def ncbifam_table_url(release):
    return f"{NCBIFAM_BASE}/{release}/hmm_PGAP.tsv"


def ncbifam_profile_url(release, accession):
    return f"{NCBIFAM_BASE}/{release}/hmm_PGAP.HMM/{accession}.HMM"


def dbcan_urls(release):
    base = f"{DBCAN_BASE}/{release}"
    return f"{base}/dbCAN.hmm", f"{base}/dbCAN_sub.hmm"
