"""`giftag build`: compile the profile database gifter's markers need.

For each source the builder streams the pinned upstream release, keeps the
profiles gifter's marker list asks for, and records for every marker whether
it can be searched and, if not, why. A marker giftag cannot search is reported,
never silently dropped, because to gifter an unsearched marker looks exactly
like an absent gene.

The one deliberate exception to "keep only what gifter asks for" is dbCAN.
run_dbcan resolves overlapping domains by keeping the best-scoring profile, so
the profiles a requested family or subfamily competes against must be searched
too, or a weaker requested profile would win a region a stronger unrequested
one owns. The full dbCAN family library is therefore kept, and for dbCAN-sub
every cluster of each requested cluster's parent family.
"""

import csv
import io
import json
import os
import shutil
import tarfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import pyhmmer

from giftag import DB_FORMAT, GiftagError, __version__
from giftag import sources as S
from giftag.fetch import Stream, iter_hmm_records, log, open_first, read_all, text_lines
from giftag.markers import load_markers

PROFILE_COLUMNS = ("source", "profile", "namespace", "accession", "rule", "threshold")
MARKER_COLUMNS = ("namespace", "accession", "name", "source", "status", "profile", "detail")


class _Ledger:
    """Searchable profiles and the status of every requested marker."""

    def __init__(self):
        self.profiles = []
        self.status = {}

    def searchable(self, marker, source, profile, rule, threshold):
        self.profiles.append(
            dict(source=source, profile=profile, namespace=marker.namespace,
                 accession=marker.accession, rule=rule, threshold=threshold)
        )
        self.status[marker] = (source, "searchable", profile, "")

    def refuse(self, marker, source, status, detail, profile=""):
        self.status[marker] = (source, status, profile, detail)


def build(db_dir, gifter_db=S.GIFTER_DB_URL, sources=S.SOURCES, *,
          kofam_release=S.KOFAM_DEFAULT_RELEASE, kofam_dir=None,
          ncbifam_release=S.NCBIFAM_DEFAULT_RELEASE, ncbifam_dir=None,
          pfam_dir=None,
          dbcan_release=S.DBCAN_DEFAULT_RELEASE, dbcan_dir=None,
          force=False, threads=8):
    db_dir = Path(db_dir)
    if (db_dir / "giftag.json").exists() and not force:
        raise GiftagError(f"{db_dir} already holds a giftag database; pass --force to rebuild it")
    unknown = set(sources) - set(S.SOURCES)
    if unknown:
        raise GiftagError(f"unknown sources: {', '.join(sorted(unknown))}")

    markers, gifter = load_markers(gifter_db)
    log(f"{len(markers)} markers from gifter database {gifter.get('gifter_db_version', gifter['source'])}")

    ledger = _Ledger()
    wanted = defaultdict(list)
    for marker in markers:
        source = S.NAMESPACE_SOURCE.get(marker.namespace)
        if source is None:
            ledger.refuse(marker, "", "unsupported_namespace",
                          f"no sequence profile source serves the {marker.namespace} namespace")
        elif source not in sources:
            ledger.refuse(marker, source, "source_skipped", "source left out of this build")
        else:
            wanted[source].append(marker)

    partial = db_dir.with_name(db_dir.name + ".partial")
    if partial.exists():
        shutil.rmtree(partial)
    partial.mkdir(parents=True)

    provenance = {}
    builders = {
        "kofam": lambda m: _build_kofam(m, partial, ledger, kofam_release, kofam_dir),
        "ncbifam": lambda m: _build_ncbifam(m, partial, ledger, ncbifam_release, ncbifam_dir, threads),
        "pfam": lambda m: _build_pfam(m, partial, ledger, pfam_dir, threads),
        "dbcan": lambda m: _build_dbcan(m, partial, ledger, dbcan_release, dbcan_dir),
    }
    for source in S.SOURCES:
        if wanted.get(source):
            log(f"building {source} ({len(wanted[source])} markers)")
            provenance[source] = builders[source](wanted[source])
            provenance[source]["terms"] = S.TERMS[source]

    _write_tsv(partial / "profiles.tsv", PROFILE_COLUMNS, ledger.profiles)
    marker_rows = []
    for marker in markers:
        source, status, profile, detail = ledger.status[marker]
        marker_rows.append(dict(namespace=marker.namespace, accession=marker.accession,
                                name=marker.name, source=source, status=status,
                                profile=profile, detail=detail))
    _write_tsv(partial / "markers.tsv", MARKER_COLUMNS, marker_rows)

    counts = Counter(row["status"] for row in marker_rows)
    manifest = {
        "format": DB_FORMAT,
        "giftag_version": __version__,
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "gifter": gifter,
        "sources": provenance,
        "marker_status": dict(sorted(counts.items())),
    }
    with open(partial / "giftag.json", "w") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")

    if db_dir.exists():
        shutil.rmtree(db_dir)
    partial.rename(db_dir)
    summary = ", ".join(f"{n} {status}" for status, n in sorted(counts.items()))
    log(f"database written to {db_dir}: {summary}")
    return manifest


def _write_tsv(path, columns, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _local(directory, *names):
    for name in names:
        path = Path(directory) / name
        if path.exists():
            return path
    raise GiftagError(f"{directory}: none of {', '.join(names)} found")


def _read_hmm(data, where):
    try:
        with pyhmmer.plan7.HMMFile(io.BytesIO(data)) as handle:
            return handle.read()
    except Exception as error:  # pyhmmer raises several parser error types
        raise GiftagError(f"{where}: not a readable HMM ({error})") from None


def _fetch_many(locations, threads):
    """Fetch `{key: location}` concurrently; a missing file maps to None."""

    def fetch(item):
        key, location = item
        try:
            return key, read_all(location, key)
        except GiftagError as error:
            if "404" in str(error) or "no such file" in str(error):
                return key, None
            raise

    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        return dict(pool.map(fetch, locations.items()))


# -- KOfam -------------------------------------------------------------------

def _build_kofam(markers, out, ledger, release, local_dir):
    if local_dir:
        ko_list_location = _local(local_dir, "ko_list", "ko_list.gz")
        profiles_dir = Path(local_dir) / "profiles"
        if profiles_dir.is_dir():
            archive = None
        else:
            profiles_dir, archive = None, _local(local_dir, "profiles.tar.gz")
    else:
        ko_list_location, archive = S.kofam_urls(release)
        profiles_dir = None

    thresholds = {}
    with open_first(ko_list_location, "KOfam ko_list") as stream:
        lines = text_lines(stream)
        next(lines, None)
        for line in lines:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                value = None if fields[1] in ("-", "") else float(fields[1])
                thresholds[fields[0]] = (value, fields[2])
        stream.drain()
        meta = {"release": release if not local_dir else "local",
                "ko_list": stream.location, "ko_list_sha256": stream.sha256,
                "last_modified": stream.last_modified}

    need = set()
    for marker in markers:
        if marker.accession not in thresholds:
            ledger.refuse(marker, "kofam", "missing_from_source", "K number not in ko_list")
        elif thresholds[marker.accession][0] is None:
            ledger.refuse(marker, "kofam", "no_threshold",
                          "KOfam publishes no score threshold for this KO, so KofamScan never assigns it")
        else:
            need.add(marker.accession)

    found = {}
    if profiles_dir is not None:
        for knum in need:
            path = profiles_dir / f"{knum}.hmm"
            if path.exists():
                found[knum] = path.read_bytes()
        meta["profiles"] = str(profiles_dir)
    else:
        with open_first(archive, "KOfam profiles") as stream:
            with tarfile.open(fileobj=stream.reader, mode="r|*") as tar:
                for member in tar:
                    base = os.path.basename(member.name)
                    if member.isfile() and base.endswith(".hmm") and base[:-4] in need:
                        found[base[:-4]] = tar.extractfile(member).read()
            stream.drain()
            meta.update(profiles=stream.location, profiles_sha256=stream.sha256)

    with open(out / "kofam.hmm", "wb") as handle:
        for knum in sorted(found):
            handle.write(found[knum])
    for marker in markers:
        if marker.accession not in need:
            continue
        if marker.accession in found:
            value, score_type = thresholds[marker.accession]
            rule = "kofam_domain" if score_type == "domain" else "kofam_full"
            ledger.searchable(marker, "kofam", marker.accession, rule, value)
        else:
            ledger.refuse(marker, "kofam", "missing_from_source", "no profile in the KOfam archive")
    meta["profiles_kept"] = len(found)
    return meta


# -- NCBIfam and TIGRFAM ---------------------------------------------------------

def _build_ncbifam(markers, out, ledger, release, local_dir, threads):
    table = _local(local_dir, "hmm_PGAP.tsv") if local_dir else S.ncbifam_table_url(release)
    versions = defaultdict(list)
    with Stream(table, "NCBIfam table") as stream:
        for line in text_lines(stream):
            if line.startswith("#") or not line.strip():
                continue
            accession = line.split("\t", 1)[0]
            base, _, version = accession.partition(".")
            versions[base].append((int(version or 0), accession))
        stream.drain()
        meta = {"release": release if not local_dir else "local", "table": str(table),
                "table_sha256": stream.sha256}

    def latest(base):
        return max(versions[base])[1] if versions.get(base) else None

    target = {}
    for marker in markers:
        if marker.namespace == "NCBIFAM":
            base = marker.accession.partition(".")[0]
            if any(acc == marker.accession for _, acc in versions.get(base, ())):
                target[marker] = marker.accession
            else:
                have = latest(base)
                detail = f"release has {have}, not {marker.accession}" if have else "accession not in release"
                ledger.refuse(marker, "ncbifam", "missing_from_source", detail)
        else:
            # gifter records TIGRFAM accessions unversioned, as InterPro does;
            # the profile is the version this NCBIfam release publishes.
            have = latest(marker.accession.upper())
            if have:
                target[marker] = have
            else:
                ledger.refuse(marker, "ncbifam", "missing_from_source", "accession not in release")

    def location(accession):
        if local_dir:
            return Path(local_dir) / "hmm_PGAP.HMM" / f"{accession}.HMM"
        return S.ncbifam_profile_url(release, accession)

    fetched = _fetch_many({acc: location(acc) for acc in set(target.values())}, threads)
    kept = {}
    for accession, data in fetched.items():
        if data is not None:
            kept[accession] = (data, _read_hmm(data, accession))

    with open(out / "ncbifam.hmm", "wb") as handle:
        for accession in sorted(kept):
            data, hmm = kept[accession]
            if hmm.cutoffs.trusted_available():
                handle.write(data)
    for marker, accession in target.items():
        if accession not in kept:
            ledger.refuse(marker, "ncbifam", "missing_from_source", f"no profile file for {accession}")
            continue
        hmm = kept[accession][1]
        if not hmm.cutoffs.trusted_available():
            ledger.refuse(marker, "ncbifam", "no_threshold", f"{accession} has no trusted cutoff", accession)
            continue
        ledger.searchable(marker, "ncbifam", accession, "trusted_cutoff", hmm.cutoffs.trusted[0])
    meta["profiles_kept"] = sum(1 for _, hmm in kept.values() if hmm.cutoffs.trusted_available())
    return meta


# -- Pfam ---------------------------------------------------------------------

def _build_pfam(markers, out, ledger, local_dir, threads):
    meta = {"release": "local" if local_dir else _interpro_pfam_release()}

    def location(accession):
        if local_dir:
            return _local(local_dir, f"{accession}.hmm", f"{accession}.hmm.gz")
        return S.PFAM_ENTRY_URL.format(accession=accession)

    bases = {marker.accession.upper().partition(".")[0] for marker in markers}
    fetched = _fetch_many({base: location(base) for base in bases}, threads)
    kept = {}
    with open(out / "pfam.hmm", "wb") as handle:
        for base in sorted(fetched):
            data = fetched[base]
            if data is None:
                continue
            hmm = _read_hmm(data, base)
            if hmm.cutoffs.gathering_available():
                handle.write(data)
            kept[base] = hmm
    for marker in markers:
        hmm = kept.get(marker.accession.upper().partition(".")[0])
        if hmm is None:
            ledger.refuse(marker, "pfam", "missing_from_source", "no Pfam profile for this accession")
        elif not hmm.cutoffs.gathering_available():
            ledger.refuse(marker, "pfam", "no_threshold", "profile has no gathering cutoff", hmm.accession)
        else:
            ledger.searchable(marker, "pfam", hmm.accession, "gathering_cutoff", hmm.cutoffs.gathering[0])
    meta["profiles_kept"] = sum(1 for hmm in kept.values() if hmm.cutoffs.gathering_available())
    return meta


def _interpro_pfam_release():
    try:
        info = json.loads(read_all(S.INTERPRO_API_URL, "InterPro release"))
        pfam = info["databases"]["pfam"]
        return f"Pfam {pfam['version']} ({pfam['releaseDate'][:10]})"
    except Exception:
        return "unknown"


# -- dbCAN --------------------------------------------------------------------

def dbcan_family_key(name):
    """Map a dbCAN family profile name to the CAZy accession it reports."""
    key = name[:-4] if name.endswith(".hmm") else name
    # run_dbcan reports every GT2_* sub-model as plain GT2.
    return "GT2" if key.startswith("GT2_") else key


def dbcan_cluster_key(name):
    """Map a dbCAN-sub profile name (`GH5_e12.hmm|GH5:41|3.2.1.4:9`) to its cluster."""
    return name.split("|", 1)[0].split(".hmm", 1)[0]


def cluster_family(cluster):
    return cluster.rsplit("_e", 1)[0]


def _build_dbcan(markers, out, ledger, release, local_dir):
    if local_dir:
        family_location = _local(local_dir, "dbCAN.hmm")
        sub_location = _local(local_dir, "dbCAN_sub.hmm") if any(
            "_e" in m.accession for m in markers) else None
    else:
        family_location, sub_location = S.dbcan_urls(release)
    clusters = [m for m in markers if "_e" in m.accession]
    families = [m for m in markers if "_e" not in m.accession]
    meta = {"release": release if not local_dir else "local",
            "evalue": S.DBCAN_EVALUE, "coverage": S.DBCAN_COVERAGE, "overlap": S.DBCAN_OVERLAP}

    family_keys = set()
    z_family = 0
    with Stream(family_location, "dbCAN families") as stream, open(out / "dbcan.hmm", "wb") as handle:
        for name, record in iter_hmm_records(stream.reader):
            z_family += 1
            family_keys.add(dbcan_family_key(name))
            handle.write(record)
        meta.update(family_library=str(family_location), family_sha256=stream.sha256,
                    z_family=z_family)

    for marker in families:
        if marker.accession in family_keys:
            ledger.searchable(marker, "dbcan", marker.accession, "dbcan_family", S.DBCAN_EVALUE)
        else:
            ledger.refuse(marker, "dbcan", "missing_from_source", "no dbCAN family profile")

    if clusters:
        parents = {cluster_family(m.accession) for m in clusters}
        seen = set()
        z_sub = 0
        kept = 0
        with Stream(sub_location, "dbCAN-sub") as stream, open(out / "dbcan_sub.hmm", "wb") as handle:
            for name, record in iter_hmm_records(stream.reader):
                z_sub += 1
                cluster = dbcan_cluster_key(name)
                if cluster_family(cluster) in parents:
                    seen.add(cluster)
                    handle.write(record)
                    kept += 1
            meta.update(sub_library=str(sub_location), sub_sha256=stream.sha256,
                        z_sub=z_sub, sub_profiles_kept=kept)
        for marker in clusters:
            if marker.accession not in seen:
                ledger.refuse(marker, "dbcan", "missing_from_source", "no dbCAN-sub profile")
            elif cluster_family(marker.accession) not in family_keys:
                ledger.refuse(marker, "dbcan", "missing_from_source",
                              "parent family has no dbCAN family profile to gate the subfamily search")
            else:
                ledger.searchable(marker, "dbcan_sub", marker.accession, "dbcan_sub", S.DBCAN_EVALUE)
    return meta
