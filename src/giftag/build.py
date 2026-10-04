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
import re
import shutil
import tarfile
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from importlib import resources
from pathlib import Path

import pyhmmer

from rich.markup import escape

from giftag import DB_FORMAT, GiftagError, __version__, ui
from giftag import sources as S
from giftag.fetch import RangeFile, Stream, iter_hmm_records, open_first, read_all, text_lines
from giftag.markers import load_markers

_TITLES = {"kofam": "KOfam", "ncbifam": "NCBIfam", "pfam": "Pfam", "dbcan": "dbCAN"}
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
    version = gifter.get("gifter_db_version")
    ui.info(f"gifter database {version or escape(Path(gifter['source']).name)}: "
            f"{len(markers):,} markers")

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
            ui.info(f"[bold cyan]{_TITLES[source]}[/bold cyan] · {len(wanted[source]):,} marker{'s' if len(wanted[source]) != 1 else ''}")
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


def _fetch_many(locations, threads, label):
    """Fetch `{key: location}` concurrently; a missing file maps to None."""
    with ui.Task(label, total=len(locations),
                 done=f"{label}: fetched {{amount}} in {{elapsed}}") as task:

        def fetch(item):
            key, location = item
            try:
                return key, read_all(location, key)
            except GiftagError as error:
                if "404" in str(error) or "no such file" in str(error):
                    return key, None
                raise
            finally:
                task.advance()

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
    with open_first(ko_list_location, "KOfam ko_list", progress=False) as stream:
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

    fetched = _fetch_many({acc: location(acc) for acc in set(target.values())}, threads,
                          "NCBIfam profiles")
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
    fetched = _fetch_many({base: location(base) for base in bases}, threads, "Pfam profiles")
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


def cazy_family(key):
    """The CAZy family of a family-library key: `GH43_18` and `GH43` are both GH43.

    dbCAN's family library holds models for official CAZy subfamilies beside
    the family models, and run_dbcan reports whichever wins a region."""
    return re.sub(r"_\d+$", "", key)


def _build_dbcan(markers, out, ledger, release, local_dir):
    if local_dir:
        family_location = _local(local_dir, "dbCAN.hmm")
        sub_location = _local(local_dir, "dbCAN_sub.hmm", "dbCAN-sub.hmm") if any(
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

    family_keys |= {cazy_family(key) for key in family_keys}
    for marker in families:
        if marker.accession in family_keys:
            ledger.searchable(marker, "dbcan", marker.accession, "dbcan_family", S.DBCAN_EVALUE)
        else:
            ledger.refuse(marker, "dbcan", "missing_from_source", "no dbCAN family profile")

    if clusters:
        parents = {cluster_family(m.accession) for m in clusters}
        layout = None if local_dir else dbcan_sub_layout(release)
        records = None
        if layout is not None:
            try:
                records, ranged = _fetch_sub_families(RangeFile(sub_location), parents, layout)
                meta.update(sub_library=sub_location, sub_etag=ranged["etag"],
                            sub_size=ranged["size"], sub_bytes_fetched=ranged["bytes_read"],
                            z_sub=sum(layout.values()))
            except GiftagError as error:
                ui.warning(f"dbCAN-sub ranged download failed ({escape(str(error))}); "
                           "streaming the whole library")
        if records is None:
            records = []
            z_sub = 0
            with Stream(sub_location, "dbCAN-sub") as stream:
                for name, record in iter_hmm_records(stream.reader):
                    z_sub += 1
                    if cluster_family(dbcan_cluster_key(name)) in parents:
                        records.append((name, record))
                meta.update(sub_library=str(sub_location), sub_sha256=stream.sha256, z_sub=z_sub)
        seen = set()
        with open(out / "dbcan_sub.hmm", "wb") as handle:
            for name, record in records:
                seen.add(dbcan_cluster_key(name))
                handle.write(record)
        meta["sub_profiles_kept"] = len(records)
        for marker in clusters:
            if marker.accession not in seen:
                ledger.refuse(marker, "dbcan", "missing_from_source", "no dbCAN-sub profile")
            elif cluster_family(marker.accession) not in family_keys:
                ledger.refuse(marker, "dbcan", "missing_from_source",
                              "parent family has no dbCAN family profile to gate the subfamily search")
            else:
                ledger.searchable(marker, "dbcan_sub", marker.accession, "dbcan_sub", S.DBCAN_EVALUE)
    return meta


# -- dbCAN-sub by byte range ---------------------------------------------------
#
# dbCAN_sub.hmm is 5 GB, and gifter needs the clusters of a few dozen of its
# 500 families. The file lists profiles in natural order of their names (AA1,
# AA3, AA11, ... CBM, CE, GH, GT, PL), so each family is one contiguous block
# that a binary search over byte offsets can find. For a release giftag knows,
# a packaged table gives each family's profile count, which checks every
# fetched block and supplies run_dbcan's Z without reading the whole file. Any
# inconsistency falls back to streaming the full library.

_PROBE_SPAN = 1 << 16
_SCAN_SPAN = 1 << 18
_BLOCK_CHUNK = 1 << 23
_COARSE_PROBES = 256
_PROBE_THREADS = 16


def dbcan_sub_layout(release):
    """`{family: profile count}` for a dbCAN release giftag knows, else None."""
    try:
        text = resources.files("giftag").joinpath("data", f"dbcan_sub_{release}.tsv").read_text()
    except (FileNotFoundError, OSError):
        return None
    layout = {}
    for line in text.splitlines():
        if line and not line.startswith(("#", "family\t")):
            family, count = line.split("\t")
            layout[family] = int(count)
    return layout


def natural_key(name):
    """Sort key matching dbCAN-sub's order: numbers compare as numbers."""
    return tuple((0, int(part), "") if part.isdigit() else (1, 0, part)
                 for part in re.split(r"(\d+)", name) if part)


def _headers(data, base):
    """`(start, cluster)` for every profile that starts inside `data`, which
    begins at byte `base` (`base` 0 is the file start)."""
    found = []
    # The NAME line must be complete: a read that ends inside it would
    # otherwise yield a truncated name and a wrong sort key.
    for match in re.finditer(rb"(?:^|\n)HMMER3/[^\n]*\nNAME +(\S+)[^\n]*\n", data):
        start = base + match.start() + (0 if match.group(0).startswith(b"HMMER3") else 1)
        if start > 0 or base == 0:
            found.append((start, dbcan_cluster_key(match.group(1).decode())))
    return found


def _record_at_or_after(rangefile, offset, cache):
    """`(start, cluster)` of the first profile starting at or after `offset`,
    or `(size, None)` past the last one."""
    if offset in cache:
        return cache[offset]
    base = max(0, offset - 1)
    data = b""
    while True:
        end = base + len(data) + _PROBE_SPAN
        data += rangefile.read(base + len(data), end)
        hits = [h for h in _headers(data, base) if h[0] >= offset]
        if hits:
            cache[offset] = hits[0]
            return hits[0]
        if end >= rangefile.size:
            cache[offset] = (rangefile.size, None)
            return cache[offset]


def _first_at_or_above(rangefile, target, cache, lo=0, hi=None):
    """Byte offset of the first profile whose name sorts at or above `target`.

    Every profile starting before `lo` sorts below `target`, and the first
    profile starting at or after `hi` sorts at or above it. Bisect until the
    gap is small, then read it whole and scan its headers."""
    hi = rangefile.size if hi is None else hi
    while hi - lo > _SCAN_SPAN:
        mid = (lo + hi) // 2
        start, cluster = _record_at_or_after(rangefile, mid, cache)
        if cluster is None or natural_key(cluster) >= target:
            hi = mid
        else:
            lo = start + 1
    candidate = _record_at_or_after(rangefile, hi, cache)
    base = max(0, lo - 1)
    data = rangefile.read(base, max(candidate[0], lo) + 4096)
    for start, cluster in _headers(data, base):
        if lo <= start < candidate[0] and natural_key(cluster) >= target:
            return start
    return candidate[0]


def _bracket(samples, target, size):
    """Offsets `[lo, hi]` known to contain the boundary for `target`, from a
    coarse map of `(probe_offset, start, cluster)` samples in file order."""
    lo, hi = 0, size
    for offset, start, cluster in samples:
        if cluster is not None and natural_key(cluster) < target:
            lo = start + 1
        else:
            hi = offset
            break
    return lo, hi


def _fetch_sub_families(rangefile, families, layout):
    """Return the `(name, record)` profiles of `families`, in file order."""
    cache = {}
    wanted = sorted((f for f in families if f in layout), key=natural_key)
    # A coarse map of the file, then a short bisection inside each bracket.
    # Every request is latency-bound, so both steps run in parallel.
    with ThreadPoolExecutor(max_workers=_PROBE_THREADS) as pool, ui.Task(
            "dbCAN-sub · locating families", total=_COARSE_PROBES + 2 * len(wanted)) as task:

        def stepped(function):
            def run(argument):
                try:
                    return function(argument)
                finally:
                    task.advance()
            return run

        offsets = [i * rangefile.size // _COARSE_PROBES for i in range(_COARSE_PROBES)]
        probe = stepped(lambda o: _record_at_or_after(rangefile, o, cache))
        samples = [(o, *pair) for o, pair in zip(offsets, pool.map(probe, offsets))]
        targets = []
        for family in wanted:
            lower = natural_key(f"{family}_e")
            targets += [lower, lower + ((2, 0, ""),)]
        bounds = list(pool.map(stepped(
            lambda t: _first_at_or_above(rangefile, t, cache, *_bracket(samples, t, rangefile.size))),
            targets))
    blocks = [(bounds[2 * i], bounds[2 * i + 1], family) for i, family in enumerate(wanted)]
    probed = rangefile.bytes_read
    records = []
    with ui.Task("dbCAN-sub", total=sum(end - start for start, end, _ in blocks), kind="bytes",
                 done="dbCAN-sub: downloaded {amount} in {elapsed}") as task:
        for start, end, family in sorted(blocks):
            data = bytearray()
            for chunk in range(start, end, _BLOCK_CHUNK):
                data += rangefile.read(chunk, min(end, chunk + _BLOCK_CHUNK))
                task.advance(min(end, chunk + _BLOCK_CHUNK) - chunk)
            block = list(iter_hmm_records(io.BytesIO(bytes(data))))
            wrong = [name for name, _ in block if cluster_family(dbcan_cluster_key(name)) != family]
            if wrong or len(block) != layout[family]:
                raise GiftagError(
                    f"{family}: found {len(block)} profiles where {layout[family]} were expected"
                    + (f", including {wrong[0]}" if wrong else ""))
            records.extend(block)
    ui.info(f"dbCAN-sub: {len(records):,} profiles of {len(blocks)} families, "
            f"{ui.format_bytes(rangefile.bytes_read)} of a {ui.format_bytes(rangefile.size)} library "
            f"({ui.format_bytes(probed)} to locate them)")
    return records, {"etag": rangefile.etag, "size": rangefile.size,
                     "bytes_read": rangefile.bytes_read}
