"""Read the marker list giftag is asked to search for.

The list comes from gifter itself, so giftag never searches for a marker no
GIFT uses and never misses one a GIFT needs. Either gifter's compiled SQLite
database (the `marker` table) or a TSV with `namespace` and `accession` columns
is accepted.
"""

import csv
import os
import sqlite3
import tempfile
from dataclasses import dataclass

from giftag import GiftagError
from giftag.fetch import is_url, read_all


@dataclass(frozen=True)
class Marker:
    namespace: str
    accession: str
    name: str = ""


def load_markers(location):
    """Return `(markers, provenance)` from a gifter SQLite file or a TSV."""
    location = str(location)
    provenance = {"source": location}
    if is_url(location):
        data = read_all(location, "gifter database")
        with tempfile.NamedTemporaryFile(suffix=".sqlite", delete=False) as handle:
            handle.write(data)
            path = handle.name
        try:
            markers, release = _from_sqlite(path)
        finally:
            os.unlink(path)
    elif location.endswith((".tsv", ".txt")):
        markers, release = _from_tsv(location), {}
    else:
        if not os.path.exists(location):
            raise GiftagError(f"{location}: no such file")
        markers, release = _from_sqlite(location)
    provenance.update(release)
    unique = sorted(set(markers), key=lambda m: (m.namespace, m.accession))
    if not unique:
        raise GiftagError(f"{location}: no markers found")
    provenance["markers"] = len(unique)
    return unique, provenance


def _from_sqlite(path):
    try:
        connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    except sqlite3.Error as error:
        raise GiftagError(f"{path}: {error}") from None
    try:
        rows = connection.execute("SELECT namespace, accession, name FROM marker").fetchall()
        release = {}
        try:
            row = connection.execute(
                "SELECT gifter_db_version, schema_version, build_date, source_commit "
                "FROM database_release"
            ).fetchone()
            if row:
                release = dict(
                    zip(("gifter_db_version", "schema_version", "build_date", "source_commit"), row)
                )
        except sqlite3.Error:
            pass
    except sqlite3.Error as error:
        raise GiftagError(f"{path}: not a gifter database ({error})") from None
    finally:
        connection.close()
    return [Marker(ns.upper(), acc, name or "") for ns, acc, name in rows], release


def _from_tsv(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"namespace", "accession"} - set(reader.fieldnames or ())
        if missing:
            raise GiftagError(f"{path}: missing columns {', '.join(sorted(missing))}")
        return [
            Marker(row["namespace"].strip().upper(), row["accession"].strip(), row.get("name") or "")
            for row in reader
            if row["namespace"].strip() and row["accession"].strip()
        ]
