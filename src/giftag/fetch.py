"""Streaming access to local files and URLs, with checksums and progress.

Upstream libraries are large (dbCAN-sub alone is 4.9 GB) and giftag keeps a
small fraction of each, so everything is read as a stream and filtered on the
way through. Nothing is written to disk that is not kept.
"""

import gzip
import hashlib
import io
import os
import sys
import time
import urllib.error
import urllib.request

from giftag import GiftagError, __version__

_USER_AGENT = f"giftag/{__version__} (+https://github.com/alberdilab/giftag)"


def is_url(location):
    return str(location).startswith(("http://", "https://", "ftp://"))


def log(message):
    print(f"giftag: {message}", file=sys.stderr, flush=True)


class _Tally(io.RawIOBase):
    """Counts and hashes every byte read through it, and reports progress."""

    def __init__(self, raw, label, total):
        self._raw = raw
        self.label = label
        self.total = total
        self.bytes = 0
        self.sha256 = hashlib.sha256()
        self._next_report = self._step()
        self._started = time.monotonic()

    def readable(self):
        return True

    def readinto(self, buffer):
        n = self._raw.readinto(buffer)
        if n:
            self.sha256.update(memoryview(buffer)[:n])
            self.bytes += n
            if self.bytes >= self._next_report:
                self._report()
        return n

    def _report(self):
        mb = self.bytes / 1e6
        if self.total:
            log(f"{self.label}: {mb:,.0f} / {self.total / 1e6:,.0f} MB")
        elif mb >= 1:
            log(f"{self.label}: {mb:,.0f} MB")
        self._next_report = self.bytes + self._step()

    def _step(self):
        return max(50_000_000, (self.total or 0) // 20)

    def close(self):
        try:
            self._raw.close()
        finally:
            super().close()


class Stream:
    """A buffered binary reader over a file or URL that remembers its checksum."""

    def __init__(self, location, label=None):
        self.location = str(location)
        label = label or os.path.basename(self.location.split("?")[0])
        if is_url(self.location):
            request = urllib.request.Request(self.location, headers={"User-Agent": _USER_AGENT})
            try:
                raw = urllib.request.urlopen(request, timeout=120)
            except urllib.error.HTTPError as error:
                raise GiftagError(f"{self.location}: HTTP {error.code} {error.reason}") from None
            except (urllib.error.URLError, OSError) as error:
                raise GiftagError(f"{self.location}: {getattr(error, 'reason', error)}") from None
            total = int(raw.headers.get("Content-Length") or 0)
            self.last_modified = raw.headers.get("Last-Modified")
        else:
            if not os.path.exists(self.location):
                raise GiftagError(f"{self.location}: no such file")
            raw = open(self.location, "rb")
            total = os.path.getsize(self.location)
            self.last_modified = None
        self._tally = _Tally(raw, label, total)
        self.reader = io.BufferedReader(self._tally, buffer_size=1 << 20)

    @property
    def sha256(self):
        return self._tally.sha256.hexdigest()

    @property
    def bytes_read(self):
        return self._tally.bytes

    def drain(self):
        """Read to the end, so the checksum covers the whole file."""
        while self.reader.read(1 << 20):
            pass

    def close(self):
        self.reader.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def open_first(locations, label=None):
    """Open the first of several mirrors that answers."""
    if isinstance(locations, (str, os.PathLike)):
        return Stream(locations, label)
    errors = []
    for location in locations:
        try:
            return Stream(location, label)
        except GiftagError as error:
            errors.append(str(error))
            log(f"{error}; trying the next mirror")
    raise GiftagError("; ".join(errors))


def read_all(location, label=None):
    """Return the full content of a file or URL, gunzipped if it is gzip."""
    with Stream(location, label) as stream:
        data = stream.reader.read()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data


def text_lines(stream):
    """Iterate decoded lines of a stream, gunzipping it if needed."""
    reader = stream.reader
    if reader.peek(2)[:2] == b"\x1f\x8b":
        reader = gzip.GzipFile(fileobj=reader)
    for line in reader:
        yield line.decode("utf-8", "replace")


def iter_hmm_records(reader):
    """Yield `(name, record_bytes)` for each profile in a HMMER3 text file."""
    lines = []
    name = None
    for line in reader:
        lines.append(line)
        if name is None and line.startswith(b"NAME"):
            name = line.split(None, 1)[1].strip().decode()
        elif line.startswith(b"//"):
            yield name, b"".join(lines)
            lines = []
            name = None
    if any(line.strip() for line in lines):
        raise GiftagError("HMM file ended inside a profile; the download was probably truncated")
