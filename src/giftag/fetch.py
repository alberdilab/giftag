"""Streaming access to local files and URLs, with checksums and progress.

Upstream libraries are large (dbCAN-sub alone is 4.9 GB) and giftag keeps a
small fraction of each, so everything is read as a stream and filtered on the
way through. Nothing is written to disk that is not kept.
"""

import gzip
import hashlib
import io
import os
import threading
import time
import urllib.error
import urllib.request

from rich.markup import escape

from giftag import GiftagError, __version__, ui

_USER_AGENT = f"giftag/{__version__} (+https://github.com/alberdilab/giftag)"


def is_url(location):
    return str(location).startswith(("http://", "https://", "ftp://"))


def log(message):
    ui.info(message)


class _Tally(io.RawIOBase):
    """Counts and hashes every byte read through it, and reports progress."""

    def __init__(self, raw, task):
        self._raw = raw
        self._task = task
        self.bytes = 0
        self.sha256 = hashlib.sha256()

    def readable(self):
        return True

    def readinto(self, buffer):
        n = self._raw.readinto(buffer)
        if n:
            self.sha256.update(memoryview(buffer)[:n])
            self.bytes += n
            if self._task is not None:
                self._task.advance(n)
        return n

    def close(self):
        try:
            self._raw.close()
            if self._task is not None:
                self._task.close(failed=self._task.total is not None
                                 and self.bytes < self._task.total)
        finally:
            super().close()


class Stream:
    """A buffered binary reader over a file or URL that remembers its checksum."""

    def __init__(self, location, label=None, progress=True):
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
        task = None
        if progress:
            verb = "downloaded" if is_url(self.location) else "read"
            task = ui.Task(label, total=total, kind="bytes",
                           done=f"{label}: {verb} {{amount}} in {{elapsed}}")
        self._tally = _Tally(raw, task)
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


def open_first(locations, label=None, progress=True):
    """Open the first of several mirrors that answers."""
    if isinstance(locations, (str, os.PathLike)):
        return Stream(locations, label, progress)
    errors = []
    for location in locations:
        try:
            return Stream(location, label, progress)
        except GiftagError as error:
            errors.append(str(error))
            ui.warning(f"{escape(str(error))}; trying the next mirror")
    raise GiftagError("; ".join(errors))


def read_all(location, label=None):
    """Return the full content of a file or URL, gunzipped if it is gzip."""
    with Stream(location, label, progress=False) as stream:
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


class RangeFile:
    """Random access to a local file or an HTTP URL that serves byte ranges."""

    def __init__(self, location):
        self.location = str(location)
        self.bytes_read = 0
        self._lock = threading.Lock()
        if is_url(self.location):
            request = urllib.request.Request(self.location, method="HEAD",
                                             headers={"User-Agent": _USER_AGENT})
            try:
                with urllib.request.urlopen(request, timeout=120) as response:
                    headers = response.headers
            except (urllib.error.URLError, OSError) as error:
                raise GiftagError(f"{self.location}: {getattr(error, 'reason', error)}") from None
            if headers.get("Accept-Ranges") != "bytes":
                raise GiftagError(f"{self.location}: server does not serve byte ranges")
            self.size = int(headers["Content-Length"])
            self.etag = (headers.get("ETag") or "").strip('"')
            self.last_modified = headers.get("Last-Modified")
        else:
            if not os.path.exists(self.location):
                raise GiftagError(f"{self.location}: no such file")
            self.size = os.path.getsize(self.location)
            self.etag = None
            self.last_modified = None

    def read(self, start, end):
        """Bytes `[start, end)`, clipped to the file."""
        end = min(end, self.size)
        if start >= end:
            return b""
        if not is_url(self.location):
            with open(self.location, "rb") as handle:
                handle.seek(start)
                data = handle.read(end - start)
        else:
            request = urllib.request.Request(
                self.location, headers={"User-Agent": _USER_AGENT, "Range": f"bytes={start}-{end - 1}"})
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(request, timeout=120) as response:
                        if response.status != 206:
                            raise GiftagError(f"{self.location}: range request returned {response.status}")
                        data = response.read()
                    break
                except (urllib.error.URLError, OSError) as error:
                    if attempt == 2:
                        raise GiftagError(f"{self.location}: {getattr(error, 'reason', error)}") from None
                    time.sleep(2 ** attempt)
        if len(data) != end - start:
            raise GiftagError(f"{self.location}: short read at byte {start}")
        with self._lock:
            self.bytes_read += len(data)
        return data
