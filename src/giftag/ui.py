"""Operational output: the package logger, the console, and progress tasks.

Everything giftag says while it works goes through here, so the library code
never decides how it is shown. On an interactive terminal a task is a live
Rich progress bar; in a job log or a pipe it is a handful of plain log lines;
with `--quiet`, or when giftag is used as a library and nothing was
configured, it is silent.
"""

import logging
import threading
import time

from rich.console import Console
from rich.highlighter import NullHighlighter
from rich.logging import RichHandler
from rich.progress import (
    BarColumn,
    Progress,
    ProgressColumn,
    SpinnerColumn,
    TextColumn,
    TimeRemainingColumn,
)
from rich.text import Text

logger = logging.getLogger("giftag")
logger.addHandler(logging.NullHandler())
console = Console(stderr=True, highlight=False)

_STARTED = time.monotonic()
_lock = threading.Lock()
_progress = None
_open_tasks = 0
_live = False


def format_duration(seconds):
    """`42.0s`, `3m09s`, `1h02m09s`."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, rest = divmod(int(seconds), 60)
    if minutes < 60:
        return f"{minutes}m{rest:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h{minutes:02d}m{rest:02d}s"


def format_bytes(value):
    size = float(value)
    for unit in ("B", "kB", "MB", "GB"):
        if size < 1000 or unit == "GB":
            return f"{size:,.0f} {unit}" if unit in ("B", "kB") else f"{size:,.1f} {unit}"
        size /= 1000


class _Elapsed(logging.Formatter):
    def format(self, record):
        stamp = format_duration(time.monotonic() - _STARTED)
        return f"[dim]\\[{stamp:>9}][/dim] {super().format(record)}"


class _Utc(logging.Formatter):
    converter = time.gmtime


class _PlainUtc(_Utc):
    """The file log carries no Rich markup."""

    def format(self, record):
        return Text.from_markup(super().format(record)).plain


def configure(verbose=False, quiet=False, log_file=None):
    """Attach the console handler, and a UTC-stamped file handler if asked."""
    global _live
    for handler in list(logger.handlers):
        if not isinstance(handler, logging.NullHandler):
            logger.removeHandler(handler)
            handler.close()
    handler = RichHandler(console=console, show_path=False, show_time=False,
                          show_level=False, markup=True, rich_tracebacks=False,
                          highlighter=NullHighlighter())
    handler.setFormatter(_Elapsed("%(message)s"))
    handler.setLevel(logging.WARNING if quiet else logging.DEBUG if verbose else logging.INFO)
    logger.addHandler(handler)
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
        file_handler.setFormatter(_PlainUtc("%(asctime)s %(levelname)-7s %(message)s",
                                            datefmt="%Y-%m-%dT%H:%M:%SZ"))
        file_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
        logger.addHandler(file_handler)
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    _live = console.is_terminal and not quiet


def info(message):
    logger.info(message)


def warning(message):
    logger.warning(f"[yellow]{message}[/yellow]")


class _Amount(ProgressColumn):
    """`312.4 MB / 1.6 GB · 1.1 MB/s` for byte tasks, `41 / 147` for counts."""

    def render(self, task):
        if task.fields.get("kind") == "bytes":
            text = format_bytes(task.completed)
            if task.total:
                text += f" / {format_bytes(task.total)}"
            if task.speed:
                text += f" · {format_bytes(task.speed)}/s"
        elif task.total:
            text = f"{int(task.completed):,} / {int(task.total):,}"
        else:
            text = f"{int(task.completed):,}" if task.completed else ""
        return Text(text, style="progress.download")


def _acquire():
    global _progress, _open_tasks
    with _lock:
        if _progress is None:
            _progress = Progress(
                SpinnerColumn(style="cyan"),
                TextColumn("{task.description}"),
                BarColumn(bar_width=28),
                _Amount(),
                TimeRemainingColumn(compact=True),
                console=console, transient=True, refresh_per_second=8,
            )
            _progress.start()
        _open_tasks += 1
        return _progress


def _release():
    global _progress, _open_tasks
    with _lock:
        _open_tasks -= 1
        if _open_tasks <= 0 and _progress is not None:
            _progress.stop()
            _progress = None
            _open_tasks = 0


class Task:
    """One unit of long work. Use as a context manager, or call `close()`.

    `kind` is `bytes` or `count`. `done` is the line logged when the task ends;
    it may use `{amount}` and `{elapsed}`. Off a terminal the task logs a line
    every tenth of the way instead of drawing a bar, unless `steps` is false.
    """

    def __init__(self, description, total=None, kind="count", done=None, steps=True):
        self.description = description
        self.total = total or None
        self.kind = kind
        self.completed = 0
        self._done = done
        self._started = time.monotonic()
        self._closed = False
        self._lock = threading.Lock()
        if _live:
            self._bar = _acquire()
            self._id = self._bar.add_task(description, total=self.total, kind=kind)
        else:
            self._bar = None
            # Off a terminal: a line every tenth of the way, or every 250 MB.
            # `steps=False` is for a task whose caller logs each step itself.
            self._step = 0 if not steps else (
                (self.total / 10) if self.total else (250e6 if kind == "bytes" else 0))
            self._next = self._step

    def advance(self, amount=1):
        with self._lock:
            self.completed += amount
            if self._bar is not None:
                self._bar.update(self._id, advance=amount)
            elif self._step and self.completed >= self._next and not self._finished():
                self._next = self.completed + self._step
                logger.info(f"{self.description}: {self._amount()}")

    def _finished(self):
        return self.total is not None and self.completed >= self.total

    def update(self, description=None, total=None):
        if description is not None:
            self.description = description
        if total is not None:
            self.total = total
        if self._bar is not None:
            self._bar.update(self._id, description=self.description, total=self.total)

    def _amount(self):
        done = format_bytes(self.completed) if self.kind == "bytes" else f"{self.completed:,}"
        if not self.total:
            return done
        total = format_bytes(self.total) if self.kind == "bytes" else f"{self.total:,}"
        return f"{done} / {total}"

    def close(self, failed=False):
        if self._closed:
            return
        self._closed = True
        if self._bar is not None:
            self._bar.remove_task(self._id)
            _release()
        if self._done and not failed:
            amount = format_bytes(self.completed) if self.kind == "bytes" else f"{self.completed:,}"
            logger.info(self._done.format(
                amount=amount, elapsed=format_duration(time.monotonic() - self._started)))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, *exc):
        self.close(failed=exc_type is not None)
