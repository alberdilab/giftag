"""Run one benchmark command and retain reproducible resource provenance.

Run each invocation in a fresh process. RUSAGE_CHILDREN peak RSS is the
largest child-process high-water mark, not concurrent workers' total memory.
"""

import argparse
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tool", required=True)
    parser.add_argument("--scope", choices=("full", "subset"), required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--input", type=Path, action="append", default=[])
    parser.add_argument("--database", type=Path, action="append", default=[])
    parser.add_argument("--version", default="", help="exact tool version or release")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("a command is required after --")
    if args.outdir.exists() and any(args.outdir.iterdir()):
        parser.error(f"{args.outdir} is nonempty; use a new output directory")
    args.outdir.mkdir(parents=True, exist_ok=True)
    inputs = []
    for path in args.input:
        if not path.is_file():
            parser.error(f"missing input: {path}")
        inputs.append({"path": str(path.resolve()), "bytes": path.stat().st_size,
                       "sha256": sha256(path)})
    databases = []
    for path in args.database:
        if not path.exists():
            parser.error(f"missing database: {path}")
        databases.append({"path": str(path.resolve()),
                          "bytes": path.stat().st_size if path.is_file() else None,
                          "sha256": sha256(path) if path.is_file() else None})
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = datetime.now(timezone.utc).isoformat()
    began = time.perf_counter()
    with open(args.outdir / "stdout.log", "wb") as stdout, \
            open(args.outdir / "stderr.log", "wb") as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    wall = time.perf_counter() - began
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    # Linux reports ru_maxrss in KiB; macOS reports bytes.
    peak_bytes = after.ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    record = {
        "tool": args.tool, "scope": args.scope, "version": args.version,
        "threads": args.threads, "command": command, "cwd": os.getcwd(),
        "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "exit_code": result.returncode, "wall_seconds": wall,
        "user_cpu_seconds": after.ru_utime - before.ru_utime,
        "system_cpu_seconds": after.ru_stime - before.ru_stime,
        "peak_child_rss_bytes": peak_bytes, "hostname": platform.node(),
        "platform": platform.platform(), "python": platform.python_version(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
        "inputs": inputs, "databases": databases,
    }
    (args.outdir / "measure.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({k: record[k] for k in ("tool", "scope", "exit_code",
                     "wall_seconds", "peak_child_rss_bytes")}, sort_keys=True))
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
