"""Repeat complete production reads sequentially in fresh processes.

Pass --baseline LABEL ROOT for each older checkout and --python PATH for each
interpreter. See optimization.md for the measured corpus and reproduction.
"""

import argparse
import hashlib
import json
import os
import statistics
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def summarize(records):
    groups = {}
    for record in records:
        key = (
            record["versions"]["pandas"],
            record["path"],
            record["release"],
            tuple(sorted(record["versions"].items())),
        )
        groups.setdefault(key, []).append(record)
    summaries = []
    for (pandas, path, release, _versions), runs in groups.items():
        summary = {"pandas": pandas, "input": path, "release": release, "runs": len(runs)}
        summary["versions"] = runs[0]["versions"]
        for name in ("wall_seconds", "cpu_seconds", "peak_rss_mib"):
            values = [run[name] for run in runs]
            summary[name] = {
                "median": statistics.median(values),
                "min": min(values),
                "max": max(values),
            }
        summary["stage_medians"] = {
            name: statistics.median(run["stages"][name] for run in runs)
            for name in runs[0]["stages"]
        }
        summary["checksum"] = runs[0]["checksum"]
        summaries.append(summary)
    return summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument(
        "--baseline", nargs=2, action="append", default=[], metavar=("LABEL", "ROOT")
    )
    parser.add_argument("--python", action="append", dest="interpreters")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--result-type", choices=["pandas", "polars"], default="pandas")
    parser.add_argument("--features", nargs="+")
    parser.add_argument("--usecols", nargs="+")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.repeats < 1 or args.threads < 1:
        parser.error("repeats and threads must be positive")
    if len({path.name for path in args.inputs}) != len(args.inputs):
        parser.error("input basenames must be distinct")
    script = Path(__file__).with_name("compare_readers.py")
    releases = [*args.baseline, ("optimized", str(script.parent.parent))]
    expected = {}
    records = []
    environment = dict(os.environ, POLARS_MAX_THREADS=str(args.threads))
    with args.output.with_suffix(".jsonl").open("w") as raw:
        for interpreter in args.interpreters or [sys.executable]:
            for repeat in range(1, args.repeats + 1):
                # Rotate reader order and reverse corpus order to distribute ordering effects.
                offset = (repeat - 1) % len(releases)
                order = releases[offset:] + releases[:offset]
                paths = args.inputs if repeat % 2 else args.inputs[::-1]
                for path in paths:
                    for label, root in order:
                        command = [
                            interpreter,
                            str(script),
                            str(path),
                            "--backend",
                            "production",
                            "--library-root",
                            root,
                            "--threads",
                            str(args.threads),
                            "--result-type",
                            args.result_type,
                        ]
                        for name in ("features", "usecols"):
                            value = getattr(args, name)
                            if value:
                                command.extend(["--" + name, *value])
                        run = subprocess.run(
                            command, env=environment, capture_output=True, text=True, check=True
                        )
                        record = json.loads(run.stdout)
                        key = path.name
                        content = record["checksum"], record["rows"], record["columns"]
                        if key in expected and expected[key] != content:
                            raise RuntimeError(
                                f"Reader results differ: {path}, {label}, {interpreter}"
                            )
                        expected[key] = content
                        record.update(path=path.name, release=label, repeat=repeat)
                        records.append(record)
                        raw.write(json.dumps(record, sort_keys=True) + "\n")
                        raw.flush()
                        print(
                            record["versions"]["pandas"],
                            repeat,
                            path.name,
                            label,
                            round(record["wall_seconds"], 3),
                            round(record["peak_rss_mib"]),
                            flush=True,
                        )
    result = {
        "description": "Complete production readers; imports/checksums outside wall time; sequential fresh processes",
        "date": datetime.now(timezone.utc).isoformat(),
        "threads": args.threads,
        "inputs": {
            path.name: {
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in args.inputs
        },
        "summaries": summarize(records),
        "measurements": records,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
