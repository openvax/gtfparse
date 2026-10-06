"""Backend experiment, not an alternative supported reader.

Run each invocation in a fresh process. See README.md for corpus provenance,
repetition, equivalence checks, and limitations of the experimental readers.
"""

import argparse
import csv
import gzip
import hashlib
import importlib
import json
import os
import platform
import resource
import sys
from pathlib import Path
from time import perf_counter, process_time

import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv

try:
    import polars as pl
except ModuleNotFoundError as error:
    if error.name != "polars":
        raise
    pl = None

# Use the requested checkout before importing gtfparse, for baseline comparisons.
import_options = argparse.ArgumentParser(add_help=False)
import_options.add_argument("--library-root", type=Path)
import_args, _ = import_options.parse_known_args()
library_root = import_args.library_root or Path(__file__).resolve().parents[1]
if not (library_root / "gtfparse" / "__init__.py").is_file():
    import_options.error("--library-root must contain the gtfparse package")
sys.path.insert(0, str(library_root))
package = importlib.import_module("gtfparse")
__version__ = package.__version__
expand_attribute_strings = package.expand_attribute_strings
read_gtf = package.read_gtf

reader = importlib.import_module("gtfparse.read_gtf")
CATEGORIES = ("seqname", "source", "feature", "strand")


class Timings:
    def __init__(self):
        self.stages = {}
        self.starts = {}

    def callback(self, stage, completed, total):
        if completed == 0:
            self.starts[stage] = perf_counter()
        if total is not None and completed == total:
            self.stages[stage] = perf_counter() - self.starts[stage]

    def run(self, stage, fn):
        start = perf_counter()
        result = fn()
        self.stages[stage] = perf_counter() - start
        return result


def leading_comments(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        count = 0
        for line in stream:
            if not line.startswith("#"):
                break
            count += 1
    return count


def pandas_read(path, features):
    # Explicit object strings match Polars.to_pandas and avoid pandas' default
    # NA vocabulary interpreting valid names such as "NA" as missing values.
    types = dict.fromkeys(CATEGORIES, "category")
    types.update(start="int64", end="int64", score="float32", frame="float64", attribute=object)
    frame = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=reader.REQUIRED_COLUMNS,
        dtype=types,
        engine="c",
        quoting=csv.QUOTE_NONE,
        skiprows=leading_comments(path),
        keep_default_na=False,
        na_values=[".", ""],
    )
    frame = frame.loc[frame.seqname.notna()]
    frame["frame"] = frame["frame"].fillna(0).astype("uint32")
    if features:
        frame = frame.loc[frame.feature.isin(features)]
    frame["attribute"] = frame.attribute.where(frame.attribute.notna(), None)
    return frame.reset_index(drop=True)


def arrow_read(path, features, extension=False):
    types = dict.fromkeys(CATEGORIES, pa.dictionary(pa.int32(), pa.string()))
    types.update(
        start=pa.int64(),
        end=pa.int64(),
        score=pa.float32(),
        frame=pa.uint32(),
        attribute=pa.string(),
    )
    table = pacsv.read_csv(
        path,
        read_options=pacsv.ReadOptions(column_names=reader.REQUIRED_COLUMNS),
        parse_options=pacsv.ParseOptions(
            delimiter="\t",
            quote_char=False,
            invalid_row_handler=lambda row: "skip" if row.text.startswith("#") else "error",
        ),
        convert_options=pacsv.ConvertOptions(
            column_types=types, null_values=[".", ""], strings_can_be_null=True
        ),
    )
    table = table.filter(pc.is_valid(table["seqname"]))
    table = table.set_column(7, "frame", pc.fill_null(table["frame"], 0))
    if features:
        table = table.filter(pc.is_in(table["feature"], value_set=pa.array(features)))
    return table.to_pandas(types_mapper=pd.ArrowDtype if extension else None)


def prototype(path, backend, timings, features, usecols, result_type):
    if backend == "pandas":
        frame = timings.run("read", lambda: pandas_read(path, features))
    else:
        frame = timings.run(
            "read", lambda: arrow_read(path, features, backend == "arrow-extension")
        )
    attributes = frame.pop("attribute")
    # Arrow extension scalars use pd.NA; the production expander accepts None.
    values = (None if value is pd.NA else value for value in attributes)
    if backend == "arrow-extension":
        # Materialization belongs in the measured expansion time.
        class SizedValues:
            def __len__(self):
                return len(attributes)

            def __iter__(self):
                return values

        attributes_for_expansion = SizedValues()
    else:
        attributes_for_expansion = attributes
    expanded = timings.run(
        "attributes", lambda: expand_attribute_strings(attributes_for_expansion, usecols=usecols)
    )
    # Match pandas 2 object strings / pandas 3 inferred StringDtype, as returned
    # by the production reader's Arrow-to-pandas conversion on each version.
    frame = timings.run("assemble", lambda: pd.concat([frame, pd.DataFrame(expanded)], axis=1))
    attributes = attributes_for_expansion = values = expanded = None

    def finish():
        nonlocal frame
        if backend == "arrow-extension":
            for name in CATEGORIES:
                # pandas 2.3 cannot directly cast an Arrow dictionary with
                # null indices to category (its zero-copy conversion fails).
                # Include the object bridge in this prototype's measured cost.
                frame[name] = frame[name].astype(object).astype("category")
            for name, dtype in {
                "start": "int64",
                "end": "int64",
                "score": "float32",
                "frame": "uint32",
            }.items():
                frame[name] = frame[name].astype(dtype)
        frame = reader._cast_version_columns(frame)
        if usecols is not None:
            frame = frame[[name for name in usecols if name in frame]]
        return pl.from_pandas(frame) if result_type == "polars" else frame

    return timings.run("convert", finish)


def checksum(frame):
    """Hash all values and column order, ignoring categorical dictionary order."""
    if pl is not None and isinstance(frame, pl.DataFrame):
        frame = frame.to_pandas()
    digest = hashlib.sha256(json.dumps(list(frame.columns)).encode())
    for name in frame.columns:
        values = frame[name]
        if pd.api.types.is_numeric_dtype(values.dtype):
            values = values.astype("float64")
        else:
            values = values.astype(object).where(values.notna(), None)
        digest.update(pd.util.hash_pandas_object(values, index=False).values.tobytes())
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path)
    parser.add_argument("--library-root", type=Path, help="Read from another gtfparse checkout")
    parser.add_argument(
        "--backend",
        choices=["production", "pandas", "arrow", "arrow-extension"],
        default="production",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--result-type", choices=["pandas", "polars"], default="pandas")
    parser.add_argument("--features", nargs="+")
    parser.add_argument("--usecols", nargs="+")
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Compare to production output after timing (small inputs).",
    )
    args = parser.parse_args()
    if args.result_type == "polars" and pl is None:
        parser.error("Polars output requires pip install 'gtfparse[polars]'")
    if pl is not None and pl.thread_pool_size() != args.threads:
        parser.error("Set POLARS_MAX_THREADS before starting the process to match --threads")
    pa.set_cpu_count(args.threads)
    timings = Timings()
    start, cpu_start = perf_counter(), process_time()
    if args.backend == "production":
        frame = read_gtf(
            str(args.path),
            result_type=args.result_type,
            features=args.features,
            usecols=args.usecols,
            progress_callback=timings.callback,
        )
    else:
        frame = prototype(
            str(args.path), args.backend, timings, args.features, args.usecols, args.result_type
        )
    elapsed, cpu = perf_counter() - start, process_time() - cpu_start
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_mib = peak / (1024**2 if sys.platform == "darwin" else 1024)
    # Capture the high-water mark BEFORE hashing, comparison, or output sizing.
    timings.stages["other"] = elapsed - sum(timings.stages.values())
    output = {
        "path": str(args.path),
        "input_bytes": args.path.stat().st_size,
        "backend": args.backend,
        "result_type": args.result_type,
        "features": args.features,
        "usecols": args.usecols,
        "threads": args.threads,
        "wall_seconds": elapsed,
        "cpu_seconds": cpu,
        "peak_rss_mib": peak_mib,
        "stages": timings.stages,
        "rows": len(frame),
        "columns": list(frame.columns),
        "checksum": checksum(frame),
        "versions": {
            "gtfparse": __version__,
            "pandas": pd.__version__,
            "polars": pl.__version__ if pl is not None else None,
            "pyarrow": pa.__version__,
            "python": platform.python_version(),
        },
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
    }
    if args.verify:
        reference = read_gtf(
            str(args.path), result_type="pandas", features=args.features, usecols=args.usecols
        )
        actual = frame.to_pandas() if pl is not None and isinstance(frame, pl.DataFrame) else frame
        pd.testing.assert_frame_equal(reference, actual, check_categorical=False)
        assert checksum(reference) == output["checksum"]
        output["verified"] = True
    print(json.dumps(output, sort_keys=True))


if __name__ == "__main__":
    main()
