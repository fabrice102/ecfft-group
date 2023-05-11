#!/usr/bin/env python3

"""
Export the benches results from `cargo bench` into CSV and Markdown formats.
Assume it is run in the folder containing `target/criterion/...`

To run unit tests, run
$ python3 -m doctest benches_table.py --verbose
"""

import argparse
from datetime import datetime, timedelta
import json
from pathlib import Path
from typing import IO

DEFAULT_CRITERION_TARGET_FOLDER = "target/criterion"
DEFAULT_CRITERION_GLOB = "**/new/estimates.json"
DEFAULT_BENCH_NAMES = ",".join([
    "ed25519-poly/pt-smallDm-hornerSmall",
    "ed25519-poly/pt-ecfftDm-extend"
])


def read_estimate(full_path: Path) -> float:
    """
    Read the estimate.json at full_path and return the mean point estimate from it
    Throw an exception if the estimate does not exist or is not a floating number
    """
    # read the estimates.json file
    with full_path.open() as f:
        j = json.load(f)

    # check the estimate exists and is a float/integer
    if "mean" not in j or "point_estimate" not in j["mean"] or \
            (not isinstance(j["mean"]["point_estimate"], float) and not isinstance(j["mean"]["point_estimate"],
                                                                                   int)):
        raise ValueError(f"bench estimate '{full_path.as_posix()}' does "
                         f"not contain a float 'mean.point_estimate'")
    estimate = float(j["mean"]["point_estimate"])

    return estimate


def parse_benches_results(
        criterion_target_folder: str,
        criterion_glob: str,
        bench_names: list[str]
) -> tuple[dict[int, dict[str, float]], datetime, datetime]:
    """
    Parse benches results from criterion_target_folder that satisfy criterion_glob pattern
    See default in global constant

    :param criterion_target_folder: folder where the criterion benchmark results are
    :param criterion_glob: glob pattern to find the criterion benchmarks (must match files estimates.json)
    :param bench_names: list of bench names that are considered
    :return: (results, earliest_bench, latest_earliest) where 
        results[bench_param][bench_name] = point estimate number of ns
        earliest/latest_bench are the datetime of the earliest/latest bench made
    """
    results = {}  # results[bench_param][bench_name] = point estimate number of ns
    latest_bench = 0  # timestamp of the latest bench that is seen
    earliest_bench = 3e10  # timestamp of the earliest bench that is seen (large number that will never be reached)

    criterion_results = Path(criterion_target_folder)
    for full_path in criterion_results.glob(criterion_glob):
        bench_name = full_path.relative_to(criterion_results).parents[2].as_posix()

        # ignore benchmarks that are not in bench_names
        if bench_name not in bench_names:
            continue

        bench_param_str = full_path.parents[1].name
        if not bench_param_str.isdigit():
            raise ValueError("'{}' last part is not an integer".format(full_path.as_posix()))

        bench_param = int(bench_param_str)

        # update earliest/latest bench
        bench_timestamp = full_path.stat().st_mtime
        earliest_bench = min(bench_timestamp, earliest_bench)
        latest_bench = max(bench_timestamp, latest_bench)

        estimate = read_estimate(full_path)

        # store it
        if bench_param not in results:
            results[bench_param] = {}
        if bench_name in results[bench_param]:
            # the same benchmark name/param can only appear once
            # otherwise there is an issue somewhere
            raise ValueError(f"bench '{bench_name}' with param '{bench_param}' seen twice")

        results[bench_param][bench_name] = estimate

    return results, datetime.fromtimestamp(earliest_bench), datetime.fromtimestamp(latest_bench)


def format_ns(ns: float) -> str:
    """
    Nicely format a ns time into a string
    :param ns:
    :return:

    >>> format_ns(10)
    '10.0 ns'
    >>> format_ns(10.5)
    '10.5 ns'
    >>> format_ns(100002)
    '100 us'
    >>> format_ns(2.126e6)
    '2.13 ms'
    >>> format_ns(1.234e10)
    '12.3 s'
    """
    if ns < 1e3:
        return f"{ns:#.3g} ns".replace(". ", "  ")
    if ns < 1e6:
        us = ns / 1e3
        return f"{us:#.3g} us".replace(". ", "  ")
    if ns < 1e9:
        ms = ns / 1e6
        return f"{ms:#.3g} ms".replace(". ", "  ")
    s = ns / 1e9
    if s > 1000:
        return f"{round(s):,} s"
    else:
        return f"{s:#.3g} s".replace(". ", "  ")


def format_md(bench_names: list[str], short_bench_names: list[str], results: dict[int, dict[str, float]]):
    out = ""

    # rows = params = criterion parameters
    params = sorted(results)

    # columns = benchmark names
    # compute column size (max length of each column)
    # note: min of length 3 for Github table
    column_sizes = {
        bench_name: max(
            max(
                (
                    len(format_ns(results[param][bench_name]))
                    for param in results
                    if bench_name in results[param]
                ),
                default=0
            ),
            len(short_bench_name),
            3
        ) for short_bench_name, bench_name in zip(short_bench_names, bench_names)
    }
    # parameter line length
    param_column_size = max(3, max((len(f"{param}") for param in params), default=0))

    # header line with the short bench names
    out += "| " + "".ljust(param_column_size) + " | "
    out += " | ".join(
        short_bench_name.ljust(column_sizes[bench_name])
        for short_bench_name, bench_name in zip(short_bench_names, bench_names)
    )
    out += " |\n"

    # separator line
    out += "|-" + "-" * param_column_size + "-|-"
    out += "-|-".join("-" * column_sizes[bench_name] for bench_name in bench_names)
    out += "-|\n"

    for param in params:
        out += "| " + f"{param}".ljust(param_column_size) + " | "
        columns = (
            format_ns(results[param][bench_name]) if bench_name in results[param] else ""
            for short_bench_name, bench_name in zip(short_bench_names, bench_names)
        )
        out += " | ".join(
            col.ljust(column_sizes[bench_name])
            for bench_name, col in zip(bench_names, columns)
        )
        out += " |\n"

    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=["md"], default="md")
    parser.add_argument("--criterion-target-folder", default=DEFAULT_CRITERION_TARGET_FOLDER)
    parser.add_argument("--criterion-glob", default=DEFAULT_CRITERION_GLOB)
    parser.add_argument("--bench-names",
                        help="comma-separated list of bench names in the order they will appear",
                        default=DEFAULT_BENCH_NAMES)
    parser.add_argument("--short-bench-names",
                        help="comma-separated short names for benches used for columns,"
                             "if empty, use same as --bench-names",
                        default="")
    args = parser.parse_args()

    bench_names = [name.strip() for name in args.bench_names.split(",")]
    if args.short_bench_names == "":
        short_bench_names = bench_names
    else:
        short_bench_names = [name.strip() for name in args.short_bench_names.split(",")]
        if len(short_bench_names) != len(bench_names):
            raise ValueError("--short-bench-names and --bench-names need to have same length")
    results, earliest_bench, latest_bench = parse_benches_results(
        args.criterion_target_folder,
        args.criterion_glob,
        bench_names
    )

    print(f"earliest bench: {earliest_bench.isoformat()}")
    print(f"latest bench:   {latest_bench.isoformat()}")
    print()
    if latest_bench - earliest_bench > timedelta(days=1):
        print("WARNING: The latest bench was made more than 1 day after the earliest bench")
        print("         Are you sure they have been made with the same software version?")
        print()

    if args.format == "md":
        print(format_md(bench_names, short_bench_names, results))
    else:
        raise NotImplementedError()


if __name__ == "__main__":
    main()
