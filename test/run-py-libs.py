import time
import quicksect
from quicksect import Interval
from ncls import NCLS
import cgranges as cr

from superintervals import IntervalMap

import numpy as np
import pandas as pd
import sys


def load_intervals(intervals_path, queries_path):
    queries = []
    intervals = []
    with open(intervals_path, "r") as f:
        for line in f:
            l = line.split("\t")
            intervals.append((int(l[1]), int(l[2])))
    with open(queries_path, "r") as f:
        for line in f:
            l = line.split("\t")
            queries.append((int(l[1]), int(l[2])))
    return np.array(intervals), np.array(queries)


def to_micro(t0):
    return int((time.time() - t0) * 1000000)


def run_tools(intervals, queries):
    # superintervals - Updated for new API
    t0 = time.time()
    sitv = IntervalMap()
    for i, (s, e) in enumerate(intervals):
        sitv.add(s, e, i)  # Using index as data
    sitv.build()
    build_time = to_micro(t0)
    print(f"SuperIntervals-py,{build_time},", end='')

    # Test search_values (returns data objects)
    t0 = time.time()
    v = 0
    for start, end in queries:
        a = sitv.search_values(start, end)
        v += len(a)
    search_time = to_micro(t0)
    print(f'{search_time},{v},', end='')

    # Test count method
    t0 = time.time()
    v = 0
    for start, end in queries:
        v += sitv.count(start, end)
    count_time = to_micro(t0)
    print(f'{count_time},{v}')

    # quicksect
    t0 = time.time()
    tree = quicksect.IntervalTree()
    for s, e in intervals:
        tree.add(s, e)
    quicksect_build = to_micro(t0)
    print(f"Quicksect,{quicksect_build},", end='')

    t0 = time.time()
    v = 0
    for start, end in queries:
        a = tree.find(Interval(start, end))
        v += len(a)
    quicksect_search = to_micro(t0)
    print(f'{quicksect_search},{v}')

    # cgranges
    t0 = time.time()
    cg = cr.cgranges()
    for s, e in intervals:
        cg.add("1", s, e + 1, 0)
    cg.index()
    cgranges_build = to_micro(t0)
    print(f"Cgranges-py,{cgranges_build},", end='')

    t0 = time.time()
    v = 0
    for start, end in queries:
        a = list(cg.overlap("1", start, end + 1))
        v += len(a)
    cgranges_search = to_micro(t0)
    print(f'{cgranges_search},{v}')

    # ncls
    t0 = time.time()
    starts = pd.Series(intervals[:, 0])
    ends = pd.Series(intervals[:, 1])
    treencls = NCLS(starts, ends, starts)
    ncls_build = to_micro(t0)
    print(f"NCLS-py,{ncls_build},", end='')

    t0 = time.time()
    v = 0
    for start, end in queries:
        a = list(treencls.find_overlap(start - 1, end + 1))
        v += len(a)
    ncls_search = to_micro(t0)
    print(f'{ncls_search},{v}')


def run_superintervals_detailed_benchmark(intervals, queries):
    """
    Run a more detailed benchmark of superintervals methods
    """
    print("\n=== SuperIntervals Detailed Benchmark ===")

    # Setup
    sitv = IntervalMap()
    for i, (s, e) in enumerate(intervals):
        sitv.add(s, e, f"data_{i}")

    t0 = time.time()
    sitv.build()
    print(f"Index time: {to_micro(t0)} microseconds")

    # Test different search methods
    methods = [
        ("search_values", lambda start, end: len(sitv.search_values(start, end))),
        ("search_items", lambda start, end: len(sitv.search_items(start, end))),
        ("search_keys", lambda start, end: len(sitv.search_keys(start, end))),
        ("search_idxs", lambda start, end: len(sitv.search_idxs(start, end))),
        ("count", lambda start, end: sitv.count(start, end)),
        ("has_overlaps", lambda start, end: int(sitv.has_overlaps(start, end))),
        ("iter_values", lambda start, end: sum(1 for _ in sitv.iter_values(start, end))),
        ("iter_items", lambda start, end: sum(1 for _ in sitv.iter_items(start, end))),
        ("iter_idxs", lambda start, end: sum(1 for _ in sitv.iter_idxs(start, end))),
        ("iter_keys", lambda start, end: sum(1 for _ in sitv.iter_keys(start, end))),
    ]

    # Verify streaming iterators return the same hits, in the same descending
    # order, as the list-based searches (reverse position-sorted, per the
    # algorithm; see the main README).
    check_queries = [(int(s), int(e)) for s, e in queries[: min(200, len(queries))]]
    lo = int(intervals[:, 0].min())
    hi = int(intervals[:, 1].max())
    mid = (lo + hi) // 2
    check_queries += [(hi + 1000, hi + 2000), (lo, hi), (mid, mid)]
    for start, end in check_queries:
        assert list(sitv.iter_idxs(start, end)) == sitv.search_idxs(start, end), (start, end)
        assert list(sitv.iter_keys(start, end)) == sitv.search_keys(start, end), (start, end)
        assert list(sitv.iter_values(start, end)) == sitv.search_values(start, end), (start, end)
        assert list(sitv.iter_items(start, end)) == sitv.search_items(start, end), (start, end)
    print("Iterator correctness checks passed")

    # has_overlaps must agree with the full search (regression: it used to check
    # only the single candidate at upper_bound(end), missing a long interval
    # hidden behind a shorter nested one).
    for start, end in check_queries:
        assert sitv.has_overlaps(start, end) == bool(sitv.search_idxs(start, end)), (start, end)
    nested_case = IntervalMap()
    nested_case.add(0, 99, "long")
    nested_case.add(10, 19, "nested")
    nested_case.build()
    assert nested_case.has_overlaps(30, 39)
    assert not nested_case.has_overlaps(100, 200)
    print("has_overlaps consistency checks passed")

    # Mutation during iteration invalidates outstanding iterators
    probe = IntervalMap()
    probe.add(0, 10, "a")
    probe.add(5, 15, "b")
    probe.build()
    pit = probe.iter_items(0, 20)
    assert sorted(pit) == [(0, 10, "a"), (5, 15, "b")]
    try:
        next(pit)
        raise AssertionError("expected StopIteration from exhausted iterator")
    except StopIteration:
        pass
    pit = probe.iter_items(0, 20)
    next(pit)
    probe.add(20, 30, "c")
    try:
        next(pit)
        raise AssertionError("expected RuntimeError after modifying map during iteration")
    except RuntimeError:
        pass
    print("Iterator mutation-guard checks passed")

    for method_name, method_func in methods:
        t0 = time.time()
        total_results = 0
        for start, end in queries:
            total_results += method_func(start, end)
        elapsed = to_micro(t0)
        print(f"{method_name}: {elapsed} microseconds, {total_results} total results")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python benchmark.py <intervals_file> <queries_file>")
        print("Both files should be in BED format (tab-separated: chr start end ...)")
        sys.exit(1)

    intervals_file = sys.argv[1]
    queries_file = sys.argv[2]

    print("Tool,BuildTime(μs),SearchTime(μs),ResultCount")

    try:
        intervals, queries = load_intervals(intervals_file, queries_file)
        print(f"# Testing with {len(intervals)} intervals and {len(queries)} queries", file=sys.stderr)

        run_tools(intervals, queries)

        # Run detailed superintervals benchmark
        run_superintervals_detailed_benchmark(intervals, queries)

    except FileNotFoundError as e:
        print(f"Error: Could not find file - {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)