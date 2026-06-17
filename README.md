SuperIntervals
==============

A fast, memory-efficient data structure for interval intersection queries.
SuperIntervals uses a novel superset-index approach that maintains 
intervals in position-sorted order, enabling cache-friendly searches and SIMD-optimized counting.

Available for [C++](#cpp), [Rust](#rust), [Python](#python).
The R package is hosted at https://github.com/kcleal/superintervalsr.

### Features:

- Linear-time index construction from sorted intervals
- Cache-friendly querying
- SIMD acceleration (AVX2/Neon) for counting operations
- Small memory overhead (one size_t per interval)
- Optional Eytzinger memory layout for slightly faster queries (C++ only)
- No dependencies, header only

### Notes:

- Intervals are considered end-inclusive 
- The build() function must be called before any queries
- Found intervals are returned in **reverse** position-sorted order
- Set operations (merge_overlaps, union, intersection, difference, symmetric_difference, gaps, span) are available in all three languages — see the per-language API docs

## Python

Install using `pip install superintervals`

```python
from superintervals import IntervalMap

imap = IntervalMap()
imap.add(10, 20, 'A')
imap.build()
results = imap.search_values(8, 20)  # ['A']
```

Python API documentation can be found here:
https://github.com/kcleal/superintervals/blob/main/src/superintervals/README.md


## Cpp

Header only implementation, copy to your include directory.

```cpp
#include "SuperIntervals.hpp"

si::IntervalMap<int, std::string> imap;
imap.add(1, 5, "A");
imap.build();
std::vector<std::string> results;
imap.search_values(4, 9, results);
```

C++ API documentation can be found here:
https://github.com/kcleal/superintervals/blob/main/src/README_cpp.md


## Rust

Add to your project using `cargo add superintervals`

```rust
use superintervals::IntervalMap;

let mut imap = IntervalMap::new();
imap.add(1, 5, "A");
imap.build();
let mut results = Vec::new();
imap.search_values(4, 11, &mut results);
```

Rust API documentation can be found here:
https://github.com/kcleal/superintervals/blob/main/src/README_rust.md



## Test programs
Test programs expect plain text BED files and only assess chr1 records - other chromosomes are ignored.

C++ program compares SuperIntervals, ImplicitIntervalTree, IntervalTree and NCLS:
```
cd test; make
./run-tests
./run-cpp-libs a.bed b.bed
```

Rust program:
```
RUSTFLAGS="-Ctarget-cpu=native" cargo run --release --example bed-intersect-si
cargo run --release --example bed-intersect-si a.bed b.bed
```

Python program:
```
python test/run-py-libs.py a.bed b.bed
```

R program:
```
Rscript src/R/benchmark.R
```

## Algorithm

SuperIntervals keeps intervals in a flat array sorted by start position. Alongside it
a single auxiliary "branch" array (one `size_t` per interval) records, for each
interval, the index to jump back to in order to skip over a contiguous block of
intervals that cannot overlap the query — effectively an implicit, branch-free index
built in a single linear pass over the sorted data.

A query first binary-searches for the last interval whose start is `<= query_end`,
then walks backwards. At each step the branch array lets it leap past whole runs of
non-overlapping intervals instead of testing them one by one. Because the data is
stored contiguously and accessed in descending order, the walk is highly
cache-friendly, and the overlap count can be computed with SIMD (AVX2/Neon) over the
packed end coordinates. This is why results come back in **reverse** position-sorted
order — it falls naturally out of the backward walk.

The result is a structure with very low memory overhead (one extra integer per
interval, no pointers or tree nodes) that is fast to build and fast to query.

## Benchmark

SuperIntervals was compared against Coitrees, Implicit Interval Tree, Interval Tree
and NCLS across a range of genomic datasets. Finding interval intersections was on
average ~1.5–3x faster than the next-best library (Coitrees for Rust, Implicit
Interval Tree for C++), and SIMD counting performance was on par with Coitrees.

Full methodology, datasets and per-dataset timing tables are in
[test/benchmark.md](test/benchmark.md).

### Third-party confirmation

An independent September 2025 benchmark by the [polars-bio](https://biodatageeks.org)
team ([report](https://biodatageeks.org/polars-bio/blog/2025/09/05/interval-operations-benchmark--update-september-2025))
evaluated SuperIntervals as a backend data structure for polars-bio. They found it to
be **consistently the fastest or tied for fastest**, delivering 1.25–1.44x speedups
over the default COITrees implementation across small, medium and large datasets,
and concluded it offers reliable, well-rounded performance without worst-case
degradation — making it an ideal default for general-purpose interval operations.


## Acknowledgements

- The rust test program borrows heavily from the coitrees package
- The superset-index implemented here exploits a similar interval ordering as described in
Schmidt 2009 "Interval Stabbing Problems in Small Integer Ranges". However, the superset-index has several advantages including
  1. An implicit memory layout
  1. General purpose implementation (not just small integer ranges)
  1. SIMD counting algorithm 
- The Eytzinger layout was adapted from Sergey Slotin, Algorithmica
