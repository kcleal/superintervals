Superintervals
==============

This a python or c++ header-only library for finding
interval intersections, with genomics data in mind.

A novel superset-index is utilised to speed up intersection
queries. Counting of interval intersections uses AVX2/Neon 
intrinsics if available.

Compile test programs
---------------------

Cpp programs

`cd test; make`

Rust program:

`RUSTFLAGS="-Ctarget-cpu=native" cargo run --release --example bed-intersect-si`

Python usage
------------
Install using `pip install superintervals`

```
from superintervals import IntervalSet

iset = IntervalSet()

# Add interval start, end, identifier. Integer values are supported
iset.add(10, 20, 0)
iset.add(19, 18, 1)
iset.add(8, 11, 2)

# Index method must be called before queries
iset.index()

iset.any_overlaps(8, 20)
# >>> True

iset.count_overlaps(8, 20)
# >>> 3

iset.find_overlaps(8, 20)
# >>> [1, 0, 2]

iset.set_search_interval(8, 20)
for itv in iset:
    print(itv)

# >>> (19, 18, 1) 
# >>> (10, 20, 0) 
# >>> (8, 11, 2)

```