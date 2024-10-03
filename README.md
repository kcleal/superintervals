Superintervals
==============

This a static data structure for find interval intersection, fast! 
There are implementations for Python, C, C++ (header-only) and Rust.

A novel superset-index is utilised to speed up intersection
queries. Counting of interval intersections uses AVX2/Neon 
intrinsics if available.

Benchmark
---------

C++ program compares SuperIntervals, ImplicitIntervalTree, IntervalTree, NCLS and FastStabbing algorithm:
```
cd test; make
./run-cpp-libs a.bed b.bed
```

Rust program:
```
RUSTFLAGS="-Ctarget-cpu=native" cargo run --release --example bed-intersect-si
./target/release/examples/bed-intersect-si a.bed b.bed
```

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

C++ usage
---------
```cpp
#include <iostream>
#include <vector>
#include "SuperIntervals.hpp"

int main() {
    // Create a SuperIntervals instance for integer intervals with string data
    SuperIntervals<int, std::string> intervals;

    // Add some intervals
    intervals.add(1, 5, "Interval A");
    intervals.add(3, 7, "Interval B");
    intervals.add(6, 10, "Interval C");
    intervals.add(8, 12, "Interval D");

    // Index the intervals (must be called before querying)
    intervals.index();

    // Find overlaps for the range [4, 9]
    std::vector<std::string> overlaps;
    intervals.findOverlaps(4, 9, overlaps);

    // Print the overlapping intervals
    for (const auto& interval : overlaps) {
        std::cout << interval << std::endl;
    }
    
    // Count the intervals instead
    std::cout << "Count: " << intervals.countOverlaps(4, 9) << std::endl;
    
    // Count stabbed intervals at point 7
    std::cout << "Number of intervals containing point 7: " << intervals.countStabbed(7) << std::endl;

    return 0;
}
```
There are also `SuperIntervalsEytz` and `SuperIntervalsDense` subclasses that can be used. `SuperIntervalsEytz` 
uses an Eytzinger memory layout that can sometimes offer faster query times at the cost of higher memory
usage and slower indexing time. The `SuperIntervalsDense` class does not use binary search, and instead uses a lookup
table. For this reason it is much slower to index and uses a lot of memory, but is usually faster for queries.

Rust usage
----------

```
use super_intervals::SuperIntervals;

fn main() {
    // Create a new instance of SuperIntervals
    let mut intervals = SuperIntervals::new();

    // Add some intervals with associated data
    intervals.add(1, 5, "Interval A");
    intervals.add(10, 15, "Interval B");
    intervals.add(7, 12, "Interval C");

    // Call index() to prepare the intervals for queries
    intervals.index();

    // Query for overlapping intervals with a range (4, 11)
    let mut found_intervals = Vec::new();
    intervals.find_overlaps(4, 11, &mut found_intervals);
    
    // Display found intervals
    for interval in found_intervals {
        println!("Found overlapping interval: {}", interval);
    }

    // Count overlaps with a range (4, 11)
    let overlap_count = intervals.count_overlaps(4, 11);
    println!("Number of overlapping intervals: {}", overlap_count);
}
```
There is also `SuperIntervalsEytz` implementation. `SuperIntervalsEytz` 
uses an Eytzinger memory layout that can sometimes offer faster query times at the cost of higher memory
usage and slower indexing time.