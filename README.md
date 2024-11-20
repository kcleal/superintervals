SuperIntervals
==============

A fast, memory-efficient data structure for interval intersection queries.
SuperIntervals uses a novel superset-index approach that maintains 
intervals in position-sorted order, enabling cache-friendly searches and SIMD-optimized counting.

### Features:

- Linear-time index construction from sorted intervals
- Cache-friendly querying
- SIMD acceleration (AVX2/Neon) for counting operations
- Minimal memory overhead (one size_t per interval)
- Available for C++, Rust, Python, and C
- Optional Eytzinger memory layout for slightly faster queries (C++/Rust only)
- No dependencies, header only


## Quick Start

- Intervals are considered end-inclusive 
- The index() function must be called before any queries
- Found intervals are returned in reverse position-sorted order

### 🐍 Python

```python
from superintervals import IntervalSet

iset = IntervalSet()
iset.add(10, 20, 'A')
iset.index()
overlaps = iset.find_overlaps(8, 20)
```

### ⚙️ C++
```cpp
#include "SuperIntervals.hpp"

SuperIntervals<int, std::string> intervals;
intervals.add(1, 5, "A");
intervals.index();
std::vector<std::string> results;
intervals.findOverlaps(4, 9, results);
```

### 🦀 Rust

```rust
use super_intervals::SuperIntervals;

let mut intervals = SuperIntervals::new();
intervals.add(1, 5, "A");
intervals.index();
let mut results = Vec::new();
intervals.find_overlaps(4, 11, &mut results);
```


## Test programs
Test programs expect plain text BED files and only assess chr1 records - other chromosomes are ignored.

C++ program compares SuperIntervals, ImplicitIntervalTree, IntervalTree and NCLS:
```
cd test; make
./run-cpp-libs a.bed b.bed
```

Rust program:
```
RUSTFLAGS="-Ctarget-cpu=native" cargo run --release --example bed-intersect-si
cargo run --release --example bed-intersect-si a.bed b.bed
```

## Benchmark

Update when issue #2 is resolved.

## Python

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

## Cpp

```cpp
#include <iostream>
#include <vector>
#include "SuperIntervals.hpp"

int main() {
    // Create a SuperIntervals instance for integer intervals with string data
    // Specify with S, T template types
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
There is also a `SuperIntervalsEytz` subclasses that can be used. `SuperIntervalsEytz` 
uses an Eytzinger memory layout that can sometimes offer faster query times at the cost of higher memory
usage and slower indexing time.

## Rust

Fetch using cargo add.

```
use super_intervals::SuperIntervals;

fn main() {
    // Create a new instance of SuperIntervals
    let mut intervals = SuperIntervals::new();

    // Add some intervals with associated data of type T
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

## Acknowledgements

- The rust test program borrows heavily from the coitrees package
- The superset-index implemented here exploits a similar interval ordering as described in
Schmidt 2009 "Interval Stabbing Problems in Small Integer Ranges". However, the superset-index has several advantages including
  1. An implicit memory layout
  1. General purpose implementation (not just small integer ranges)
  1. SIMD counting algorithm 
