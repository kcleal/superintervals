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

SuperIntervals (SI) was compared with:
- Coitrees (Rust: https://github.com/dcjones/coitrees)
- Implicit Interval Tree (C++: https://github.com/lh3/cgranges)
- Interval Tree (C++: https://github.com/ekg/intervaltree)
- Nested Containment List (C: https://github.com/pyranges/ncls/tree/master/ncls/src)

Main results:
- Roughly ~2-3x faster than the next best library (Coitrees for Rust, Implicit Interval Tree for C++)

### Datasets:
1. Random regions generated using bedtools
2. RNA-seq reads and annotations from cgranges repository
3. ONT reads from sample PAO33946 (chr1, chrM)
4. Paired-end reads from sample DB53, NCBI BioProject PRJNA417592, (chr1, chrM)
5. UCSC genes from hg19

Test programs use internal timers and print data to stdout, measuring the 
index time, and time to find all intersections. Other steps such as file IO are ignored. Test programs also 
only assess chr1 bed records - other chromosomes are ignored. For 'chrM' records,
the M was replaced with 1 using sed. Data were assessed in position sorted and random order.
Datasets can be found on the Releases page, and the `test/run_tools.sh` script has instructions
for how to repeat the benchmark.

Timings were in microseconds using an i9-11900K, 64 GB, 2TB NVMe machine.


### 1. Finding interval intersections

- Coitrees-s uses the `SortedQuerent` version of coitrees
- SI = superintervals. Eytz refers to the eytzinger layout. `-rs` is the Rust implementation.

#### Intervals in sorted order

|                       | Coitrees | Coitrees-s | SI-rs       | SI-rs     | ImplicitITree-C++ | IntervalTree-C++ | NCLS-C   | SI-C++  | SI-Eytz-C++ |
| --------------------- | -------- | ---------- |-------------|-----------| ----------------- | ---------------- | -------- |---------|-------------|
| DB53 reads, ONT reads | 1649.6   | 3169       | 732         | **729**   | 3802.6            | 46393.8          | 10833.6  | 1391.6  | **1365.6**  |
| DB53 reads, genes     | 54.2     | 82.8       | **21**      | **21**    | 121.6             | 108              | 292.8    | 43      | **40.2**    |
| ONT reads, DB53 reads | 6487.2   | 3437.2     | 534.6       | **533.6** | 18067.4           | 12448            | 31466.2  | 5333.2  | **4545.2**  |
| anno, rna             | 49.6     | 33.6       | 17.2        | **17**    | 127.2             | 91.2             | 210.6    | 31.2    | **21.2**    |
| genes, DB53 reads     | 1171     | 992.8      | 270         | **269.2** | 3141              | 1339.8           | 1768     | 441.8   | **315**     |
| mito-b, mito-a        | 35046.2  | 35134      | **13115.2** | 13117.2   | 95137.4           | 108567.8         | 250671.8 | 33703.8 | **33298.6** |
| rna, anno             | 31.8     | 22.6       | **4**       | **4**     | 71.2              | 54               | 238.8    | 29.4    | **27.2**    |

#### Intervals in random order

|                       | Coitrees | Coitrees-s | SI-rs     | SI-Eytz-rs | ImplicitITree-C++ | IntervalTree-C++ | NCLS-C   | SI-C++     | SI-Eytz-C++ |
| --------------------- | -------- | ---------- |-----------|------------| ----------------- | ---------------- | -------- |------------|-------------|
| DB53 reads, ONT reads | 2939.6   | 4746.6     | 1323      | **1273**   | 6654.6            | 46771.8          | 12082.4  | 2544.4     | **2180.2**  |
| DB53 reads, genes     | 75.2     | 131        | 26.6      | **26**     | 168.2             | 122.8            | 308.2    | 56.4       | **51.4**    |
| ONT reads, DB53 reads | 17100.6  | 19309.2    | 3815      | **3714.6** | 40490.8           | 28633.2          | 55317.6  | 24047      | **23664**   |
| anno, rna             | 89.6     | 110        | 42.2      | **41.8**   | 188.8             | 150.2            | 299.4    | **58**     | **58**      |
| genes, DB53 reads     | 2217.6   | 2448.8     | 1343.8    | **1331.6** | 4495.8            | 2747.2           | 3632.2   | **1265.2** | 1730.8      |
| mito-b, mito-a        | 39002.8  | 88901.8    | **13540** | 13541.8    | 128507.2          | 120712           | 261409.2 | 43682      | **42576.8** |
| rna, anno             | 51       | 69.2       | 12        | **11.8**   | 140.4             | 84.4             | 323.8    | 54.2       | **53**      |

### 2. Counting interval intersections

#### Intervals in sorted order

|                       | Coitrees | SI-rs     | SI-Eytz-rs | SI-C++    | SI-Eytz-C++ |
| --------------------- | -------- |-----------|------------|-----------|-------------|
| DB53 reads, ONT reads | 551.4    | 337.6     | 338        | **239.4** | 265         |
| DB53 reads, genes     | 26       | 10.6      | 10.8       | 8         | **7**       |
| ONT reads, DB53 reads | 2517.2   | **795.4** | 796.6      | 2234.2    | 1414.2      |
| anno, rna             | 26.8     | 13.4      | 13.2       | 22.6      | **12**      |
| genes, DB53 reads     | 737.4    | **292.6** | 294.6      | 459.6     | 338.2       |
| mito-b, mito-a        | 7030     | 6634.6    | 6633.4     | 3065.6    | **2991.8**  |
| rna, anno             | 9        | **4**     | **4**      | 12        | 10          |

#### Intervals in random order

|                       | Coitrees | SI-rs  | SI-Eytz-rs | SI-C++ | SI-Eytz-C++ |
| --------------------- | -------- | ------ | ---------- | ------ | ----------- |
| DB53 reads, ONT reads | 1990     | 937.2  | 883.4      | 1018.8 | **789.6**       |
| DB53 reads, genes     | 49.2     | 16     | 15         | 15.2   | **13.4**        |
| ONT reads, DB53 reads | 6835     | 4037.8 | **3964.4**     | 8547.8 | 10153.8     |
| anno, rna             | 52       | 39     | **38.6**       | 47     | 46          |
| genes, DB53 reads     | 1523.6   | 1261   | 1269       | **1119.4** | 1519.6      |
| mito-b, mito-a        | 15001.2  | 7290.6 | 7298.4     | 4493.6 | **4452.4**      |
| rna, anno             | 22       | **12**     | **12**         | 25.2   | 25.4        |

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