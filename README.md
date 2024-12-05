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

Benchmark

SuperIntervals (SI) was compared with:

    Coitrees (Rust: https://github.com/dcjones/coitrees)
    Implicit Interval Tree (C++: https://github.com/lh3/cgranges)
    Interval Tree (C++: https://github.com/ekg/intervaltree)
    Nested Containment List (C: https://github.com/pyranges/ncls/tree/master/ncls/src)

Main results:

- Finding interval intersections is roughly ~1.5-3x faster than the next best library (Coitrees for Rust, Implicit Interval Tree for C++), with some 
exceptions. Coitrees-s was faster for one test (ONT reads, sorted DB53 reads).
- The SIMD counting performance of coitrees and superintervals is similar.

Datasets:

    1. `rna / anno` RNA-seq reads and annotations from cgranges repository
    2. `ONT reads` nanopore alignments from sample PAO33946 chr1, converted to bed format
    3. `DB53 reads` paired-end reads from sample DB53, NCBI BioProject PRJNA417592, chr1, converted to bed format
    4. `mito-b, mito-a` paired-end reads from sample DB53 chrM, converted to bed format (mito-b and mito-a are the same)
    5. `genes` UCSC genes from hg19

Test programs use internal timers and print data to stdout, measuring the index time, and time to find all intersections. Other steps such as file IO are ignored. Test programs also only assess chr1 bed records - other chromosomes are ignored. For 'chrM' records, the M was replaced with 1 using sed. Data were assessed in position sorted and random order. Datasets can be found on the Releases page, and the test/run_tools.sh script has instructions for how to repeat the benchmark.

Timings were in microseconds using an i9-11900K, 64 GB, 2TB NVMe machine.
## Finding interval intersections

    Coitrees-s uses the SortedQuerent version of coitrees
    SI = superintervals. Eytz refers to the eytzinger layout. -rs is the Rust implementation.

### Intervals in sorted order

|                       | Coitrees | Coitrees-s | SuperIntervals-rs | SuperIntervalsEytz-rs | ImplicitITree-C++ | IntervalTree-C++ | NCLS-C | SuperIntervals-C++ | SuperIntervalsEytz-C++ |
| --------------------- | -------- | ---------- |-------------------| --------------------- | ----------------- | ---------------- | ------ | ------------------ | ---------------------- |
| DB53 reads, ONT reads | 1668     | 3179       | **757**             | **757**                   | 3831              | 44404            | 10642  | **1315**               | 1358                   |
| DB53 reads, genes     | 55       | 84         | **21**                | **21**                    | 122               | 109              | 291    | 42                 | **40**                     |
| ONT reads, DB53 reads | 6504     | **3354**       | 3859              | 3854                  | 17949             | 12280            | 30772  | 5290               | **4462**                   |
| anno, rna             | 50       | 35         | **18**                | **18**                    | 127               | 90               | 208    | 29                 | **22**                     |
| genes, DB53 reads     | 1171     | 1018       | 301               | **296**                   | 3129              | 1315             | 1780   | 442                | **323**                    |
| mito-b, mito-a        | 34769    | 34594      | 16971             | **16952**                 | 93900             | 107660           | 251707 | 33177              | **32985**                  |
| rna, anno             | 31       | 23         | 21                | **20**                    | 70                | 55               | 233    | 28                 | **27**                     |

### Intervals in random order

|                       | Coitrees | Coitrees-s | SuperIntervals-rs | SuperIntervalsEytz-rs | ImplicitITree-C++ | IntervalTree-C++ | NCLS-C | SuperIntervals-C++ | SuperIntervalsEytz-C++ |
| --------------------- | -------- | ---------- | ----------------- | --------------------- | ----------------- | ---------------- | ------ | ------------------ | ---------------------- |
| DB53 reads, ONT reads | 2943     | 4663       | 1356              | **1355**                  | 6505              | 46743            | 11947  | 2491               | **2169**                   |
| DB53 reads, genes     | 78       | 130        | 27                | **26**                    | 170               | 125              | 305    | 58                 | **51**                     |
| ONT reads, DB53 reads | 16650    | 18931      | 16116             | **16037**                 | 38677             | 27832            | 53452  | **23003**              | 23232                  |
| anno, rna             | 89       | 105        | **54**                | **54**                    | 188               | 143              | 294    | **58**                 | 60                     |
| genes, DB53 reads     | 2222     | 2424       | 1693              | **1684**                  | 4490              | 2701             | 3605   | **1251**               | 1749                   |
| mito-b, mito-a        | 38030    | 86309      | **18326**             | 18368                 | 125336            | 118321           | 256293 | 42195              | **41695**                  |
| rna, anno             | 53       | 73         | **45**                | **45**                    | 137               | 83               | 311    | **52**                 | **52**                     |

## Counting interval intersections

### Intervals in sorted order

|                       | Coitrees | SuperIntervals-rs | SuperIntervalsEytz-rs | SuperIntervals-C++ | SuperIntervalsEytz-C++ |
| --------------------- | -------- | ----------------- | --------------------- | ------------------ | ---------------------- |
| DB53 reads, ONT reads | 551      | 370               | 371                   | **241**                | 263                    |
| DB53 reads, genes     | 28       | 12                | 12                    | 8                  | **7**                      |
| ONT reads, DB53 reads | 2478     | 1909              | 1890                  | 2209               | **1312**                   |
| anno, rna             | 26       | 14                | 14                    | 22                 | **11**                     |
| genes, DB53 reads     | 747      | 321               | 336                   | 446                | **290**                    |
| mito-b, mito-a        | 6894     | 6727              | 6746                  | 3088               | **2966**                   |
| rna, anno             | **9**        | 13                | 13                    | 12                 | 10                     |

### Intervals in random order

|                       | Coitrees | SuperIntervals-rs | SuperIntervalsEytz-rs | SuperIntervals-C++ | SuperIntervalsEytz-C++ |
| --------------------- | -------- | ----------------- | --------------------- | ------------------ | ---------------------- |
| DB53 reads, ONT reads | 1988     | 972               | 969                   | 1016               | **778**                    |
| DB53 reads, genes     | 53       | 20                | 20                    | 16                 | **13**                     |
| ONT reads, DB53 reads | 6692     | 8864              | 8733                  | **8182**               | 9523                   |
| anno, rna             | 52       | 49                | 48                    | **47**                 | 50                     |
| genes, DB53 reads     | 1503     | 1628              | 1592                  | **1120**               | 1623                   |
| mito-b, mito-a        | 14354    | 7579              | 7600                  | 4442               | **4383**                   |
| rna, anno             | 22       | 30                | 29                    | **25**                 | **25**                     |

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
- The Eytzinger layout was adapted from Sergey Slotin, Algorithmica
