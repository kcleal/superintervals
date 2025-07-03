SuperIntervals
==============

A fast, memory-efficient data structure for interval intersection queries.
SuperIntervals uses a novel superset-index approach that maintains 
intervals in position-sorted order, enabling cache-friendly searches and SIMD-optimized counting.

Available for [C++](#cpp), [Rust](#rust), [Python](#python), [R](#r).

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

## Python

Install using `pip install superintervals`

```python
from superintervals import IntervalMap

imap = IntervalMap()
imap.add(10, 20, 'A')
imap.build()
results = imap.search_values(8, 20)  # ['A']
```

### API Reference

**IntervalMap** class:

- `add(start, end, value=None)`  
  Add interval with associated value


- `build()`  
  Build index (required before queries)


- `clear()`  
  Remove all intervals


- `reserve(n)`  
  Reserve space for n intervals


- `size()`  
  Get number of intervals


- `at(index)`  
  Get interval at index as (start, end, value)


- `starts_at(index)`  
  Get start position at index


- `ends_at(index)`  
  Get end position at index


- `data_at(index)`  
  Get value at index


- `has_overlaps(start, end)`  
  Check if any intervals overlap range


- `count(start, end)`  
  Count overlapping intervals


- `search_values(start, end)`  
  Get values of overlapping intervals


- `search_idxs(start, end)`  
  Get indices of overlapping intervals


- `search_keys(start, end)`  
  Get (start, end) pairs of overlapping intervals


- `search_items(start, end)`  
  Get (start, end, value) tuples of overlapping intervals


- `coverage(start, end)`  
  Get (count, total_coverage) for range

## R

```r
library(superintervals)

imap <- IntervalMap()
add(imap, 10, 20, "A")
build(imap)
results <- search_values(imap, 8, 20)
```

### API Reference

**IntervalMap()** - Create new IntervalMap object

- `add(x, start, end, value=NULL)`  
  Add interval with associated value


- `build(x)`  
  Build index (required before queries)


- `clear(x)`  
  Remove all intervals


- `reserve(x, n)`  
  Reserve space for n intervals


- `size(x)`  
  Get number of intervals


- `at(x, idx)`  
  Get interval at index as list(start, end, value)


- `starts_at(x, idx)`  
  Get start position at index


- `ends_at(x, idx)`  
  Get end position at index


- `data_at(x, idx)`  
  Get value at index


- `has_overlaps(x, start, end)`  
  Check if any intervals overlap range


- `count(x, start, end)`  
  Count overlapping intervals


- `search_values(x, start, end)`  
  Get values of overlapping intervals


- `search_idxs(x, start, end)`  
  Get indices of overlapping intervals


- `search_keys(x, start, end)`  
  Get list(start, end) of overlapping intervals


- `search_items(x, start, end)`  
  Get list(start, end, value) of overlapping intervals


- `coverage(x, start, end)`  
  Get list(count, total_coverage) for range


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

### API Reference

**IntervalMap<S, T>** class (also **IntervalMapEytz<S, T>** for Eytzinger layout):

- `void add(S start, S end, const T& value)`  
  Add interval with associated value


- `void build()`  
  Build index (required before queries)


- `void clear()`  
  Remove all intervals


- `void reserve(size_t n)`  
  Reserve space for n intervals


- `size_t size()`  
  Get number of intervals


- `const Interval<S,T>& at(size_t index)`  
  Get Interval<S,T> at index


- `void at(size_t index, Interval<S,T>& interval)`  
  Fill provided Interval object


- `bool has_overlaps(S start, S end)`  
  Check if any intervals overlap range


- `size_t count(S start, S end)`  
  Count overlapping intervals (SIMD optimized)


- `size_t count_linear(S start, S end)`  
  Count overlapping intervals (linear)


- `size_t count_large(S start, S end)`  
  Count optimized for large ranges


- `void search_values(S start, S end, std::vector<T>& found)`  
  Fill vector with values of overlapping intervals


- `void search_values_large(S start, S end, std::vector<T>& found)`  
  Search optimized for large ranges


- `void search_idxs(S start, S end, std::vector<size_t>& found)`  
  Fill vector with indices of overlapping intervals


- `void search_keys(S start, S end, std::vector<std::pair<S,S>>& found)`  
  Fill vector with (start,end) pairs


- `void search_items(S start, S end, std::vector<Interval<S,T>>& found)`  
  Fill vector with Interval<S,T> objects


- `void search_point(S point, std::vector<T>& found)`  
  Find intervals containing single point


- `void coverage(S start, S end, std::pair<size_t,S>& result)`  
  Get pair(count, total_coverage) for range


- `IndexRange search_idxs(S start, S end)`  
  Returns IndexRange for range-based loops over indices


- `ItemRange search_items(S start, S end)`  
  Returns ItemRange for range-based loops over intervals


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

### API Reference

**IntervalMap<T>** struct (also **IntervalMapEytz<T>** for Eytzinger layout):

- `fn new() -> Self`  
  Create new IntervalMap


- `fn add(&mut self, start: i32, end: i32, value: T)`  
  Add interval with associated value


- `fn build(&mut self)`  
  Build index (required before queries)


- `fn clear(&mut self)`  
  Remove all intervals


- `fn reserve(&mut self, n: usize)`  
  Reserve space for n intervals


- `fn size(&self) -> usize`  
  Get number of intervals


- `fn at(&self, index: usize) -> Interval<T>`  
  Get Interval<T> at index


- `fn has_overlaps(&mut self, start: i32, end: i32) -> bool`  
  Check if any intervals overlap range


- `fn count(&mut self, start: i32, end: i32) -> usize`  
  Count overlapping intervals (SIMD optimized)


- `fn count_linear(&mut self, start: i32, end: i32) -> usize`  
  Count overlapping intervals (linear)


- `fn count_large(&mut self, start: i32, end: i32) -> usize`  
  Count optimized for large ranges


- `fn search_values(&mut self, start: i32, end: i32, found: &mut Vec<T>)`  
  Fill vector with values of overlapping intervals


- `fn search_values_large(&mut self, start: i32, end: i32, found: &mut Vec<T>)`  
  Search optimized for large ranges


- `fn search_idxs(&mut self, start: i32, end: i32, found: &mut Vec<usize>)`  
  Fill vector with indices of overlapping intervals


- `fn search_keys(&mut self, start: i32, end: i32, found: &mut Vec<(i32, i32)>)`  
  Fill vector with (start,end) pairs


- `fn search_items(&mut self, start: i32, end: i32, found: &mut Vec<Interval<T>>)`  
  Fill vector with Interval<T> objects


- `fn search_stabbed(&mut self, point: i32, found: &mut Vec<T>)`  
  Find intervals containing single point


- `fn coverage(&mut self, start: i32, end: i32) -> (usize, i32)`  
  Get (count, total_coverage) for range


- `fn search_idxs_iter(&mut self, start: i32, end: i32) -> IndexIterator<T>`  
  Iterator over indices


- `fn search_items_iter(&mut self, start: i32, end: i32) -> ItemIterator<T>`  
  Iterator over intervals


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

## Benchmark

SuperIntervals (SI) was compared with:

- Coitrees (Rust: https://github.com/dcjones/coitrees)
- Implicit Interval Tree (C++: https://github.com/lh3/cgranges)
- Interval Tree (C++: https://github.com/ekg/intervaltree)
- Nested Containment List (C: https://github.com/pyranges/ncls/tree/master/ncls/src)

Main results:

- Finding interval intersections is on average ~1.5-3x faster than other libraries (Coitrees for Rust, Implicit Interval Tree for C++), with some 
exceptions. Coitrees-s was faster for one test (ONT reads, sorted DB53 reads).
- The SIMD counting performance of coitrees and superintervals was similar.

Datasets https://github.com/kcleal/superintervals/releases/download/v0.2.0/data.tar.gz:

- `rna / anno` RNA-seq reads and annotations from cgranges repository
- `ONT reads` nanopore alignments from sample PAO33946 chr1, converted to bed format
- `DB53 reads` paired-end reads from sample DB53, NCBI BioProject PRJNA417592, chr1, converted to bed format
- `mito-b, mito-a` paired-end reads from sample DB53 chrM, converted to bed format (mito-b and mito-a are the same)
- `genes` UCSC genes from hg19

Test programs use internal timers and print data to stdout, measuring the index time, and time to find all intersections. Other steps such as file IO are ignored. Test programs also only assess chr1 bed records - other chromosomes are ignored. For 'chrM' records, the M was replaced with 1 using sed. Data were assessed in position sorted and random order. Datasets can be found on the Releases page, and the test/run_tools.sh script has instructions for how to repeat the benchmark.

Timings were in microseconds using an i9-11900K, 64 GB, 2TB NVMe machine.

## Finding interval intersections

- Coitrees-s uses the SortedQuerent version of coitrees
- SI = superintervals. Eytz refers to the eytzinger layout. -rs is the Rust implementation.

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


## Acknowledgements

- The rust test program borrows heavily from the coitrees package
- The superset-index implemented here exploits a similar interval ordering as described in
Schmidt 2009 "Interval Stabbing Problems in Small Integer Ranges". However, the superset-index has several advantages including
  1. An implicit memory layout
  1. General purpose implementation (not just small integer ranges)
  1. SIMD counting algorithm 
- The Eytzinger layout was adapted from Sergey Slotin, Algorithmica
