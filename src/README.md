# Superintervals c++/rust API

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
