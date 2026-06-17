# Superintervals C++ API

Superintervals implements `IntervalMap` which is a multimap like data structure for 
interval intersection queries.
'Keys' correspond to `start` and `end` coordinates of reference intervals 
of type `S`. 'Values' are your data of interest of type `T`.

To install, copy the superintervals.hpp header to your include directory.


```cpp
#include "SuperIntervals.hpp"

si::IntervalMap<int, std::string> imap;
imap.add(1, 5, "A");
imap.build();

// Collect results into a vector
std::vector<std::string> results;
imap.search_values(4, 9, results);

// Or use lazy iterator interfaces
for (const auto &value : imap.search_values(query_start, query_end)) {
    std::cout << "Found: " << value << std::endl;
}
// also search_keys, search_idxs, search_items
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


- `KeyRange search_keys(S start, S end)`  
  Returns KeyRange for range-based loops over intervals

- `ValueRange search_values(S start, S end)`  
  Returns ValueRange for range-based loops over intervals

### Set Operations

These build new disjoint interval sets from one or two maps. Results are returned
**unindexed** — call `build()` on the result before querying it (or before passing it
as `other` to another set operation). The folding operations accept an optional
`combine` callable `T(const T& a, const T& b)` to decide which value survives when
intervals merge; the default keeps the first. `gaps`, `difference`, and
`symmetric_difference` use `±1` and therefore require an integral `S`.

```cpp
si::IntervalMap<int, std::string> a, b;
// ... add() intervals to both ...
a.build(); b.build();

auto merged = a.merge_overlaps();   // collapse self-overlaps
auto u      = a.union_with(b);      // a ∪ b
auto i      = a.intersection(b);    // a ∩ b
auto d      = a.difference(b);      // a \ b
merged.build();                     // build() before querying a result
```

- `IntervalMap<S,T> merge_overlaps(Combine combine = keep_first)`  
  Collapse self-overlapping intervals into a disjoint set

- `IntervalMap<S,T> gaps(S lo, S hi, const T& fill = T{})`  
  Uncovered regions within `[lo, hi]` (integral `S`)

- `IntervalMap<S,T> union_with(const IntervalMap<S,T>& other, Combine combine = keep_first)`  
  Union of the two sets

- `IntervalMap<S,T> intersection(const IntervalMap<S,T>& other, Combine combine = keep_first)`  
  Intersection of the two sets

- `IntervalMap<S,T> difference(const IntervalMap<S,T>& other)`  
  `this \ other` (integral `S`)

- `IntervalMap<S,T> symmetric_difference(const IntervalMap<S,T>& other)`  
  Regions in exactly one of the two sets (integral `S`)

- `bool span(std::pair<S,S>& result)`  
  Fill `result` with `(min start, max end)`; returns `false` if empty