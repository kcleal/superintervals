# Superintervals C++ API


Header only implementation, copy to your include directory.

```cpp
#include "SuperIntervals.hpp"

si::IntervalMap<int, std::string> imap;
imap.add(1, 5, "A");
imap.build();

// Collect results into a vector
std::vector<std::string> results;
imap.search_values(4, 9, results);

// Or use lazy iterator interfaces
for (const auto &value : imap.search_values_iter(query_start, query_end)) {
    std::cout << "Found: " << value << std::endl;
}
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