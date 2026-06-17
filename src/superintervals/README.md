# Superintervals Python API

Install using `pip install superintervals`

```python
from superintervals import IntervalMap
from array import array

# Method 1: Manual construction
imap = IntervalMap()
imap.add(10, 20, 'A')
imap.add(15, 25, 'B')
imap.build()
results = imap.search_values(8, 20)  # ['A', 'B']

# Method 2: Efficient construction from arrays (no build() needed)
starts = array('i', [10, 15, 30])
ends = array('i', [20, 25, 40])
values = ['A', 'B', 'C']
imap = IntervalMap.from_arrays(starts, ends, values)
results = imap.search_values(8, 20)  # ['A', 'B']

# Batch operations for high performance
query_starts = array('i', [5, 18, 35])
query_ends = array('i', [12, 22, 45])
counts = imap.count_batch(query_starts, query_ends)  # [1, 2, 1]
indices = imap.search_idxs_batch(query_starts, query_ends)  # [[0], [0, 1], [2]]
```

### API Reference

**IntervalMap** class:

#### Construction Methods
- `IntervalMap()`  
  Create empty interval map

- `IntervalMap.from_arrays(starts, ends, values=None)`  
  Create from arrays (array.array, numpy arrays, or lists)  
  Returns ready-to-use IntervalMap (no build() needed)

#### Adding Intervals
- `add(start, end, value=None)`  
  Add interval with associated value

- `build()`  
  Build index (required before queries when using add())

#### Memory Management
- `clear()`  
  Remove all intervals

- `reserve(n)`  
  Reserve space for n intervals

- `size()`  
  Get number of intervals

#### Access Methods
- `at(index)`  
  Get interval at index as (start, end, value)

- `starts_at(index)`  
  Get start position at index

- `ends_at(index)`  
  Get end position at index

- `data_at(index)`  
  Get value at index

#### Single Query Methods
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

#### Batch Query Methods (High Performance)
- `count_batch(starts, ends)`  
  Count overlaps for multiple ranges  
  Args: Memory views (array.array, numpy arrays)  
  Returns: List of counts

- `search_idxs_batch(starts, ends)`  
  Get indices for multiple ranges  
  Args: Memory views (array.array, numpy arrays)  
  Returns: List of lists containing indices

- `search_values_batch(starts, ends)`  
  Get values for multiple ranges  
  Args: Memory views (array.array, numpy arrays)  
  Returns: List of lists containing values

#### Set Operations

These build new disjoint interval maps from one or two existing maps. Unlike the C++
and Rust APIs, the returned map is already built and ready to query. When intervals
merge, the surviving value defaults to a **tuple** of the merged values (since Python
data is opaque); pass a `combine(a, b)` callable to override this.

```python
a = IntervalMap()
a.add(1, 5, 'x'); a.add(3, 8, 'y'); a.build()
b = IntervalMap()
b.add(4, 10, 'B'); b.build()

m = a.merge_overlaps()              # one interval (1, 8); value ('x', 'y')
m = a.merge_overlaps(lambda p, q: p)  # custom combiner -> keep first
u = a.union_with(b)                 # a ∪ b
i = a.intersection(b)               # a ∩ b; values paired (a_val, b_val) by default
d = a.difference(b)                 # a \ b
x = a.symmetric_difference(b)       # in exactly one of a, b
g = a.gaps(0, 100)                  # uncovered regions within [0, 100]
sp = a.span()                       # (min_start, max_end) or None if empty
```

- `merge_overlaps(combine=None)`  
  Collapse self-overlapping intervals; default merges values into a tuple

- `union_with(other, combine=None)`  
  Union of the two maps

- `intersection(other, combine=None)`  
  Intersection; default pairs overlapping values as `(a_val, b_val)`

- `difference(other)`  
  `self \ other`

- `symmetric_difference(other)`  
  Regions in exactly one of the two maps

- `gaps(lo, hi, fill=None)`  
  Uncovered regions within `[lo, hi]`

- `span()`  
  `(min_start, max_end)` tuple, or `None` if empty

### Performance Tips

- Use `IntervalMap.from_arrays()` for best construction performance
- Use batch methods for multiple queries (often 5-10x faster)
- Convert lists to arrays for batch operations:
  ```python
  from array import array
  starts = array('i', [1, 5, 10])  # For batch methods
  ```
- Reserve space with `reserve(n)` when adding many intervals manually