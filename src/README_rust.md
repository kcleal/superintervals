# Superintervals rust API

Add to your project using `cargo add superintervals`

```rust
use superintervals::IntervalMap;

let mut imap = IntervalMap::new();
imap.add(1, 5, "A");
imap.build();

// Collect results into a Vec
let mut results = Vec::new();
imap.search_values(4, 11, &mut results);

// Or use lazy iterator interfaces
for value in imap.search_values_iter(query_start, query_end) {
    println!("Found: {}", value);
}
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

- `fn search_keys_iter(&mut self, start: i32, end: i32) -> KeyIterator<T>`  
  Iterator over keys

- `fn search_values_iter(&mut self, start: i32, end: i32) -> ValueIterator<T>`  
  Iterator over values