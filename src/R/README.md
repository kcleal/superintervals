# SuperIntervals

The R Bioconductor package `superintervalsr` provides a fast, memory-efficient data structure for interval intersection queries in R.
Built on a novel superset-index approach that maintains intervals in position-sorted order,
enabling cache-friendly searches and SIMD-optimized counting operations.

## Key Features

- Linear-time index construction from sorted intervals
- Cache-friendly querying with minimal memory overhead
- SIMD acceleration (AVX2/Neon) for counting operations
- High performance - often faster than existing interval libraries
- Minimal memory overhead - one size_t per interval
- End-inclusive intervals - intervals include both start and end positions

## Installation

Build and install:
```r
bash ./build.sh

# Run a small benchmark vs IRanges
Rscript benchmark.R
```

## Quick Start

### Creating an IntervalMap

The core data structure is `IntervalMap`, which stores intervals with associated values:

```r
# Create interval map
imap <- IntervalMap()

# Add genomic intervals (e.g., genes) - intervals are end-inclusive
add(imap, 1000, 5000, list(gene_id="ENSG001", name="GENE_A"))
add(imap, 3000, 8000, list(gene_id="ENSG002", name="GENE_B")) 
add(imap, 10000, 15000, list(gene_id="ENSG003", name="GENE_C"))

# The build() function must be called before any queries
build(imap)

# Check basic properties
cat("Number of intervals:", length(imap), "\n")
cat("Total intervals:", size(imap), "\n")
```

### Accessing Intervals

```r
# Access intervals by index (1-based in R)
interval1 <- imap[1]
print(interval1)

# Alternative syntax
interval2 <- at(imap, 2)
print(interval2)

# Get specific components
start_pos <- starts_at(imap, 1)
end_pos <- ends_at(imap, 1)
value <- data_at(imap, 1)
```

### Querying Overlaps

SuperIntervals provides multiple ways to query interval overlaps:

```r
# Define query region
query_start <- 4000
query_end <- 6000

# Check for any overlaps
has_overlaps_result <- has_overlaps(imap, query_start, query_end)
cat("Has overlaps:", has_overlaps_result, "\n")

# Count overlapping intervals (SIMD-optimized)
overlap_count <- count(imap, query_start, query_end)
cat("Number of overlaps:", overlap_count, "\n")

# Find indices of overlapping intervals
indices <- search_idxs(imap, query_start, query_end)
print("Overlapping interval indices:")
print(indices)

# Find overlapping interval coordinates
keys <- search_keys(imap, query_start, query_end)
print("Overlapping coordinates:")
print(keys)

# Find values associated with overlapping intervals
values <- search_values(imap, query_start, query_end)
print("Associated values:")
print(values)

# Find complete overlapping intervals (coordinates + values)
overlaps <- search_items(imap, query_start, query_end)
print("Complete overlaps:")
print(overlaps)
```

### Coverage Analysis

```r
# Get coverage statistics for a range
coverage_stats <- coverage(imap, 1000, 15000)
print("Coverage statistics:")
print(coverage_stats)
```

### Memory Management

```r
# Reserve space for better performance when adding many intervals
reserve(imap, 1000)  # Reserve space for 1000 intervals

# Clear all intervals
clear(imap)
```

## Example: Genomic Analysis

```r
library(superintervalsr)

# Create interval map for gene annotations
genes <- IntervalMap()

# Add gene intervals
add(genes, 1000, 5000, "GENE_A")
add(genes, 3000, 8000, "GENE_B")
add(genes, 10000, 15000, "GENE_C")
add(genes, 12000, 18000, "GENE_D")

# Build index for efficient querying
build(genes)

# Query for genes overlapping a region of interest
roi_start <- 4000
roi_end <- 12500

# Find all overlapping genes
overlapping_genes <- search_values(genes, roi_start, roi_end)
cat("Genes overlapping region", roi_start, "-", roi_end, ":\n")
print(unlist(overlapping_genes))

# Get detailed overlap information
overlap_details <- search_items(genes, roi_start, roi_end)
for (i in seq_along(overlap_details$start)) {
  cat(sprintf("Gene: %s, Position: %d-%d\n", 
              overlap_details$value[[i]], 
              overlap_details$start[i], 
              overlap_details$end[i]))
}
```

## Performance Notes

- Call `build()` after adding all intervals and before performing queries
- Use `reserve()` when you know the approximate number of intervals to improve performance
- The package is optimized for scenarios with many intervals and frequent queries
- SIMD optimizations are automatically used when available on your system

## API Reference

### Core Functions
- `IntervalMap()` - Create a new interval map
- `add(x, start, end, value)` - Add an interval with optional value
- `build(x)` - Build index for efficient querying
- `clear(x)` - Remove all intervals
- `reserve(x, n)` - Reserve space for n intervals

### Query Functions
- `has_overlaps(x, start, end)` - Check if any intervals overlap
- `count(x, start, end)` - Count overlapping intervals
- `search_idxs(x, start, end)` - Get indices of overlapping intervals
- `search_keys(x, start, end)` - Get coordinates of overlapping intervals  
- `search_values(x, start, end)` - Get values of overlapping intervals
- `search_items(x, start, end)` - Get complete overlap information
- `coverage(x, start, end)` - Get coverage statistics

### Access Functions
- `length(x)` - Number of intervals
- `size(x)` - Number of intervals (alias)
- `at(x, idx)` - Get interval at index
- `starts_at(x, idx)` - Get start position at index
- `ends_at(x, idx)` - Get end position at index  
- `data_at(x, idx)` - Get value at index