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
- Batch operations for efficient bulk queries

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
# Method 1: Create empty and add intervals
imap <- IntervalMap()

# Add genomic intervals (e.g., genes) - intervals are end-inclusive
add(imap, 1000, 5000, list(gene_id="ENSG001", name="GENE_A"))
add(imap, 3000, 8000, list(gene_id="ENSG002", name="GENE_B")) 
add(imap, 10000, 15000, list(gene_id="ENSG003", name="GENE_C"))

# The build() function must be called before any queries
build(imap)

# Method 2: Create directly from vectors (no build() needed)
starts <- c(1000, 3000, 10000)
ends <- c(5000, 8000, 15000)
values <- list(
  list(gene_id="ENSG001", name="GENE_A"),
  list(gene_id="ENSG002", name="GENE_B"),
  list(gene_id="ENSG003", name="GENE_C")
)
imap2 <- IntervalMap.from_vectors(starts, ends, values)

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

### Batch Operations for High Performance

For scenarios with many queries, batch operations provide significant performance improvements:

```r
# Define multiple query regions
query_starts <- c(2000, 7000, 12000)
query_ends <- c(4000, 9000, 14000)

# Batch count overlaps - much faster than individual count() calls
batch_counts <- count.batch(imap, query_starts, query_ends)
cat("Overlap counts for each query:", batch_counts, "\n")

# Batch search for indices
batch_indices <- search_idxs.batch(imap, query_starts, query_ends)
print("Indices for each query:")
for (i in seq_along(batch_indices)) {
  cat("Query", i, ":", batch_indices[[i]], "\n")
}

# Verify consistency between batch and individual operations
for (i in seq_along(query_starts)) {
  individual_count <- count(imap, query_starts[i], query_ends[i])
  batch_count <- batch_counts[i]
  cat("Query", i, "- Individual:", individual_count, "Batch:", batch_count, "\n")
}
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

# Create interval map for gene annotations using from_vectors
gene_starts <- c(1000, 3000, 10000, 12000)
gene_ends <- c(5000, 8000, 15000, 18000)
gene_names <- c("GENE_A", "GENE_B", "GENE_C", "GENE_D")

genes <- IntervalMap.from_vectors(gene_starts, gene_ends, gene_names)

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

## Example: Batch Processing for Large-Scale Analysis

```r
library(superintervalsr)

# Create a large set of genomic features
n_features <- 10000
feature_starts <- sort(sample(1:1000000, n_features))
feature_ends <- feature_starts + sample(100:5000, n_features, replace = TRUE)
feature_ids <- paste0("feature_", 1:n_features)

# Create IntervalMap efficiently
features <- IntervalMap.from_vectors(feature_starts, feature_ends, feature_ids)

# Define many query regions (e.g., from sequencing reads)
n_queries <- 1000
query_starts <- sort(sample(1:1000000, n_queries))
query_ends <- query_starts + sample(50:500, n_queries, replace = TRUE)

# Batch processing is much faster than individual queries
system.time({
  batch_counts <- count.batch(features, query_starts, query_ends)
  batch_indices <- search_idxs.batch(features, query_starts, query_ends)
})

# Analyze results
total_overlaps <- sum(batch_counts)
queries_with_overlaps <- sum(batch_counts > 0)
cat(sprintf("Total overlaps: %d\n", total_overlaps))
cat(sprintf("Queries with overlaps: %d/%d (%.1f%%)\n", 
           queries_with_overlaps, n_queries, 
           100 * queries_with_overlaps / n_queries))
```

## Performance Notes

- Call `build()` after adding all intervals and before performing queries (not needed with `IntervalMap.from_vectors`)
- Use `IntervalMap.from_vectors()` for better performance when creating from existing data
- Use batch operations (`count.batch`, `search_idxs.batch`) for multiple queries - often 2-10x faster
- Use `reserve()` when you know the approximate number of intervals to improve performance
- The package is optimized for scenarios with many intervals and frequent queries
- SIMD optimizations are automatically used when available on your system

## API Reference

### Core Functions
- `IntervalMap()` - Create a new empty interval map
- `IntervalMap.from_vectors(starts, ends, values)` - Create interval map from vectors (ready to use)
- `add(x, start, end, value)` - Add an interval with optional value
- `build(x)` - Build index for efficient querying (required after adding intervals)
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

### Batch Query Functions
- `count.batch(x, starts, ends)` - Count overlaps for multiple query ranges
- `search_idxs.batch(x, starts, ends)` - Get indices for multiple query ranges

### Access Functions
- `length(x)` - Number of intervals
- `size(x)` - Number of intervals (alias)
- `at(x, idx)` - Get interval at index
- `starts_at(x, idx)` - Get start position at index
- `ends_at(x, idx)` - Get end position at index  
- `data_at(x, idx)` - Get value at index