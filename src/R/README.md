# SuperIntervals

The R Bioconductor package `superintervalsr` provides very fast interval intersection queries in R.
Built on a novel superset-index approach with batch operations in mind.

## Why SuperIntervals?

- Faster interval intersection queries
- Low memory overhead
- Bioconductor friendly: GRanges → IntervalMap → query → back to GRanges

## Installation

```r
remotes::install_github("kcleal/superintervalsr")
```

## Core Workflow

### 1. Create IntervalMap from Your Data

```r
# From vectors (recommended - no build() needed)
imap <- IntervalMap.from_vectors(starts, ends, values)

# From GRanges (Bioconductor integration)
library(GenomicRanges)
gr <- GRanges("chr1", IRanges(starts, ends), gene_id = values)
imap <- IntervalMap(gr, value_column = "gene_id", seqname = "chr1")

# From data.frame
df <- data.frame(start = starts, end = ends, gene = values)
imap <- IntervalMap(df, value_column = "gene")
```

### 2. Batch Query Operations (The Fast Way)

```r
# Multiple query regions
query_starts <- c(2000, 7000, 12000, 50000)
query_ends <- c(4000, 9000, 14000, 55000)

# Batch counting - orders of magnitude faster
counts <- count.batch(imap, query_starts, query_ends)
#> [1] 3 2 1 4

# Batch index search - get all overlapping interval indices  
indices_list <- search_idxs.batch(imap, query_starts, query_ends)
#> [[1]]
#> [1] 12 15 23
#> [[2]] 
#> [1] 45 67
#> ... etc

# Use indices to get back to original data
for (i in seq_along(indices_list)) {
  if (length(indices_list[[i]]) > 0) {
    overlapping_features <- values[indices_list[[i]]]
    cat("Query", i, "overlaps:", paste(overlapping_features, collapse = ", "), "\n")
  }
}
```

### 3. Single Query Operations (When You Need Just One)

```r
# For occasional single queries
single_count <- count(imap, 4000, 6000)
single_indices <- search_idxs(imap, 4000, 6000)
single_values <- search_values(imap, 4000, 6000)

# Get complete overlap information
overlaps <- search_items(imap, 4000, 6000)
# Returns: list(start = c(...), end = c(...), value = list(...))
```

## Real-World Genomics Example

```r
library(superintervalsr)
library(GenomicRanges)

# Load your gene annotations
genes <- GRanges("chr1", IRanges(c(1000, 5000, 10000), c(4000, 8000, 15000)),
                 gene_id = c("ENSG001", "ENSG002", "ENSG003"),
                 gene_name = c("GENE_A", "GENE_B", "GENE_C"))

# Convert to IntervalMap for fast queries
gene_imap <- IntervalMap(genes, value_column = "gene_id", seqname = "chr1")

# You have thousands of ChIP-seq peaks to annotate
n_peaks <- 5000
peak_starts <- sort(sample(1:20000, n_peaks))
peak_ends <- peak_starts + sample(200:1000, n_peaks, replace = TRUE)

# Batch annotate all peaks (lightning fast!)
peak_gene_counts <- count.batch(gene_imap, peak_starts, peak_ends)
peak_gene_indices <- search_idxs.batch(gene_imap, peak_starts, peak_ends)

# Analyze results
peaks_with_genes <- sum(peak_gene_counts > 0)
cat("Peaks overlapping genes:", peaks_with_genes, "out of", n_peaks, "\n")

# Get detailed annotations for peaks with genes
for (i in which(peak_gene_counts > 0)[1:5]) {  # Show first 5
  overlapping_genes <- genes[peak_gene_indices[[i]]]
  cat("Peak", i, "overlaps:", 
      paste(mcols(overlapping_genes)$gene_name, collapse = ", "), "\n")
}
```

## Performance Comparison

```
Rscript ./benchmark.R

>>> 4. SUMMARY
----------
Performance vs IRanges (NCLS algorithm):
- Overlap counting: 48.9x faster
- Overlap finding:  10.6x faster
```


## Complete API Reference

### Core Data Structure
```r
# Create IntervalMap
IntervalMap()                                     # Empty map
IntervalMap.from_vectors(starts, ends, values)    # From vectors (recommended)
IntervalMap(granges, value_column, seqname)       # From GRanges
IntervalMap(dataframe, start_column, end_column)  # From data.frame
```

### Batch Operations (Use These!)
```r
count.batch(imap, starts, ends)        # Count overlaps for multiple queries
search_idxs.batch(imap, starts, ends)  # Get indices for multiple queries
```

### Single Query Operations
```r
has_overlaps(imap, start, end)    # TRUE/FALSE for any overlap
count(imap, start, end)           # Count overlapping intervals
search_idxs(imap, start, end)     # Get indices of overlapping intervals
search_values(imap, start, end)   # Get values of overlapping intervals
search_keys(imap, start, end)     # Get coordinates of overlapping intervals
search_items(imap, start, end)    # Get complete overlap information (coords + values)
coverage(imap, start, end)        # Coverage statistics
```

### Bioconductor Integration
```r
as_granges(search_items_result, seqname, valuename)  # Convert back to GRanges
```

### Utilities
```r
length(imap)           # Number of intervals
at(imap, index)        # Get interval at specific index
clear(imap)            # Remove all intervals
reserve(imap, n)       # Pre-allocate space for n intervals
```

## Advanced: Manual Interval Addition

```r
# For dynamic construction (less common)
imap <- IntervalMap()
add(imap, 1000, 5000, "interval1")
add(imap, 3000, 8000, "interval2")
build(imap)  # Required before queries when using add()

# Batch operations work the same way
counts <- count.batch(imap, query_starts, query_ends)
```

## Key Technical Details

- **End-inclusive intervals**: Both start and end positions are included
- **Automatic sorting**: Intervals are maintained in sorted order
- **SIMD optimization**: Automatic vectorization on modern CPUs  
- **Memory efficient**: One size_t per interval overhead
- **Thread safe**: Read operations can be parallelized
- **Cache friendly**: Optimized memory layout for fast iteration