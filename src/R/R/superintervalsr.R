#' @useDynLib superintervalsr, .registration=TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Create a new IntervalMap
#'
#' Creates a new IntervalMap object for storing and querying intervals efficiently.
#'
#' @return An IntervalMap object
#' @export
#' @examples
#' im <- IntervalMap()
#' add(im, 1, 10, "gene1")
#' add(im, 15, 25, "gene2")
#' build(im)
IntervalMap <- function() {
  structure(create_intervalmap(), class = "IntervalMap")
}

#' Create IntervalMap from vectors
#'
#' Creates a new IntervalMap object from vectors of start positions, end positions,
#' and optional values. The resulting IntervalMap is ready to use (no need to call build() afterwards).
#'
#' @param starts Integer vector of start positions (inclusive)
#' @param ends Integer vector of end positions (inclusive)
#' @param values Optional list of values to associate with each interval. If NULL, no values are stored.
#' @return An IntervalMap object
#' @export
#' @examples
#' # Create with values
#' im <- IntervalMap.from_vectors(c(1, 10), c(5, 15), c("gene1", "gene2"))
#'
#' # Create without values
#' im2 <- IntervalMap.from_vectors(c(1, 10), c(5, 15))
IntervalMap.from_vectors <- function(starts, ends, values = NULL) {
    if (is.null(values)) {
        values <- list()
    }
    im <- create_intervalmap_from_vectors(
        as.integer(starts),
        as.integer(ends),
        as.list(values)
    )
    structure(im, class = "IntervalMap")
}

#' Get the length of an IntervalMap
#'
#' @param x IntervalMap object
#' @return Number of intervals in the map
#' @export
length.IntervalMap <- function(x) {
  get_size(x)
}

#' Print method for IntervalMap
#'
#' @param x IntervalMap object
#' @param ... Additional arguments (ignored)
#' @export
print.IntervalMap <- function(x, ...) {
  n <- length(x)
  cat("IntervalMap with", n, "intervals\n")
  if (n > 0 && n <= 10) {
    for (i in 1:n) {
      interval <- at(x, i)
      cat(sprintf("  [%d] %d-%d", i, interval$start, interval$end))
      if (!is.null(interval$value)) {
        cat(sprintf(" -> %s", as.character(interval$value)))
      }
      cat("\n")
    }
  } else if (n > 10) {
    for (i in 1:3) {
      interval <- at(x, i)
      cat(sprintf("  [%d] %d-%d", i, interval$start, interval$end))
      if (!is.null(interval$value)) {
        cat(sprintf(" -> %s", as.character(interval$value)))
      }
      cat("\n")
    }
    cat("  ...\n")
    for (i in (n-2):n) {
      interval <- at(x, i)
      cat(sprintf("  [%d] %d-%d", i, interval$start, interval$end))
      if (!is.null(interval$value)) {
        cat(sprintf(" -> %s", as.character(interval$value)))
      }
      cat("\n")
    }
  }
  invisible(x)
}

#' Access interval at specific index
#'
#' @param x IntervalMap object
#' @param i Index (1-based in R)
#' @return List with start, end, and value components
#' @export
`[.IntervalMap` <- function(x, i) {
  at(x, i)
}

#' Add an interval to the map
#'
#' @param x IntervalMap object
#' @param start Start position (inclusive)
#' @param end End position (inclusive)
#' @param value Associated value (optional)
#' @return The IntervalMap object (invisibly)
#' @export
#' @examples
#' im <- IntervalMap()
#' add(im, 1, 10, "gene1")
#' add(im, 15, 25, "gene2")
add <- function(x, start, end, value = NULL) {
  add_interval(x, as.integer(start), as.integer(end), value)
  invisible(x)
}

#' Build index for efficient queries
#'
#' Must be called after adding intervals and before performing queries.
#'
#' @param x IntervalMap object
#' @return The IntervalMap object (invisibly)
#' @export
#' @examples
#' im <- IntervalMap()
#' add(im, 1, 10, "gene1")
#' build(im)
build <- function(x) {
  build_index(x)
  invisible(x)
}

#' Get interval at specific index
#'
#' @param x IntervalMap object
#' @param idx Index position (1-based)
#' @return List with start, end, and value components
#' @export
#' @examples
#' im <- IntervalMap()
#' add(im, 1, 10, "gene1")
#' build(im)
#' at(im, 1)
at <- function(x, idx) {
  get_interval_at(x, as.integer(idx))
}

#' Get start position at specific index
#'
#' @param x IntervalMap object
#' @param idx Index position (1-based)
#' @return Start position
#' @export
starts_at <- function(x, idx) {
  get_start_at(x, as.integer(idx))
}

#' Get end position at specific index
#'
#' @param x IntervalMap object
#' @param idx Index position (1-based)
#' @return End position
#' @export
ends_at <- function(x, idx) {
  get_end_at(x, as.integer(idx))
}

#' Get data at specific index
#'
#' @param x IntervalMap object
#' @param idx Index position (1-based)
#' @return Associated data value
#' @export
data_at <- function(x, idx) {
  get_data_at(x, as.integer(idx))
}

#' Clear all intervals
#'
#' @param x IntervalMap object
#' @return The IntervalMap object (invisibly)
#' @export
clear <- function(x) {
  clear_intervals(x)
  invisible(x)
}

#' Reserve space for intervals
#'
#' @param x IntervalMap object
#' @param n Number of intervals to reserve space for
#' @return The IntervalMap object (invisibly)
#' @export
reserve <- function(x, n) {
  reserve_intervals(x, as.integer(n))
  invisible(x)
}

#' Get number of intervals
#'
#' @param x IntervalMap object
#' @return Number of intervals
#' @export
size <- function(x) {
  get_size(x)
}

#' Check if any intervals overlap with a range
#'
#' @param x IntervalMap object
#' @param start Start of query range (inclusive)
#' @param end End of query range (inclusive)
#' @return TRUE if any intervals overlap, FALSE otherwise
#' @export
#' @examples
#' im <- IntervalMap()
#' add(im, 1, 10, "gene1")
#' build(im)
#' has_overlaps(im, 5, 15)  # TRUE
has_overlaps <- function(x, start, end) {
  cpp_has_overlaps(x, as.integer(start), as.integer(end))
}

#' Count overlapping intervals
#'
#' @param x IntervalMap object
#' @param start Start of query range (inclusive)
#' @param end End of query range (inclusive)
#' @return Number of overlapping intervals
#' @export
#' @examples
#' im <- IntervalMap()
#' add(im, 1, 10, "gene1")
#' add(im, 5, 15, "gene2")
#' build(im)
#' count(im, 7, 12)  # 2
count <- function(x, start, end) {
  count_overlaps(x, as.integer(start), as.integer(end))
}

#' Count overlaps for multiple ranges
#'
#' @param x IntervalMap object
#' @param starts Vector of query start positions
#' @param ends Vector of query end positions
#' @return Integer vector of overlap counts, one per query
#' @export
#' @examples
#' im <- IntervalMap()
#' add(im, 1, 10, "gene1")
#' add(im, 5, 15, "gene2")
#' add(im, 20, 30, "gene3")
#' build(im)
#' count.batch(im, c(7, 25), c(12, 35))  # c(2, 1)
count.batch <- function(x, starts, ends) {
  count_batch(x, as.integer(starts), as.integer(ends))
}

#' Search for values in overlapping intervals
#'
#' @param x IntervalMap object
#' @param start Start of query range (inclusive)
#' @param end End of query range (inclusive)
#' @return List of values from overlapping intervals
#' @export
#' @examples
#' im <- IntervalMap()
#' add(im, 1, 10, "gene1")
#' add(im, 5, 15, "gene2")
#' build(im)
#' search_values(im, 7, 12)  # list("gene1", "gene2")
search_values <- function(x, start, end) {
  cpp_search_values(x, as.integer(start), as.integer(end))
}

#' Search for indices of overlapping intervals
#'
#' @param x IntervalMap object
#' @param start Start of query range (inclusive)
#' @param end End of query range (inclusive)
#' @return Integer vector of 1-based indices
#' @export
#' @examples
#' im <- IntervalMap()
#' add(im, 1, 10, "gene1")
#' add(im, 5, 15, "gene2")
#' build(im)
#' search_idxs(im, 7, 12)  # c(1, 2)
search_idxs <- function(x, start, end) {
  search_indexes(x, as.integer(start), as.integer(end))
}

#' Batch search for indices
#'
#' Search for indices of overlapping intervals for multiple query ranges simultaneously.
#' This is more efficient than calling search_idxs multiple times.
#'
#' @param x IntervalMap object
#' @param starts Integer vector of query start positions (inclusive)
#' @param ends Integer vector of query end positions (inclusive)
#' @return List of integer vectors, where each element contains the 1-based indices
#'   of intervals that overlap with the corresponding query range
#' @export
#' @examples
#' im <- IntervalMap()
#' add(im, 1, 10, "gene1")
#' add(im, 5, 15, "gene2")
#' add(im, 20, 30, "gene3")
#' build(im)
#'
#' # Search multiple ranges at once
#' result <- search_idxs.batch(im, c(7, 25), c(12, 35))
#' # Returns list(c(1, 2), c(3))
search_idxs.batch <- function(x, starts, ends) {
    search_idxs_batch(x, as.integer(starts), as.integer(ends))
}

#' Search for keys (start, end pairs) of overlapping intervals
#'
#' @param x IntervalMap object
#' @param start Start of query range (inclusive)
#' @param end End of query range (inclusive)
#' @return List with start and end vectors
#' @export
#' @examples
#' im <- IntervalMap()
#' add(im, 1, 10, "gene1")
#' add(im, 5, 15, "gene2")
#' build(im)
#' search_keys(im, 7, 12)  # list(start=c(1,5), end=c(10,15))
search_keys <- function(x, start, end) {
  cpp_search_keys(x, as.integer(start), as.integer(end))
}

#' Search for complete items (start, end, value) of overlapping intervals
#'
#' @param x IntervalMap object
#' @param start Start of query range (inclusive)
#' @param end End of query range (inclusive)
#' @return List with start, end, and value vectors
#' @export
#' @examples
#' im <- IntervalMap()
#' add(im, 1, 10, "gene1")
#' add(im, 5, 15, "gene2")
#' build(im)
#' search_items(im, 7, 12)
search_items <- function(x, start, end) {
  cpp_search_items(x, as.integer(start), as.integer(end))
}

#' Get coverage statistics for a range
#'
#' @param x IntervalMap object
#' @param start Start of query range (inclusive)
#' @param end End of query range (inclusive)
#' @return List with count and total_coverage
#' @export
#' @examples
#' im <- IntervalMap()
#' add(im, 1, 10, "gene1")
#' add(im, 5, 15, "gene2")
#' build(im)
#' coverage(im, 1, 20)
coverage <- function(x, start, end) {
  get_coverage(x, as.integer(start), as.integer(end))
}