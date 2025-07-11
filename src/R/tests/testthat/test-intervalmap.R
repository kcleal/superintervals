library(testthat)
library(superintervalsr)

test_that("IntervalMap creation and basic operations work", {
  imap <- IntervalMap()
  expect_s3_class(imap, "IntervalMap")
  expect_equal(length(imap), 0)
  expect_equal(size(imap), 0)
})

test_that("Adding intervals with values works", {
  imap <- IntervalMap()
  add(imap, 1, 10, "test_value")
  add(imap, 5, 15, list(name = "gene", score = 100))

  expect_equal(length(imap), 2)

  build(imap)  # Updated from index(imap)

  # Test accessing intervals
  interval1 <- at(imap, 1)
  expect_equal(interval1$start, 1)
  expect_equal(interval1$end, 10)
  expect_equal(interval1$value, "test_value")

  interval2 <- imap[2]
  expect_equal(interval2$start, 5)
  expect_equal(interval2$end, 15)
  expect_equal(interval2$value$name, "gene")
  expect_equal(interval2$value$score, 100)
})

test_that("Adding intervals without values works", {
  imap <- IntervalMap()
  add(imap, 1, 10)  # No value (NULL)
  add(imap, 5, 15, "gene")  # With value
  add(imap, 20, 30)  # No value (NULL)

  expect_equal(length(imap), 3)

  build(imap)

  # Test accessing intervals with mixed data
  interval1 <- at(imap, 1)
  expect_equal(interval1$start, 1)
  expect_equal(interval1$end, 10)
  expect_null(interval1$value)

  interval2 <- at(imap, 2)
  expect_equal(interval2$start, 5)
  expect_equal(interval2$end, 15)
  expect_equal(interval2$value, "gene")

  interval3 <- at(imap, 3)
  expect_equal(interval3$start, 20)
  expect_equal(interval3$end, 30)
  expect_null(interval3$value)
})

test_that("Individual accessor functions work", {
  imap <- IntervalMap()
  add(imap, 1, 10, "test")
  add(imap, 5, 15)  # NULL value
  build(imap)

  # Test starts, ends, data functions
  expect_equal(starts_at(imap, 1), 1)
  expect_equal(ends_at(imap, 1), 10)
  expect_equal(data_at(imap, 1), "test")

  expect_equal(starts_at(imap, 2), 5)
  expect_equal(ends_at(imap, 2), 15)
  expect_null(data_at(imap, 2))
})

test_that("Overlap queries work correctly", {
  imap <- IntervalMap()
  add(imap, 1, 10, "A")
  add(imap, 5, 15, "B")
  add(imap, 20, 30, "C")
  build(imap)  # Updated from index(imap)

  # Test has_overlaps
  expect_true(has_overlaps(imap, 1, 5))
  expect_true(has_overlaps(imap, 8, 12))
  expect_false(has_overlaps(imap, 16, 19))
  expect_true(has_overlaps(imap, 25, 35))

  # Test count
  expect_equal(count(imap, 1, 5), 2)  # Updated from count_overlaps
  expect_equal(count(imap, 8, 12), 2)
  expect_equal(count(imap, 16, 19), 0)
  expect_equal(count(imap, 25, 35), 1)

  # Test search_values
  values1 <- search_values(imap, 1, 5)
  expect_equal(length(values1), 2)
  expect_true("A" %in% values1)
  expect_true("B" %in% values1)

  values2 <- search_values(imap, 25, 35)
  expect_equal(length(values2), 1)
  expect_equal(values2[[1]], "C")

  # Test search_idxs
  idxs1 <- search_idxs(imap, 1, 5)  # Updated from find_indexes
  expect_equal(length(idxs1), 2)
  expect_true(1 %in% idxs1)
  expect_true(2 %in% idxs1)

  # Test search_keys
  keys1 <- search_keys(imap, 1, 5)
  expect_equal(length(keys1$start), 2)
  expect_equal(length(keys1$end), 2)
  expect_true(1 %in% keys1$start)
  expect_true(5 %in% keys1$start)

  # Test search_items
  items1 <- search_items(imap, 1, 5)
  expect_equal(length(items1$start), 2)
  expect_equal(length(items1$end), 2)
  expect_equal(length(items1$value), 2)
  expect_true("A" %in% items1$value)
  expect_true("B" %in% items1$value)
})

test_that("Coverage function works", {
  imap <- IntervalMap()
  add(imap, 1, 10, "A")
  add(imap, 5, 15, "B")
  add(imap, 20, 30, "C")
  build(imap)

  # Test coverage
  cov1 <- coverage(imap, 1, 15)
  expect_equal(cov1$count, 2)
  expect_gt(cov1$total_coverage, 0)

  cov2 <- coverage(imap, 16, 19)
  expect_equal(cov2$count, 0)
  expect_equal(cov2$total_coverage, 0)
})

test_that("Clear and reserve work", {
  imap <- IntervalMap()
  add(imap, 1, 10, "A")
  add(imap, 5, 15, "B")
  expect_equal(length(imap), 2)

  clear(imap)
  expect_equal(length(imap), 0)
  expect_equal(size(imap), 0)

  # Test reserve
  reserve(imap, 1000)
  # Should not change size, just pre-allocate
  expect_equal(length(imap), 0)
})

test_that("Error handling for out of range access", {
  imap <- IntervalMap()
  add(imap, 1, 10, "A")
  build(imap)

  expect_error(at(imap, 0), "Index out of range")
  expect_error(at(imap, 2), "Index out of range")
  expect_error(starts_at(imap, 0), "Index out of range")
  expect_error(ends_at(imap, 2), "Index out of range")
  expect_error(data_at(imap, -1), "Index out of range")
})

test_that("Complex data types can be stored", {
  imap <- IntervalMap()

  # Store a data frame
  df <- data.frame(id = 1:3, name = c("A", "B", "C"))
  add(imap, 1, 10, df)

  # Store a function
  my_func <- function(x) x^2
  add(imap, 15, 25, my_func)

  # Store a list
  my_list <- list(genes = c("TP53", "BRCA1"), scores = c(0.95, 0.87))
  add(imap, 30, 40, my_list)

  build(imap)

  # Test retrieval
  interval1 <- at(imap, 1)
  expect_s3_class(interval1$value, "data.frame")
  expect_equal(nrow(interval1$value), 3)

  interval2 <- at(imap, 2)
  expect_true(is.function(interval2$value))
  expect_equal(interval2$value(3), 9)

  interval3 <- at(imap, 3)
  expect_true(is.list(interval3$value))
  expect_equal(interval3$value$genes[1], "TP53")
})

test_that("Print method works", {
  imap <- IntervalMap()
  add(imap, 1, 10, "A")
  add(imap, 5, 15, "B")

  # Test that print doesn't error
  expect_output(print(imap), "IntervalMap with 2 intervals")

  # Test empty IntervalMap print
  empty_imap <- IntervalMap()
  expect_output(print(empty_imap), "IntervalMap with 0 intervals")
})


test_that("IntervalMap.from_vectors works with all parameters", {
  # Test with starts, ends, and values
  starts <- c(1, 10, 20)
  ends <- c(5, 15, 25)
  values <- c("A", "B", "C")

  imap <- IntervalMap.from_vectors(starts, ends, values)
  expect_s3_class(imap, "IntervalMap")
  expect_equal(length(imap), 3)

  # Test accessing intervals
  interval1 <- at(imap, 1)
  expect_equal(interval1$start, 1)
  expect_equal(interval1$end, 5)
  expect_equal(interval1$value, "A")

  interval2 <- at(imap, 2)
  expect_equal(interval2$start, 10)
  expect_equal(interval2$end, 15)
  expect_equal(interval2$value, "B")

  interval3 <- at(imap, 3)
  expect_equal(interval3$start, 20)
  expect_equal(interval3$end, 25)
  expect_equal(interval3$value, "C")
})

test_that("IntervalMap.from_vectors works without values", {
  # Test with just starts and ends (no values)
  starts <- c(1, 10, 20)
  ends <- c(5, 15, 25)

  imap <- IntervalMap.from_vectors(starts, ends)
  expect_s3_class(imap, "IntervalMap")
  expect_equal(length(imap), 3)

  # Test that values are NULL
  interval1 <- at(imap, 1)
  expect_equal(interval1$start, 1)
  expect_equal(interval1$end, 5)
  expect_null(interval1$value)

  interval2 <- at(imap, 2)
  expect_equal(interval2$start, 10)
  expect_equal(interval2$end, 15)
  expect_null(interval2$value)
})

test_that("IntervalMap.from_vectors works with complex values", {
  # Test with complex data types as values
  starts <- c(1, 10)
  ends <- c(5, 15)
  values <- list(
    list(name = "gene1", score = 0.95),
    data.frame(id = 1:2, name = c("A", "B"))
  )

  imap <- IntervalMap.from_vectors(starts, ends, values)
  expect_equal(length(imap), 2)

  interval1 <- at(imap, 1)
  expect_equal(interval1$value$name, "gene1")
  expect_equal(interval1$value$score, 0.95)

  interval2 <- at(imap, 2)
  expect_s3_class(interval2$value, "data.frame")
  expect_equal(nrow(interval2$value), 2)
})

test_that("IntervalMap.from_vectors handles empty inputs", {
  # Test with empty vectors
  imap <- IntervalMap.from_vectors(integer(0), integer(0))
  expect_s3_class(imap, "IntervalMap")
  expect_equal(length(imap), 0)

  # Test with empty vectors but with values parameter
  imap2 <- IntervalMap.from_vectors(integer(0), integer(0), list())
  expect_s3_class(imap2, "IntervalMap")
  expect_equal(length(imap2), 0)
})

test_that("count.batch works correctly", {
  # Create test IntervalMap
  imap <- IntervalMap()
  add(imap, 1, 10, "A")
  add(imap, 5, 15, "B")
  add(imap, 20, 30, "C")
  add(imap, 25, 35, "D")
  build(imap)

  # Test batch counting
  query_starts <- c(1, 8, 16, 27)
  query_ends <- c(5, 12, 19, 32)

  counts <- count.batch(imap, query_starts, query_ends)
  expect_equal(length(counts), 4)
  expect_equal(counts[1], 2)  # overlaps with A and B
  expect_equal(counts[2], 2)  # overlaps with A and B
  expect_equal(counts[3], 0)  # no overlaps
  expect_equal(counts[4], 2)  # overlaps with C and D
})

test_that("count.batch handles edge cases", {
  imap <- IntervalMap()
  add(imap, 1, 10, "A")
  build(imap)

  # Test with empty query vectors
  counts <- count.batch(imap, integer(0), integer(0))
  expect_equal(length(counts), 0)

  # Test with single query
  counts <- count.batch(imap, 5, 8)
  expect_equal(length(counts), 1)
  expect_equal(counts[1], 1)

  # Test with queries that don't overlap
  counts <- count.batch(imap, c(15, 20), c(18, 25))
  expect_equal(length(counts), 2)
  expect_equal(counts[1], 0)
  expect_equal(counts[2], 0)
})

test_that("search_idxs.batch works correctly", {
  # Create test IntervalMap
  imap <- IntervalMap()
  add(imap, 1, 10, "A")
  add(imap, 5, 15, "B")
  add(imap, 20, 30, "C")
  add(imap, 25, 35, "D")
  build(imap)

  # Test batch index searching
  query_starts <- c(1, 8, 16, 27)
  query_ends <- c(5, 12, 19, 32)

  result <- search_idxs.batch(imap, query_starts, query_ends)
  expect_true(is.list(result))
  expect_equal(length(result), 4)

  # Check first query (should find intervals 1 and 2)
  expect_equal(length(result[[1]]), 2)
  expect_true(1 %in% result[[1]])
  expect_true(2 %in% result[[1]])

  # Check second query (should find intervals 1 and 2)
  expect_equal(length(result[[2]]), 2)
  expect_true(1 %in% result[[2]])
  expect_true(2 %in% result[[2]])

  # Check third query (should find no intervals)
  expect_equal(length(result[[3]]), 0)

  # Check fourth query (should find intervals 3 and 4)
  expect_equal(length(result[[4]]), 2)
  expect_true(3 %in% result[[4]])
  expect_true(4 %in% result[[4]])
})

test_that("search_idxs.batch handles edge cases", {
  imap <- IntervalMap()
  add(imap, 1, 10, "A")
  build(imap)

  # Test with empty query vectors
  result <- search_idxs.batch(imap, integer(0), integer(0))
  expect_true(is.list(result))
  expect_equal(length(result), 0)

  # Test with single query
  result <- search_idxs.batch(imap, 5, 8)
  expect_true(is.list(result))
  expect_equal(length(result), 1)
  expect_equal(length(result[[1]]), 1)
  expect_equal(result[[1]][1], 1)
})

test_that("Batch functions work with large datasets", {
  # Create a larger test dataset
  n_intervals <- 100
  starts <- seq(1, 1000, length.out = n_intervals)
  ends <- starts + 50
  values <- paste0("interval_", 1:n_intervals)

  imap <- IntervalMap.from_vectors(as.integer(starts), as.integer(ends), values)
  expect_equal(length(imap), n_intervals)

  # Test batch operations on larger dataset
  n_queries <- 20
  query_starts <- seq(10, 900, length.out = n_queries)
  query_ends <- query_starts + 30

  # Test count.batch
  counts <- count.batch(imap, as.integer(query_starts), as.integer(query_ends))
  expect_equal(length(counts), n_queries)
  expect_true(all(counts >= 0))

  # Test search_idxs.batch
  result <- search_idxs.batch(imap, as.integer(query_starts), as.integer(query_ends))
  expect_equal(length(result), n_queries)
  expect_true(is.list(result))

  # Verify consistency between count.batch and search_idxs.batch
  for (i in 1:n_queries) {
    expect_equal(counts[i], length(result[[i]]))
  }
})

test_that("Integration test: from_vectors with batch operations", {
  # Create IntervalMap using from_vectors
  starts <- c(1, 10, 20, 30, 40)
  ends <- c(8, 18, 28, 38, 48)
  values <- c("A", "B", "C", "D", "E")

  imap <- IntervalMap.from_vectors(starts, ends, values)

  # Test that it works with batch operations without needing build()
  query_starts <- c(5, 15, 25, 35)
  query_ends <- c(12, 22, 32, 42)

  counts <- count.batch(imap, query_starts, query_ends)
  expect_equal(length(counts), 4)
  expect_true(all(counts > 0))  # All queries should find overlaps

  idxs <- search_idxs.batch(imap, query_starts, query_ends)
  expect_equal(length(idxs), 4)

  # Verify consistency
  for (i in 1:4) {
    expect_equal(counts[i], length(idxs[[i]]))
  }
})

test_that("Error handling for batch functions", {
  imap <- IntervalMap()
  add(imap, 1, 10, "A")
  build(imap)

  # Test mismatched vector lengths
  expect_error(count.batch(imap, c(1, 5), c(3, 8, 12)), "length")
  expect_error(search_idxs.batch(imap, c(1, 5), c(3, 8, 12)), "length")
})

test_that("Type conversion works correctly", {
  # Test that the functions handle different numeric types
  imap <- IntervalMap()
  add(imap, 1.0, 10.0, "A")  # doubles should be converted to integers
  add(imap, 5L, 15L, "B")    # integers should work
  build(imap)

  # Test with double inputs to batch functions
  counts <- count.batch(imap, c(1.0, 8.0), c(5.0, 12.0))
  expect_equal(length(counts), 2)
  expect_equal(counts[1], 2)
  expect_equal(counts[2], 2)

  # Test from_vectors with doubles
  imap2 <- IntervalMap.from_vectors(c(1.0, 10.0), c(5.0, 15.0), c("X", "Y"))
  expect_equal(length(imap2), 2)
  interval1 <- at(imap2, 1)
  expect_equal(interval1$start, 1)
  expect_equal(interval1$end, 5)
})