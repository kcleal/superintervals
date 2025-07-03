library(testthat)
library(superintervals)

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