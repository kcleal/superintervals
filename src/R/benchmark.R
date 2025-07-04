#!/usr/bin/env Rscript
# benchmark.R - SuperIntervals-r vs IRanges (NCLS algorithm)

suppressPackageStartupMessages({
  library(superintervalsr)
  library(IRanges)
  library(microbenchmark)
})

cat("SuperIntervals vs IRanges (NCLS) Benchmark\n")
cat("===========================================\n")



create_test_intervals <- function(n, max_pos = 1000000, min_width = 100, max_width = 50000) {
  starts <- sort(sample(1:(max_pos - max_width), n))
  widths <- sample(min_width:max_width, n, replace = TRUE)
  ends <- starts + widths - 1
  list(starts = starts, ends = ends)
}

create_superintervals <- function(intervals) {
  si <- IntervalMap()
  reserve(si, length(intervals$starts))
  for (i in seq_along(intervals$starts)) {
    add(si, intervals$starts[i], intervals$ends[i], i)  # Updated from add_int_value
  }
  build(si)  # Updated from index(si)
  return(si)
}

create_iranges <- function(intervals) {
  IRanges(start = intervals$starts, end = intervals$ends)
}

create_test_queries <- function(intervals, n_queries = 10) {
  range_start <- min(intervals$starts)
  range_end <- max(intervals$ends)
  range_width <- range_end - range_start

  queries <- list()

  # Create queries of varying sizes that are guaranteed to hit some intervals
  for (i in 1:n_queries) {
    # Position queries to sample different parts of the data
    q_center <- range_start + (i - 1) * range_width / (n_queries - 1)

    # Vary query sizes: some small, some medium, some large
    if (i <= n_queries/3) {
      # Small queries
      q_width <- range_width * runif(1, 0.001, 0.01)
    } else if (i <= 2*n_queries/3) {
      # Medium queries
      q_width <- range_width * runif(1, 0.01, 0.05)
    } else {
      # Large queries
      q_width <- range_width * runif(1, 0.05, 0.10)
    }

    q_start <- max(range_start, q_center - q_width/2)
    q_end <- min(range_end, q_center + q_width/2)

    queries[[i]] <- list(start = round(q_start), end = round(q_end))
  }

  return(queries)
}

# Improved build time measurement
measure_build_time <- function(build_func, n_runs = 5) {
  times <- numeric(n_runs)

  for (i in 1:n_runs) {
    gc()  # Force garbage collection before timing
    start_time <- Sys.time()
    result <- build_func()
    end_time <- Sys.time()
    times[i] <- as.numeric(end_time - start_time, units = "secs")

    if (i == 1) final_result <- result  # Keep first result
  }
  list(time = median(times), object = final_result)
}

# Validation function to ensure both libraries return identical results
validate_results <- function(si, ir, test_queries = 20, max_pos = 1000000) {
  cat("Validating results consistency...\n")

  all_passed <- TRUE
  failed_tests <- list()

  # Generate random test queries
  set.seed(123)  # For reproducible validation
  for (i in 1:test_queries) {
    # Generate random query interval
    query_start <- sample(1:(max_pos - 1000), 1)
    query_end <- query_start + sample(100:2000, 1)

    # Get results from both libraries
    si_indexes <- search_idxs(si, query_start, query_end)

    query_ir <- IRanges(query_start, query_end)
    ir_hits <- findOverlaps(query_ir, ir)
    ir_indexes <- subjectHits(ir_hits)

    # SuperIntervals already returns 1-based indexing
    si_sorted <- sort(as.integer(si_indexes))
    ir_sorted <- sort(as.integer(ir_indexes))

    # Compare results
    if (!identical(si_sorted, ir_sorted)) {
      all_passed <- FALSE
      failed_tests[[length(failed_tests) + 1]] <- list(
        test_num = i,
        query_start = query_start,
        query_end = query_end,
        si_result = si_sorted,
        ir_result = ir_sorted,
        si_length = length(si_sorted),
        ir_length = length(ir_sorted)
      )

      cat(sprintf("  FAILED Test %d: Query [%d, %d] - SI: %d hits, IR: %d hits\n",
                  i, query_start, query_end, length(si_sorted), length(ir_sorted)))
    }
  }

  if (all_passed) {
    return(TRUE)
  } else {
    cat(sprintf("  ✗ %d out of %d validation tests FAILED\n",
                length(failed_tests), test_queries))

    # Print details of first few failures
    for (i in seq_len(min(3, length(failed_tests)))) {
      failure <- failed_tests[[i]]
      cat(sprintf("    Failure %d: Query [%d, %d]\n",
                  failure$test_num, failure$query_start, failure$query_end))
      cat(sprintf("      SuperIntervals: %s\n",
                  paste(failure$si_result[1:min(10, length(failure$si_result))], collapse = ", ")))
      cat(sprintf("      IRanges:        %s\n",
                  paste(failure$ir_result[1:min(10, length(failure$ir_result))], collapse = ", ")))
    }

    return(FALSE)
  }
}

# Comprehensive validation function
run_validation <- function(si, ir) {
  cat("\nRUNNING VALIDATION TESTS\n")
  cat("========================\n")

  # Basic validation
  overall_passed <- validate_results(si, ir, test_queries = 10)
  if (overall_passed) {
    cat("\n✓ ALL VALIDATION TESTS PASSED\n")
  } else {
    cat("\n✗ VALIDATION FAILED\n")
    stop("Benchmark aborted due to validation failure.")
  }
  return(overall_passed)
}

run_detailed_benchmark <- function(si, ir, intervals) {
  queries <- create_test_queries(intervals, n_queries = 5)

  cat(sprintf("  Testing %d diverse queries...\n", length(queries)))

  # Pre-build all query IRanges objects to ensure fair comparison
  query_iranges <- lapply(queries, function(q) {
    IRanges(q$start, q$end)
  })

  all_si_count_times <- c()
  all_ir_count_times <- c()
  all_si_find_times <- c()
  all_ir_find_times <- c()

  total_si_hits <- 0
  total_ir_hits <- 0

  for (i in seq_along(queries)) {
    query <- queries[[i]]
    query_ir <- query_iranges[[i]]

    # Validate this specific query first
    si_indexes <- search_idxs(si, query$start, query$end)
    ir_hits <- findOverlaps(query_ir, ir)
    ir_indexes <- subjectHits(ir_hits)

    si_count <- length(si_indexes)
    ir_count <- length(ir_indexes)
    total_si_hits <- total_si_hits + si_count
    total_ir_hits <- total_ir_hits + ir_count

    # Ensure they return the same number of results
    if (si_count != ir_count) {
      stop(sprintf("Query %d validation failed: SI=%d hits, IR=%d hits", i, si_count, ir_count))
    }

    # Benchmark counting operations
    count_bench <- microbenchmark(
      SuperIntervals = count(si, query$start, query$end),
      IRanges = length(findOverlaps(query_ir, ir)),
      times = 5,
      unit = "milliseconds"
    )

    # Benchmark finding operations
    find_bench <- microbenchmark(
      SuperIntervals = search_idxs(si, query$start, query$end),  # Changed from search_values to search_idxs
      IRanges = {
        hits <- findOverlaps(query_ir, ir)
        subjectHits(hits)
      },
      times = 5,
      unit = "milliseconds"
    )

    # Extract median times and convert to microseconds for display
    count_summary <- summary(count_bench)
    si_count_time <- count_summary[count_summary$expr == "SuperIntervals", "median"] * 1000
    ir_count_time <- count_summary[count_summary$expr == "IRanges", "median"] * 1000

    find_summary <- summary(find_bench)
    si_find_time <- find_summary[find_summary$expr == "SuperIntervals", "median"] * 1000
    ir_find_time <- find_summary[find_summary$expr == "IRanges", "median"] * 1000

    all_si_count_times <- c(all_si_count_times, si_count_time)
    all_ir_count_times <- c(all_ir_count_times, ir_count_time)
    all_si_find_times <- c(all_si_find_times, si_find_time)
    all_ir_find_times <- c(all_ir_find_times, ir_find_time)

    if (i <= 5) {
      cat(sprintf("    Query %d [%d, %d]: %d hits\n",
                  i, query$start, query$end, si_count))
    }
  }

  return(list(
    si_count_median = median(all_si_count_times),
    ir_count_median = median(all_ir_count_times),
    si_find_median = median(all_si_find_times),
    ir_find_median = median(all_ir_find_times),
    total_queries = length(queries),
    avg_hits = total_si_hits / length(queries)
  ))
}

# Test different dataset sizes
dataset_sizes <- c(1000, 5000, 10000)
results <- list()

for (size in dataset_sizes) {
  cat("Testing with", size, "intervals...\n")

  # Create test data
  set.seed(42)
  intervals <- create_test_intervals(size)

  # Measure build times properly with multiple runs
  cat("  Building SuperIntervals...")
  si_build <- measure_build_time(function() create_superintervals(intervals))
  cat(" done (", round(si_build$time, 4), "s)\n")

  cat("  Building IRanges...")
  ir_build <- measure_build_time(function() create_iranges(intervals))
  cat(" done (", round(ir_build$time, 4), "s)\n")

  validation_passed <- run_validation(si_build$object, ir_build$object)

  if (validation_passed) {
    cat("  Running comprehensive benchmarks...\n")
    bench_results <- run_detailed_benchmark(si_build$object, ir_build$object, intervals)

    # Store results
    results[[paste0("n_", size)]] <- list(
      size = size,
      validation_passed = TRUE,
      build_time_si = si_build$time,
      build_time_ir = ir_build$time,
      count_si_median = bench_results$si_count_median,
      count_ir_median = bench_results$ir_count_median,
      find_si_median = bench_results$si_find_median,
      find_ir_median = bench_results$ir_find_median,
      avg_hits_per_query = bench_results$avg_hits,
      total_queries_tested = bench_results$total_queries
    )
  } else {
    results[[paste0("n_", size)]] <- list(
      size = size,
      validation_passed = FALSE
    )
  }
  cat("  Complete!\n\n")
}

cat("BENCHMARK RESULTS: SuperIntervals vs IRanges (NCLS)\n")
cat("====================================================\n\n")

cat("0. VALIDATION SUMMARY\n")
cat("--------------------\n")
cat(sprintf("%-10s %-15s\n", "Size", "Validation"))
cat(sprintf("%-10s %-15s\n", "----", "----------"))

all_validations_passed <- TRUE
for (name in names(results)) {
  result <- results[[name]]
  validation_status <- if (result$validation_passed) "✓ PASSED" else "✗ FAILED"
  cat(sprintf("%-10s %-15s\n", result$size, validation_status))

  if (!result$validation_passed) {
    all_validations_passed <- FALSE
  }
}

if (all_validations_passed) {
  cat("\n✓ All validations passed - benchmark results are reliable!\n")
} else {
  cat("\n⚠ Some validations failed - benchmark results may not be meaningful!\n")
}
cat("\n")

if (all_validations_passed) {

  cat("1. BUILD TIME COMPARISON\n")
  cat("------------------------\n")
  cat(sprintf("%-10s %-15s %-15s %-10s\n", "Size", "SuperIntervals", "IRanges", "Speedup"))
  cat(sprintf("%-10s %-15s %-15s %-10s\n", "----", "-------------", "--------", "-------"))

  for (name in names(results)) {
    result <- results[[name]]
    if (result$validation_passed) {
      speedup <- result$build_time_ir / result$build_time_si
      cat(sprintf("%-10s %-15.4f %-15.4f %-10.2f\n",
                  result$size, result$build_time_si, result$build_time_ir, speedup))
    }
  }

  cat("\n2. OVERLAP COUNTING PERFORMANCE\n")
  cat("-------------------------------\n")
  cat(sprintf("%-10s %-15s %-15s %-10s %-12s\n", "Size", "SuperIntervals", "IRanges", "Speedup", "Avg Hits"))
  cat(sprintf("%-10s %-15s %-15s %-10s %-12s\n", "----", "-------------", "--------", "-------", "--------"))

  for (name in names(results)) {
    result <- results[[name]]
    if (result$validation_passed) {
      speedup <- result$count_ir_median / result$count_si_median
      cat(sprintf("%-10s %-15.1f %-15.1f %-10.1f %-12.1f\n",
                  result$size, result$count_si_median, result$count_ir_median,
                  speedup, result$avg_hits_per_query))
    }
  }

  cat("\n3. OVERLAP FINDING PERFORMANCE\n")
  cat("------------------------------\n")
  cat(sprintf("%-10s %-15s %-15s %-10s %-12s\n", "Size", "SuperIntervals", "IRanges", "Speedup", "Avg Hits"))
  cat(sprintf("%-10s %-15s %-15s %-10s %-12s\n", "----", "-------------", "--------", "-------", "--------"))

  for (name in names(results)) {
    result <- results[[name]]
    if (result$validation_passed) {
      speedup <- result$find_ir_median / result$find_si_median
      cat(sprintf("%-10s %-15.1f %-15.1f %-10.1f %-12.1f\n",
                  result$size, result$find_si_median, result$find_ir_median,
                  speedup, result$avg_hits_per_query))
    }
  }

  # Calculate average improvements
  cat("\n4. SUMMARY\n")
  cat("----------\n")

  all_count_speedups <- numeric()
  all_find_speedups <- numeric()

  for (name in names(results)) {
    result <- results[[name]]
    if (result$validation_passed) {
      all_count_speedups <- c(all_count_speedups, result$count_ir_median / result$count_si_median)
      all_find_speedups <- c(all_find_speedups, result$find_ir_median / result$find_si_median)
    }
  }

  cat("Performance vs IRanges (NCLS algorithm):\n")
  cat(sprintf("- Overlap counting: %.1fx faster\n", mean(all_count_speedups)))
  cat(sprintf("- Overlap finding:  %.1fx faster\n", mean(all_find_speedups)))

} else {
  cat("Performance results not shown due to validation failures.\n")
}

cat("\nBenchmark completed!\n")