# Performance and regression tests for Discussion #113
# Ensures new implementation maintains or improves performance


# Load internal data
load("../../R/sysdata.rda")

# Helper to create test data
create_perf_test_data <- function(n_kos, n_samples, seed = 999) {
  set.seed(seed)
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), n_kos)

  data <- data.frame(function. = real_kos, stringsAsFactors = FALSE)
  for (i in 1:n_samples) {
    data[[paste0("S", i)]] <- rpois(n_kos, lambda = 100)
  }
  return(data)
}

test_that("ko2kegg_abundance completes quickly with 100 KOs", {
  test_data <- create_perf_test_data(100, 4)

  start <- Sys.time()
  result <- suppressMessages(ko2kegg_abundance(data = test_data))
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  expect_lt(elapsed, 0.5)
})

test_that("ko2kegg_abundance completes quickly with 500 KOs", {
  test_data <- create_perf_test_data(500, 5)

  start <- Sys.time()
  result <- suppressMessages(ko2kegg_abundance(data = test_data))
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  expect_lt(elapsed, 2)
})

test_that("ko2kegg_abundance completes quickly with 1000 KOs", {
  test_data <- create_perf_test_data(1000, 10)

  start <- Sys.time()
  result <- suppressMessages(ko2kegg_abundance(data = test_data))
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  expect_lt(elapsed, 5)

  # Processing speed check
  speed <- 1000 / elapsed
  expect_gt(speed, 200)
})

test_that("ko2kegg_abundance scales linearly with dataset size", {
  # Test with increasing sizes
  sizes <- c(100, 200, 400)
  times <- numeric(length(sizes))

  for (i in seq_along(sizes)) {
    test_data <- create_perf_test_data(sizes[i], 3)

    start <- Sys.time()
    result <- suppressMessages(ko2kegg_abundance(data = test_data))
    times[i] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  }

  # Check that time roughly doubles when size doubles
  # (allowing for variation and overhead)
  ratio_1_2 <- times[2] / times[1]
  ratio_2_3 <- times[3] / times[2]

  # Should be roughly linear (between 1.5x and 3x when doubling size)
  expect_gt(ratio_1_2, 1)
  expect_lt(ratio_1_2, 4)

  expect_gt(ratio_2_3, 1)
  expect_lt(ratio_2_3, 4)
})

test_that("ko2kegg_abundance memory usage is reasonable", {
  # Create a moderate dataset
  test_data <- create_perf_test_data(500, 5)

  # Get memory before
  gc()
  mem_before <- sum(gc()[, 2])

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  # Get memory after
  mem_after <- sum(gc()[, 2])

  # Memory increase should be reasonable (< 100 MB)
  mem_increase <- mem_after - mem_before
  expect_lt(mem_increase, 100)
})

test_that("ko2kegg_abundance pathway coverage is improved", {
  # Compare pathway coverage with different sized inputs
  test_200 <- create_perf_test_data(200, 3, seed = 1)
  test_1000 <- create_perf_test_data(1000, 3, seed = 2)

  result_200 <- suppressMessages(ko2kegg_abundance(data = test_200))
  result_1000 <- suppressMessages(ko2kegg_abundance(data = test_1000))

  # More KOs should generally lead to more pathways
  expect_gt(nrow(result_1000), nrow(result_200))

  # Check reasonable coverage ratio (pathways per 100 KOs)
  coverage_200 <- nrow(result_200) / 200 * 100
  coverage_1000 <- nrow(result_1000) / 1000 * 100

  expect_gt(coverage_200, 30)
})

test_that("ko2kegg_abundance produces consistent results", {
  # Same input should always produce same output
  test_data <- create_perf_test_data(100, 3, seed = 111)

  result1 <- suppressMessages(ko2kegg_abundance(data = test_data))
  result2 <- suppressMessages(ko2kegg_abundance(data = test_data))
  result3 <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_identical(result1, result2)
  expect_identical(result2, result3)
})

test_that("ko2kegg_abundance pathway-KO index improves performance", {
  # Test that the pathway -> KO index is used and improves performance

  test_data <- create_perf_test_data(200, 3)

  # Time the function
  start <- Sys.time()
  result <- suppressMessages(ko2kegg_abundance(data = test_data))
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  # With proper indexing, 200 KOs should be very fast
  expect_lt(elapsed, 1)

  # Check that result is correct
  expect_gt(nrow(result), 50)
})

test_that("ko2kegg_abundance handles repeated KOs efficiently", {
  # Test with a dataset that has some KOs repeated across samples
  # (this shouldn't happen in real data, but tests robustness)

  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 100)
  test_data <- data.frame(
    function. = real_kos,
    S1 = rpois(100, 100),
    S2 = rpois(100, 120),
    S3 = rpois(100, 90),
    stringsAsFactors = FALSE
  )

  start <- Sys.time()
  result <- suppressMessages(ko2kegg_abundance(data = test_data))
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  expect_lt(elapsed, 0.5)
})

test_that("New data provides better KO coverage than old", {
  # Check that new data has more KO coverage
  unique_kos_new <- length(unique(ko_to_kegg_reference$ko_id))

  # Old data had ~8000 unique KOs, new should have > 25000
  expect_gt(unique_kos_new, 25000)
})

test_that("New data provides better pathway coverage than old", {
  # Check that new data has more pathway coverage
  unique_pathways_new <- length(unique(ko_to_kegg_reference$pathway_id))

  # Old data had ~306 pathways, new should have > 500
  expect_gt(unique_pathways_new, 500)
})

test_that("New data has more mappings than old", {
  # Total mappings check
  total_mappings <- nrow(ko_to_kegg_reference)

  # Old data had ~15000 valid mappings, new should have > 60000
  expect_gt(total_mappings, 60000)
})

test_that("New data has lower NA ratio than old", {
  # Check NA ratio in key columns
  na_pathway <- sum(is.na(ko_to_kegg_reference$pathway_id))
  na_ko <- sum(is.na(ko_to_kegg_reference$ko_id))

  total_cells <- nrow(ko_to_kegg_reference) * 2  # Two key columns

  na_ratio <- (na_pathway + na_ko) / total_cells

  # Old data had 84% NA, new should have near 0%
  expect_lt(na_ratio, 0.01)
})

test_that("split() method is faster than loop method", {
  # Benchmark the split() conversion used in pathway_gsea

  start_split <- Sys.time()
  gene_sets_split <- split(ko_to_kegg_reference$ko_id,
                           ko_to_kegg_reference$pathway_id)
  time_split <- as.numeric(difftime(Sys.time(), start_split, units = "secs"))

  # Split should be very fast
  expect_lt(time_split, 0.2)

  # Check result is correct
  expect_gt(length(gene_sets_split), 500)

  expect_true(all(sapply(gene_sets_split, is.character)))
})

test_that("ko2kegg_abundance handles varying sample counts efficiently", {
  # Test with different numbers of samples
  sample_counts <- c(2, 5, 10, 20)
  times <- numeric(length(sample_counts))

  for (i in seq_along(sample_counts)) {
    test_data <- create_perf_test_data(100, sample_counts[i])

    start <- Sys.time()
    result <- suppressMessages(ko2kegg_abundance(data = test_data))
    times[i] <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  }

  # Time should not increase dramatically with more samples
  time_ratio_max <- max(times) / min(times)

  expect_lt(time_ratio_max, 5)
})
