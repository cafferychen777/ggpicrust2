# Complete tests for Discussion #113 implementation
# Covers data structure, functionality, performance, and edge cases

# Load package - sysdata.rda will be loaded automatically by the package
library(ggpicrust2)
devtools::load_all("../..", quiet = TRUE)

# Access internal data from package namespace for structure tests
ko_to_kegg_reference <- getFromNamespace("ko_to_kegg_reference", "ggpicrust2")
ko_pathway_index <- getFromNamespace("ko_pathway_index", "ggpicrust2")

# ==============================================================================
# DATA STRUCTURE TESTS
# ==============================================================================

test_that("ko_to_kegg_reference has correct long-format structure", {
  expect_s3_class(ko_to_kegg_reference, "data.frame")
  expect_gt(nrow(ko_to_kegg_reference), 10000)
  expect_lt(ncol(ko_to_kegg_reference), 20)

  # Check required columns
  required_cols <- c("pathway_id", "ko_id", "pathway_name")
  expect_true(all(required_cols %in% names(ko_to_kegg_reference)))
})

test_that("ko_to_kegg_reference has no NA in key columns", {
  expect_equal(sum(is.na(ko_to_kegg_reference$pathway_id)), 0)
  expect_equal(sum(is.na(ko_to_kegg_reference$ko_id)), 0)
})

test_that("ko_to_kegg_reference IDs have correct format", {
  expect_true(all(grepl("^ko\\d", ko_to_kegg_reference$pathway_id)))
  expect_true(all(grepl("^K\\d", trimws(ko_to_kegg_reference$ko_id))))
})

test_that("ko_to_kegg_reference has sufficient coverage", {
  unique_pathways <- length(unique(ko_to_kegg_reference$pathway_id))
  unique_kos <- length(unique(ko_to_kegg_reference$ko_id))

  expect_gte(unique_pathways, 500)
  expect_gte(unique_kos, 25000)
  expect_gte(nrow(ko_to_kegg_reference), 60000)
})

test_that("ko_pathway_index is correctly structured", {
  expect_type(ko_pathway_index, "list")
  expect_gt(length(ko_pathway_index), 25000)

  # Verify a random mapping
  test_ko <- names(ko_pathway_index)[1]
  expected <- ko_to_kegg_reference$pathway_id[ko_to_kegg_reference$ko_id == test_ko]
  actual <- ko_pathway_index[[test_ko]]
  expect_setequal(actual, expected)
})

# ==============================================================================
# FUNCTIONALITY TESTS
# ==============================================================================

create_test_data <- function(n_kos, n_samples, seed = 123) {
  set.seed(seed)
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), n_kos)

  data <- data.frame(function. = real_kos, stringsAsFactors = FALSE)
  for (i in 1:n_samples) {
    data[[paste0("S", i)]] <- rpois(n_kos, lambda = 100)
  }
  return(data)
}

test_that("ko2kegg_abundance works with small dataset", {
  test_data <- create_test_data(10, 2)
  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
  expect_equal(ncol(result), 2)
  expect_true(all(result >= 0))
  expect_false(any(is.na(result)))
})

test_that("ko2kegg_abundance works with medium dataset", {
  test_data <- create_test_data(200, 4)
  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 50)
  expect_equal(ncol(result), 4)
  expect_equal(colnames(result), colnames(test_data)[-1])
})

test_that("ko2kegg_abundance works with large dataset", {
  test_data <- create_test_data(1000, 10)

  start_time <- Sys.time()
  result <- suppressMessages(ko2kegg_abundance(data = test_data))
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 100)
  expect_lt(elapsed, 5)
})

test_that("ko2kegg_abundance handles single KO", {
  single_ko <- head(unique(ko_to_kegg_reference$ko_id), 1)
  test_data <- data.frame(function. = single_ko, S1 = 100, stringsAsFactors = FALSE)

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)

  # Verify result pathways are a subset of expected pathways
  expected_pathways <- unique(ko_to_kegg_reference$pathway_id[ko_to_kegg_reference$ko_id == single_ko])
  # All results should be in expected pathways (allowing for some to be filtered out due to zero abundance)
  expect_true(all(rownames(result) %in% expected_pathways))
})

test_that("ko2kegg_abundance handles non-existent KOs", {
  fake_data <- data.frame(function. = c("K99999", "K88888"), S1 = c(100, 200), stringsAsFactors = FALSE)
  result <- suppressMessages(ko2kegg_abundance(data = fake_data))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

test_that("ko2kegg_abundance preserves sample names", {
  test_data <- data.frame(
    function. = head(unique(ko_to_kegg_reference$ko_id), 10),
    Ctrl_1 = rpois(10, 100),
    Treat_1 = rpois(10, 100),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))
  expect_equal(colnames(result), c("Ctrl_1", "Treat_1"))
})

test_that("ko2kegg_abundance row names are valid pathway IDs", {
  test_data <- create_test_data(100, 3)
  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  # Pathway IDs can be 4-5 digits (e.g., ko00010 or ko9980)
  expect_true(all(grepl("^ko[0-9]{4,5}$", rownames(result))))
})

test_that("ko2kegg_abundance removes zero-abundance pathways", {
  test_data <- create_test_data(100, 3)
  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  row_sums <- rowSums(result)
  expect_true(all(row_sums > 0))
})

test_that("ko2kegg_abundance is deterministic", {
  test_data <- create_test_data(50, 3, seed = 456)

  result1 <- suppressMessages(ko2kegg_abundance(data = test_data))
  result2 <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_identical(result1, result2)
})

# ==============================================================================
# PERFORMANCE TESTS
# ==============================================================================

test_that("ko2kegg_abundance meets performance targets", {
  sizes <- c(100, 500)
  targets <- c(0.5, 2)  # seconds

  for (i in seq_along(sizes)) {
    test_data <- create_test_data(sizes[i], 3, seed = i)

    start <- Sys.time()
    result <- suppressMessages(ko2kegg_abundance(data = test_data))
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    expect_lt(elapsed, targets[i])
  }
})

test_that("ko2kegg_abundance processing speed is acceptable", {
  test_data <- create_test_data(500, 5)

  start <- Sys.time()
  result <- suppressMessages(ko2kegg_abundance(data = test_data))
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  speed <- 500 / elapsed
  expect_gt(speed, 100)
})

test_that("New data has better coverage than old", {
  expect_gt(length(unique(ko_to_kegg_reference$ko_id)), 25000)
  expect_gt(length(unique(ko_to_kegg_reference$pathway_id)), 500)
  expect_gt(nrow(ko_to_kegg_reference), 60000)
})

test_that("New data has low NA ratio", {
  na_pathway <- sum(is.na(ko_to_kegg_reference$pathway_id))
  na_ko <- sum(is.na(ko_to_kegg_reference$ko_id))
  total_cells <- nrow(ko_to_kegg_reference) * 2
  na_ratio <- (na_pathway + na_ko) / total_cells

  expect_lt(na_ratio, 0.01)
})

# ==============================================================================
# PATHWAY_GSEA INTEGRATION TESTS
# ==============================================================================

test_that("pathway_gsea gene_sets created correctly from long-format", {
  gene_sets <- split(ko_to_kegg_reference$ko_id, ko_to_kegg_reference$pathway_id)

  expect_type(gene_sets, "list")
  expect_gt(length(gene_sets), 500)

  # Check a random pathway
  test_pathway <- names(gene_sets)[100]
  expected_kos <- sort(ko_to_kegg_reference$ko_id[ko_to_kegg_reference$pathway_id == test_pathway])
  actual_kos <- sort(gene_sets[[test_pathway]])

  expect_equal(actual_kos, expected_kos)
})

test_that("pathway_gsea gene_sets size distribution is reasonable", {
  gene_sets <- split(ko_to_kegg_reference$ko_id, ko_to_kegg_reference$pathway_id)
  set_sizes <- sapply(gene_sets, length)

  expect_gt(min(set_sizes), 0)
  expect_lt(max(set_sizes), 5000)
  expect_gt(mean(set_sizes), 10)
  expect_lt(mean(set_sizes), 500)
})

test_that("pathway_gsea conversion is fast", {
  start <- Sys.time()
  gene_sets <- split(ko_to_kegg_reference$ko_id, ko_to_kegg_reference$pathway_id)
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  expect_lt(elapsed, 0.2)
  expect_gt(length(gene_sets), 500)
})

# ==============================================================================
# EDGE CASES TESTS
# ==============================================================================

test_that("ko2kegg_abundance handles all-zero abundances", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 20)
  test_data <- data.frame(function. = real_kos, S1 = rep(0, 20), S2 = rep(0, 20), stringsAsFactors = FALSE)

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
  expect_equal(ncol(result), 2)
})

test_that("ko2kegg_abundance handles single sample", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 20)
  test_data <- data.frame(function. = real_kos, OnlySample = rpois(20, 100), stringsAsFactors = FALSE)

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "OnlySample")
})

test_that("ko2kegg_abundance handles many samples", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 50)
  test_data <- data.frame(function. = real_kos, stringsAsFactors = FALSE)

  for (i in 1:50) {
    test_data[[paste0("S", i)]] <- rpois(50, 100)
  }

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 50)
})

test_that("ko2kegg_abundance handles fractional abundances", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 10)
  test_data <- data.frame(
    function. = real_kos,
    S1 = runif(10, 0.1, 100.5),
    S2 = runif(10, 0.1, 100.5),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
})

test_that("ko2kegg_abundance validates input format", {
  # Missing function column
  bad_data <- data.frame(ko_id = c("K00001", "K00002"), S1 = c(100, 200))
  expect_error(ko2kegg_abundance(data = bad_data))

  # Empty data
  empty_data <- data.frame(function. = character(0), S1 = numeric(0))
  expect_error(ko2kegg_abundance(data = empty_data))
})

test_that("ko2kegg_abundance handles special characters in sample names", {
  test_data <- create_test_data(20, 3)
  colnames(test_data) <- c("function.", "Sample-1", "Sample.2", "Sample_3")

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_equal(colnames(result), c("Sample-1", "Sample.2", "Sample_3"))
})
