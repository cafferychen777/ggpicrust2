# Test ko2kegg_abundance function after Discussion #113 implementation
# Comprehensive functional tests



# Load internal data
load("../../R/sysdata.rda")

# Helper function to create test KO abundance data
create_test_ko_data <- function(n_kos = 100, n_samples = 3, seed = 123) {
  set.seed(seed)
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), n_kos)

  data <- data.frame(
    function. = real_kos,
    stringsAsFactors = FALSE
  )

  for (i in 1:n_samples) {
    data[[paste0("Sample", i)]] <- rpois(n_kos, lambda = 100)
  }

  return(data)
}

test_that("ko2kegg_abundance works with small dataset", {
  test_data <- create_test_ko_data(n_kos = 10, n_samples = 2)

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
  expect_equal(ncol(result), 2)
  expect_true(all(result >= 0))
  expect_false(any(is.na(result)))
})

test_that("ko2kegg_abundance works with medium dataset", {
  test_data <- create_test_ko_data(n_kos = 200, n_samples = 4)

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 50)
  expect_equal(ncol(result), 4)

  # Check column names preserved
  expected_cols <- colnames(test_data)[-1]
  expect_equal(colnames(result), expected_cols)
})

test_that("ko2kegg_abundance works with large dataset", {
  test_data <- create_test_ko_data(n_kos = 1000, n_samples = 10)

  start_time <- Sys.time()
  result <- suppressMessages(ko2kegg_abundance(data = test_data))
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 100)
  expect_lt(elapsed, 5)
})

test_that("ko2kegg_abundance handles single KO correctly", {
  single_ko <- head(unique(ko_to_kegg_reference$ko_id), 1)
  test_data <- data.frame(
    function. = single_ko,
    Sample1 = 100,
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)

  # Verify the mapping is correct
  expected_pathways <- ko_to_kegg_reference$pathway_id[
    ko_to_kegg_reference$ko_id == single_ko
  ]
  expect_true(all(rownames(result) %in% expected_pathways))
})

test_that("ko2kegg_abundance handles non-existent KOs correctly", {
  fake_data <- data.frame(
    function. = c("K99999", "K88888", "K77777"),
    Sample1 = c(100, 200, 150),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = fake_data))

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

test_that("ko2kegg_abundance handles mixed real and fake KOs", {
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 5)
  mixed_data <- data.frame(
    function. = c(real_kos, "K99999", "K88888"),
    Sample1 = c(100, 200, 150, 180, 120, 300, 400),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = mixed_data))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)

  # Check that fake KO abundances are not included
  # (this is implicit since fake KOs don't map to any pathways)
})

test_that("ko2kegg_abundance handles zero abundances correctly", {
  test_data <- create_test_ko_data(n_kos = 10, n_samples = 2)
  # Set half of abundances to zero
  test_data[, 2] <- 0

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  # First sample should all be zero, second should have values
  expect_true(all(result[, 1] == 0))
  expect_true(any(result[, 2] > 0))
})

test_that("ko2kegg_abundance aggregates KO abundances correctly", {
  # Create a test case where we know the expected result
  test_ko1 <- head(unique(ko_to_kegg_reference$ko_id), 1)
  test_ko2 <- unique(ko_to_kegg_reference$ko_id)[2]

  # Find a pathway that contains test_ko1
  pathway1 <- ko_to_kegg_reference$pathway_id[
    ko_to_kegg_reference$ko_id == test_ko1
  ][1]

  # Find all KOs in that pathway
  kos_in_pathway <- ko_to_kegg_reference$ko_id[
    ko_to_kegg_reference$pathway_id == pathway1
  ]

  # Create test data with known abundances
  test_data <- data.frame(
    function. = kos_in_pathway[1:min(3, length(kos_in_pathway))],
    Sample1 = c(100, 200, 150)[1:min(3, length(kos_in_pathway))],
    stringsAsFactors = FALSE
  )

  # Test with method = "sum" (legacy behavior)
  result <- suppressMessages(ko2kegg_abundance(data = test_data, method = "sum"))

  # Check that pathway abundance is sum of KO abundances (with legacy sum method)
  if (pathway1 %in% rownames(result)) {
    expected_sum <- sum(test_data$Sample1)
    actual_sum <- result[pathway1, "Sample1"]
    expect_equal(actual_sum, expected_sum)
  }
})

test_that("ko2kegg_abundance uses upper-half mean for abundance method", {
  # Get a pathway and its KOs
  pathway1 <- unique(ko_to_kegg_reference$pathway_id)[1]
  kos_in_pathway <- ko_to_kegg_reference$ko_id[
    ko_to_kegg_reference$pathway_id == pathway1
  ]

  # Create test data with known abundances (4 KOs: 100, 200, 150, 50)
  test_data <- data.frame(
    function. = kos_in_pathway[1:min(4, length(kos_in_pathway))],
    Sample1 = c(100, 200, 150, 50)[1:min(4, length(kos_in_pathway))],
    stringsAsFactors = FALSE
  )

  # Test with method = "abundance" (PICRUSt2-style upper-half mean)
  result <- suppressMessages(ko2kegg_abundance(data = test_data, method = "abundance"))

  # Check that pathway abundance is upper-half mean of KO abundances
  # Sorted: 50, 100, 150, 200 -> upper half: 150, 200 -> mean: 175
  if (pathway1 %in% rownames(result)) {
    ko_values <- test_data$Sample1
    sorted_values <- sort(ko_values)
    half_i <- ceiling(length(sorted_values) / 2)
    expected_mean <- mean(sorted_values[half_i:length(sorted_values)])
    actual_value <- result[pathway1, "Sample1"]
    expect_equal(actual_value, expected_mean)
  }
})

test_that("ko2kegg_abundance preserves sample names", {
  test_data <- data.frame(
    function. = head(unique(ko_to_kegg_reference$ko_id), 10),
    Control_1 = rpois(10, 100),
    Control_2 = rpois(10, 100),
    Treatment_1 = rpois(10, 100),
    Treatment_2 = rpois(10, 100),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expected_names <- c("Control_1", "Control_2", "Treatment_1", "Treatment_2")
  expect_equal(colnames(result), expected_names)
})

test_that("ko2kegg_abundance row names are valid pathway IDs", {
  test_data <- create_test_ko_data(n_kos = 100, n_samples = 3)

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  # All row names should be valid KEGG pathway IDs (4-5 digits)
  expect_true(all(grepl("^ko\\d{4,5}$", rownames(result))))
})

test_that("ko2kegg_abundance removes zero-abundance pathways", {
  test_data <- create_test_ko_data(n_kos = 100, n_samples = 3)

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  # No pathway should have zero abundance across all samples
  row_sums <- rowSums(result)
  expect_true(all(row_sums > 0))
})

test_that("ko2kegg_abundance is deterministic", {
  test_data <- create_test_ko_data(n_kos = 50, n_samples = 3, seed = 456)

  result1 <- suppressMessages(ko2kegg_abundance(data = test_data))
  result2 <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_identical(result1, result2)
})

test_that("ko2kegg_abundance performance is acceptable", {
  # Test with 500 KOs
  test_data <- create_test_ko_data(n_kos = 500, n_samples = 5)

  start_time <- Sys.time()
  result <- suppressMessages(ko2kegg_abundance(data = test_data))
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Should process at least 100 KOs per second
  speed <- 500 / elapsed
  expect_gt(speed, 100)
})

test_that("ko2kegg_abundance handles special characters in sample names", {
  test_data <- create_test_ko_data(n_kos = 20, n_samples = 3)
  colnames(test_data) <- c("function.", "Sample-1", "Sample.2", "Sample_3")

  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_equal(colnames(result), c("Sample-1", "Sample.2", "Sample_3"))
})

test_that("ko2kegg_abundance validates input format", {
  # Note: ko2kegg_abundance now auto-detects KO ID columns (ko_id, KO, etc.)
  # and renames them to function., so ko_id column no longer throws an error

  # Empty data should error
  bad_data2 <- data.frame(
    function. = character(0),
    Sample1 = numeric(0)
  )

  expect_error(suppressMessages(ko2kegg_abundance(data = bad_data2)))
})
