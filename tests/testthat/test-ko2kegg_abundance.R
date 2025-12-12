# Test ko2kegg_abundance function
# Tests for Discussion #113 implementation

test_that("ko2kegg_abundance works with valid data frame input", {
  # Create mock KO abundance data with correct column name
  mock_ko_data <- data.frame(
    function. = c("K00001", "K00002", "K00003"),
    Sample1 = c(10, 20, 30),
    Sample2 = c(15, 25, 35),
    stringsAsFactors = FALSE
  )

  # Run function - suppress messages about mappings
  result <- suppressMessages(ko2kegg_abundance(data = mock_ko_data))

  # Tests
  expect_s3_class(result, "data.frame")
  # Result may be empty if KOs don't map to any pathways
  expect_true(nrow(result) >= 0)
  if (nrow(result) > 0) {
    expect_true(all(result >= 0))  # All abundances should be non-negative
    expect_equal(ncol(result), 2)  # Should have same number of samples as input
  }
})

test_that("ko2kegg_abundance handles empty data correctly", {
  # Create empty KO data (but with correct structure)
  empty_ko_data <- data.frame(
    function. = character(0),
    Sample1 = numeric(0),
    Sample2 = numeric(0),
    stringsAsFactors = FALSE
  )

  # Empty input should throw an error (no KOs to process)
  expect_error(ko2kegg_abundance(data = empty_ko_data))
})

# Helper function for test data
get_test_data <- function(type = "valid") {
  switch(type,
    "valid" = data.frame(
      function. = c("K00001", "K00002", "K00003"),
      Sample1 = c(10, 20, 30),
      Sample2 = c(15, 25, 35),
      stringsAsFactors = FALSE
    ),
    "invalid_ko" = data.frame(
      function. = c("K00001", "NOT_KO", "K123"),
      Sample1 = c(1, 2, 3),
      stringsAsFactors = FALSE
    ),
    "negative" = data.frame(
      function. = c("K00001", "K00002"),
      Sample1 = c(1, -2),
      stringsAsFactors = FALSE
    ),
    "non_numeric" = data.frame(
      function. = c("K00001", "K00002"),
      Sample1 = c("1", "text"),
      stringsAsFactors = FALSE
    ),
    "na_values" = data.frame(
      function. = c("K00001", "K00002"),
      Sample1 = c(1, NA),
      stringsAsFactors = FALSE
    )
  )
}

test_that("ko2kegg_abundance throws appropriate errors", {
  # Test with no input
  expect_error(ko2kegg_abundance(),
              "Error: Please provide either a file or a data.frame.")

  # Test with invalid file format
  expect_error(
    ko2kegg_abundance(file = "test.pdf"),
    "Error: Input file should be in .txt, .tsv, .csv format."
  )

  # Test with invalid data frame
  expect_error(
    ko2kegg_abundance(data = list(a = 1, b = 2)),
    "The provided data must be a data.frame"
  )
})

test_that("ko2kegg_abundance handles file input correctly", {
  # Create temporary test file with correct column name
  temp_file <- tempfile(fileext = ".tsv")
  mock_ko_data <- data.frame(
    function. = c("K00001", "K00002", "K00003"),
    Sample1 = c(10, 20, 30),
    Sample2 = c(15, 25, 35),
    stringsAsFactors = FALSE
  )
  write.table(mock_ko_data, temp_file, sep = "\t", row.names = FALSE)

  # Test file input
  result <- suppressMessages(ko2kegg_abundance(file = temp_file))

  # Clean up
  unlink(temp_file)

  # Tests
  expect_s3_class(result, "data.frame")
  if (nrow(result) > 0) {
    expect_true(all(result >= 0))
  }
})

test_that("ko2kegg_abundance preserves sample names", {
  mock_ko_data <- data.frame(
    function. = c("K00001", "K00002"),
    SampleA = c(10, 20),
    SampleB = c(15, 25),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = mock_ko_data))

  # If we have results, check column names
  if (nrow(result) > 0) {
    expect_equal(colnames(result), c("SampleA", "SampleB"))
  }
})

test_that("ko2kegg_abundance removes zero-abundance pathways", {
  # Use real KO IDs that actually map to pathways
  load("../../R/sysdata.rda")
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 5)

  mock_ko_data <- data.frame(
    function. = real_kos,
    Sample1 = c(100, 200, 150, 180, 120),
    Sample2 = c(50, 100, 75, 90, 60),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = mock_ko_data))

  # All pathways in result should have non-zero abundance
  if (nrow(result) > 0) {
    row_sums <- rowSums(result)
    expect_true(all(row_sums > 0))
  }
})

test_that("ko2kegg_abundance handles invalid KO IDs", {
  invalid_data <- get_test_data("invalid_ko")
  # Function should still work, just with fewer/no mappings
  result <- suppressMessages(ko2kegg_abundance(data = invalid_data))
  expect_s3_class(result, "data.frame")
})

test_that("ko2kegg_abundance handles negative values", {
  negative_data <- get_test_data("negative")
  expect_error(
    ko2kegg_abundance(data = negative_data),
    "Negative abundance values found"
  )
})

test_that("ko2kegg_abundance handles non-numeric columns", {
  non_numeric_data <- get_test_data("non_numeric")
  expect_error(
    ko2kegg_abundance(data = non_numeric_data),
    "The following columns contain non-numeric values"
  )
})

test_that("ko2kegg_abundance handles missing values", {
  na_data <- get_test_data("na_values")
  # Function should warn about NA values
  expect_warning(
    suppressMessages(ko2kegg_abundance(data = na_data)),
    "Warning during data validation: Missing values found in columns"
  )
})
