# Create test data
test_that("ko2kegg_abundance works with valid data frame input", {
  # Create mock KO abundance data
  mock_ko_data <- data.frame(
    KO = c("K00001", "K00002", "K00003"),
    Sample1 = c(10, 20, 30),
    Sample2 = c(15, 25, 35)
  )

  # Run function
  result <- ko2kegg_abundance(data = mock_ko_data)

  # Tests
  expect_s3_class(result, "data.frame")
  expect_true(all(result >= 0))  # All abundances should be non-negative
  expect_equal(ncol(result), 2)  # Should have same number of samples as input
})

test_that("ko2kegg_abundance handles empty data correctly", {
  # Create empty KO data (but with correct structure)
  empty_ko_data <- data.frame(
    KO = character(0),
    Sample1 = numeric(0),
    Sample2 = numeric(0)
  )

  result <- ko2kegg_abundance(data = empty_ko_data)
  expect_true(nrow(result) == 0)
})

# 首先添加测试辅助函数
get_test_data <- function(type = "valid") {
  switch(type,
    "valid" = data.frame(
      KO = c("K00001", "K00002", "K00003"),
      Sample1 = c(10, 20, 30),
      Sample2 = c(15, 25, 35)
    ),
    "invalid_ko" = data.frame(
      KO = c("K00001", "NOT_KO", "K123"),
      Sample1 = c(1, 2, 3)
    ),
    "negative" = data.frame(
      KO = c("K00001", "K00002"),
      Sample1 = c(1, -2)
    ),
    "non_numeric" = data.frame(
      KO = c("K00001", "K00002"),
      Sample1 = c("1", "text")
    ),
    "na_values" = data.frame(
      KO = c("K00001", "K00002"),
      Sample1 = c(1, NA)
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
  # Create temporary test file
  temp_file <- tempfile(fileext = ".tsv")
  mock_ko_data <- data.frame(
    KO = c("K00001", "K00002", "K00003"),
    Sample1 = c(10, 20, 30),
    Sample2 = c(15, 25, 35)
  )
  write.table(mock_ko_data, temp_file, sep = "\t", row.names = FALSE)

  # Test file input
  result <- ko2kegg_abundance(file = temp_file)

  # Clean up
  unlink(temp_file)

  # Tests
  expect_s3_class(result, "data.frame")
  expect_true(all(result >= 0))
})

test_that("ko2kegg_abundance preserves sample names", {
  mock_ko_data <- data.frame(
    KO = c("K00001", "K00002"),
    SampleA = c(10, 20),
    SampleB = c(15, 25)
  )

  result <- ko2kegg_abundance(data = mock_ko_data)

  expect_equal(colnames(result), c("SampleA", "SampleB"))
})

test_that("ko2kegg_abundance removes zero-abundance pathways", {
  mock_ko_data <- data.frame(
    KO = c("K00001", "K00002"),
    Sample1 = c(0, 0),
    Sample2 = c(0, 0)
  )

  result <- ko2kegg_abundance(data = mock_ko_data)

  expect_true(nrow(result) < nrow(mock_ko_data))
})

test_that("ko2kegg_abundance handles invalid KO IDs", {
  invalid_data <- get_test_data("invalid_ko")
  expect_warning(
    ko2kegg_abundance(data = invalid_data),
    "Warning during data validation: Found .* KO IDs that don't match the expected format"
  )
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
  expect_warning(
    ko2kegg_abundance(data = na_data),
    "Warning during data validation: Missing values found in columns"
  )
})
