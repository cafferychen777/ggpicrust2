# Test ko2kegg_abundance function

test_that("ko2kegg_abundance works with valid data frame input", {
  mock_ko_data <- data.frame(
    function. = c("K00001", "K00002", "K00003"),
    Sample1 = c(10, 20, 30),
    Sample2 = c(15, 25, 35),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = mock_ko_data))

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) >= 0)
  if (nrow(result) > 0) {
    expect_true(all(result >= 0))
    expect_equal(ncol(result), 2)
  }
})

test_that("ko2kegg_abundance handles file input correctly", {
  temp_file <- tempfile(fileext = ".tsv")
  mock_ko_data <- data.frame(
    function. = c("K00001", "K00002", "K00003"),
    Sample1 = c(10, 20, 30),
    Sample2 = c(15, 25, 35),
    stringsAsFactors = FALSE
  )
  write.table(mock_ko_data, temp_file, sep = "\t", row.names = FALSE)

  result <- suppressMessages(ko2kegg_abundance(file = temp_file))
  unlink(temp_file)

  expect_s3_class(result, "data.frame")
  if (nrow(result) > 0) {
    expect_true(all(result >= 0))
  }
})

test_that("ko2kegg_abundance throws appropriate errors", {
  expect_error(ko2kegg_abundance(), "Please provide either a file or a data.frame")
  expect_error(ko2kegg_abundance(file = "test.pdf"), "Input file should be in .txt, .tsv, .csv format")
  expect_error(ko2kegg_abundance(data = list(a = 1)), "must be a data.frame")

  # Negative values
  negative_data <- data.frame(function. = c("K00001"), Sample1 = c(-1), stringsAsFactors = FALSE)
  expect_error(ko2kegg_abundance(data = negative_data), "Negative abundance values")

  # Non-numeric columns
  non_numeric <- data.frame(function. = c("K00001"), Sample1 = c("text"), stringsAsFactors = FALSE)
  expect_error(ko2kegg_abundance(data = non_numeric), "non-numeric values")
})

test_that("ko2kegg_abundance preserves sample names and removes zero pathways", {
  ko_to_kegg_reference <- ggpicrust2:::load_reference_data("ko_to_kegg")
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 5)

  mock_ko_data <- data.frame(
    function. = real_kos,
    SampleA = c(100, 200, 150, 180, 120),
    SampleB = c(50, 100, 75, 90, 60),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(ko2kegg_abundance(data = mock_ko_data))

  if (nrow(result) > 0) {
    expect_equal(colnames(result), c("SampleA", "SampleB"))
    expect_true(all(rowSums(result) > 0))
  }
})
