# Tests for pathway_annotation function

quiet_run <- function(expr) {
  suppressMessages(suppressWarnings(expr))
}

test_that("pathway_annotation basic functionality works", {
  temp_file <- tempfile(fileext = ".tsv")
  test_data <- data.frame(
    feature = c("K00001", "K00002"),
    sample1 = c(1, 2),
    sample2 = c(3, 4)
  )
  quiet_run(write.table(test_data, temp_file, sep = "\t", row.names = FALSE))

  result <- quiet_run(pathway_annotation(file = temp_file, pathway = "KO", ko_to_kegg = FALSE))

  expect_s3_class(result, "data.frame")
  expect_true("description" %in% colnames(result))
  expect_equal(nrow(result), 2)

  unlink(temp_file)
})

test_that("pathway_annotation works with daa_results_df input", {
  test_daa_df <- data.frame(
    feature = c("K00001", "K00002"),
    p_adjust = c(0.04, 0.03),
    stringsAsFactors = FALSE
  )

  result <- quiet_run(pathway_annotation(pathway = "KO", daa_results_df = test_daa_df, ko_to_kegg = FALSE))

  expect_s3_class(result, "data.frame")
  expect_true("description" %in% colnames(result))
  expect_equal(nrow(result), 2)
})

test_that("pathway_annotation works with ko_to_kegg", {
  skip_if_offline()

  test_daa_df <- data.frame(
    feature = c("K00001", "K00002"),
    p_adjust = c(0.04, 0.03),
    stringsAsFactors = FALSE
  )

  result <- quiet_run(pathway_annotation(pathway = "KO", daa_results_df = test_daa_df, ko_to_kegg = TRUE))

  expect_s3_class(result, "data.frame")
  expect_true(all(c("pathway_name", "pathway_class") %in% colnames(result)))
})

test_that("pathway_annotation works with all pathway types", {
  test_daa_df <- data.frame(
    feature = c("K00001", "EC:1.1.1.1", "PWY-7219"),
    p_adjust = c(0.04, 0.03, 0.02),
    stringsAsFactors = FALSE
  )

  ko_result <- quiet_run(pathway_annotation(pathway = "KO", daa_results_df = test_daa_df, ko_to_kegg = FALSE))
  expect_s3_class(ko_result, "data.frame")

  ec_result <- quiet_run(pathway_annotation(pathway = "EC", daa_results_df = test_daa_df, ko_to_kegg = FALSE))
  expect_s3_class(ec_result, "data.frame")

  metacyc_result <- quiet_run(pathway_annotation(pathway = "MetaCyc", daa_results_df = test_daa_df, ko_to_kegg = FALSE))
  expect_s3_class(metacyc_result, "data.frame")
})

test_that("pathway_annotation validates inputs", {
  expect_error(pathway_annotation(), "Please input")
  expect_error(pathway_annotation(file = "test.tsv", pathway = "INVALID"), "Unknown reference type")

  empty_df <- data.frame(feature = character(0), p_adjust = numeric(0))
  expect_error(quiet_run(pathway_annotation(pathway = "KO", daa_results_df = empty_df)), "empty")
})
