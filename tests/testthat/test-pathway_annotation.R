# Tests for pathway_annotation function

test_that("pathway_annotation basic functionality works", {
  temp_file <- tempfile(fileext = ".tsv")
  test_data <- data.frame(
    feature = c("K00001", "K00002"),
    sample1 = c(1, 2),
    sample2 = c(3, 4)
  )
  suppressMessages(write.table(test_data, temp_file, sep = "\t", row.names = FALSE))

  result <- suppressMessages(pathway_annotation(file = temp_file, pathway = "KO", ko_to_kegg = FALSE))

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

  result <- suppressMessages(pathway_annotation(pathway = "KO", daa_results_df = test_daa_df, ko_to_kegg = FALSE))

  expect_s3_class(result, "data.frame")
  expect_true("description" %in% colnames(result))
  expect_equal(nrow(result), 2)
})

test_that("pathway_annotation works with ko_to_kegg", {
  skip_if(
    Sys.getenv("GGPICRUST2_RUN_NETWORK_TESTS", "false") != "true",
    "Set GGPICRUST2_RUN_NETWORK_TESTS=true to run network-dependent KEGG tests."
  )
  skip_if_offline()

  test_daa_df <- data.frame(
    feature = c("K00001", "K00002"),
    p_adjust = c(0.04, 0.03),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(pathway_annotation(pathway = "KO", daa_results_df = test_daa_df, ko_to_kegg = TRUE))

  expect_s3_class(result, "data.frame")
  expect_true(all(c("pathway_name", "pathway_class") %in% colnames(result)))
})

test_that("pathway_annotation works with all pathway types", {
  # Each annotator gets its own correctly-typed features
  ko_df <- data.frame(feature = c("K00001", "K00002"), p_adjust = c(.04, .03), stringsAsFactors = FALSE)
  ec_df <- data.frame(feature = c("EC:1.1.1.1", "EC:2.7.1.1"), p_adjust = c(.04, .03), stringsAsFactors = FALSE)
  metacyc_df <- data.frame(feature = c("PWY-7219", "GLYCOLYSIS"), p_adjust = c(.04, .03), stringsAsFactors = FALSE)

  ko_result <- suppressMessages(pathway_annotation(pathway = "KO", daa_results_df = ko_df, ko_to_kegg = FALSE))
  expect_s3_class(ko_result, "data.frame")

  ec_result <- suppressMessages(pathway_annotation(pathway = "EC", daa_results_df = ec_df, ko_to_kegg = FALSE))
  expect_s3_class(ec_result, "data.frame")

  metacyc_result <- suppressMessages(pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_df, ko_to_kegg = FALSE))
  expect_s3_class(metacyc_result, "data.frame")
})

test_that("pathway_annotation validates inputs", {
  expect_error(pathway_annotation(), "Please input")

  # Test invalid pathway type through daa_results_df path (avoids file-existence check)
  test_df <- data.frame(feature = "K00001", p_adjust = 0.04, stringsAsFactors = FALSE)
  expect_error(
    suppressMessages(pathway_annotation(pathway = "INVALID", daa_results_df = test_df, ko_to_kegg = FALSE)),
    "Unknown reference type"
  )

  empty_df <- data.frame(feature = character(0), p_adjust = numeric(0))
  expect_error(suppressMessages(pathway_annotation(pathway = "KO", daa_results_df = empty_df)), "empty")
})
