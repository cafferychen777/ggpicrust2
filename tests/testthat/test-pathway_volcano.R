# Tests for pathway_volcano function

test_that("pathway_volcano creates a ggplot object", {
  skip_if_not_installed("ggrepel")

  # Create mock DAA results
  daa_results <- data.frame(
    feature = paste0("ko", sprintf("%05d", 1:20)),
    pathway_name = paste0("Pathway ", 1:20),
    log2FoldChange = rnorm(20, 0, 2),
    p_adjust = runif(20, 0, 0.1),
    stringsAsFactors = FALSE
  )

  p <- pathway_volcano(daa_results)

  expect_s3_class(p, "ggplot")
})

test_that("pathway_volcano handles custom column names", {
  skip_if_not_installed("ggrepel")

  # Create mock DAA results with custom column names
  daa_results <- data.frame(
    id = paste0("ko", sprintf("%05d", 1:20)),
    name = paste0("Pathway ", 1:20),
    fc = rnorm(20, 0, 2),
    pval = runif(20, 0, 0.1),
    stringsAsFactors = FALSE
  )

  p <- pathway_volcano(
    daa_results,
    fc_col = "fc",
    p_col = "pval",
    label_col = "name"
  )

  expect_s3_class(p, "ggplot")
})

test_that("pathway_volcano handles custom thresholds", {
  skip_if_not_installed("ggrepel")

  daa_results <- data.frame(
    feature = paste0("ko", sprintf("%05d", 1:20)),
    pathway_name = paste0("Pathway ", 1:20),
    log2FoldChange = rnorm(20, 0, 2),
    p_adjust = runif(20, 0, 0.1),
    stringsAsFactors = FALSE
  )

  p <- pathway_volcano(
    daa_results,
    fc_threshold = 0.5,
    p_threshold = 0.01,
    label_top_n = 5
  )

  expect_s3_class(p, "ggplot")
})

test_that("pathway_volcano handles custom colors", {
  skip_if_not_installed("ggrepel")

  daa_results <- data.frame(
    feature = paste0("ko", sprintf("%05d", 1:20)),
    pathway_name = paste0("Pathway ", 1:20),
    log2FoldChange = rnorm(20, 0, 2),
    p_adjust = runif(20, 0, 0.1),
    stringsAsFactors = FALSE
  )

  p <- pathway_volcano(
    daa_results,
    colors = c("Down" = "blue", "Not Significant" = "gray", "Up" = "red")
  )

  expect_s3_class(p, "ggplot")
})

test_that("pathway_volcano handles no labels", {
  skip_if_not_installed("ggrepel")

  daa_results <- data.frame(
    feature = paste0("ko", sprintf("%05d", 1:20)),
    pathway_name = paste0("Pathway ", 1:20),
    log2FoldChange = rnorm(20, 0, 2),
    p_adjust = runif(20, 0, 0.1),
    stringsAsFactors = FALSE
  )

  p <- pathway_volcano(
    daa_results,
    label_top_n = 0
  )

  expect_s3_class(p, "ggplot")
})

test_that("pathway_volcano errors on missing columns", {
  daa_results <- data.frame(
    feature = paste0("ko", sprintf("%05d", 1:20)),
    some_value = rnorm(20),
    stringsAsFactors = FALSE
  )

  expect_error(
    pathway_volcano(daa_results),
    "not found"
  )
})

test_that("pathway_volcano errors on non-data.frame input", {
  expect_error(
    pathway_volcano("not a data frame"),
    "must be a data frame"
  )
})

test_that("pathway_volcano works with real DAA workflow", {
  skip_if_not_installed("ggrepel")
  skip_on_cran()

  # Load example data
  data(ko_abundance, package = "ggpicrust2")
  data(metadata, package = "ggpicrust2")

  # Run the full DAA workflow
  kegg_abundance <- ko2kegg_abundance(data = ko_abundance)

  daa_results <- suppressMessages(pathway_daa(
    abundance = kegg_abundance,
    metadata = metadata,
    group = "Environment",
    daa_method = "LinDA"
  ))

  daa_annotated <- suppressMessages(pathway_annotation(
    pathway = "KO",
    ko_to_kegg = TRUE,
    daa_results_df = daa_results
  ))

  p <- pathway_volcano(daa_annotated, label_top_n = 5)

  expect_s3_class(p, "ggplot")
})
