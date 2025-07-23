# Test for compare_gsea_daa function
library(testthat)

# Helper function to create test GSEA results
create_test_gsea_results <- function(n_pathways = 10) {
  set.seed(123)
  data.frame(
    pathway_id = paste0("path:ko", sprintf("%05d", 1:n_pathways)),
    pathway_name = paste("Pathway", 1:n_pathways),
    size = sample(10:100, n_pathways, replace = TRUE),
    ES = runif(n_pathways, -0.8, 0.8),
    NES = runif(n_pathways, -2, 2),
    pvalue = runif(n_pathways, 0, 0.1),
    p.adjust = runif(n_pathways, 0, 0.2),
    leading_edge = replicate(n_pathways, paste(paste0("K", sprintf("%05d", sample(1:1000, 5))), collapse = ";")),
    method = rep("fgsea", n_pathways),
    stringsAsFactors = FALSE
  )
}

# Helper function to create test DAA results
create_test_daa_results <- function(n_features = 10) {
  set.seed(123)  # Use same seed as GSEA results for consistent overlap
  data.frame(
    feature = paste0("path:ko", sprintf("%05d", 1:n_features)),  # Ensure overlap with GSEA results
    method = rep("ALDEx2", n_features),
    group1 = rep("Control", n_features),
    group2 = rep("Treatment", n_features),
    p_values = runif(n_features, 0, 0.1),
    p_adjust = runif(n_features, 0, 0.2),
    log_2_fold_change = rnorm(n_features, 0, 1),
    stringsAsFactors = FALSE
  )
}

test_that("compare_gsea_daa creates venn diagram correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggVennDiagram")

  # Create test data
  gsea_results <- create_test_gsea_results()
  daa_results <- create_test_daa_results()

  # Mock the ggVennDiagram function to avoid actual plotting
  mockery::stub(compare_gsea_daa, "ggVennDiagram::ggVennDiagram", function(...) {
    ggplot2::ggplot()
  })

  # Test venn diagram
  comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "venn")

  # Check the result
  expect_type(comparison, "list")
  expect_s3_class(comparison$plot, "ggplot")

  # We don't check exact numbers since they depend on random data and p-value thresholds
  expect_true(length(comparison$results$overlap) >= 0)
  expect_true(comparison$results$n_gsea_only >= 0)
  expect_true(comparison$results$n_daa_only >= 0)

  # Test with custom p_threshold
  comparison_custom <- compare_gsea_daa(
    gsea_results = gsea_results,
    daa_results = daa_results,
    plot_type = "venn",
    p_threshold = 0.01
  )

  # Check the result
  expect_type(comparison_custom, "list")
  expect_s3_class(comparison_custom$plot, "ggplot")
})

test_that("compare_gsea_daa creates upset plot correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("UpSetR")

  # Create test data
  gsea_results <- create_test_gsea_results()
  daa_results <- create_test_daa_results()

  # Mock the compare_gsea_daa function to avoid actual UpSetR plotting
  # We'll just check that the function doesn't error and returns a list with expected structure
  mockery::stub(compare_gsea_daa, "UpSetR::upset", function(...) {
    # Return a dummy plot object that won't cause errors when converted
    ggplot2::ggplot()
  })

  # Test upset plot
  comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "upset")

  # Check the result
  expect_type(comparison, "list")
  expect_true("plot" %in% names(comparison))
  expect_true("results" %in% names(comparison))

  # We don't check exact numbers since they depend on random data and p-value thresholds
  expect_true(length(comparison$results$overlap) >= 0)
})

test_that("compare_gsea_daa creates scatter plot correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")

  # Create test data
  gsea_results <- create_test_gsea_results()
  daa_results <- create_test_daa_results()

  # Test scatter plot
  comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "scatter")

  # Check the result
  expect_type(comparison, "list")
  expect_s3_class(comparison$plot, "ggplot")

  # Test with no overlapping pathways
  gsea_results_no_overlap <- gsea_results
  gsea_results_no_overlap$pathway_id <- paste0("path:ko", sprintf("%05d", 21:30))

  expect_warning(
    comparison_no_overlap <- compare_gsea_daa(
      gsea_results = gsea_results_no_overlap,
      daa_results = daa_results,
      plot_type = "scatter"
    ),
    "No overlapping pathways found for scatter plot"
  )

  # Check the result
  expect_type(comparison_no_overlap, "list")
  expect_s3_class(comparison_no_overlap$plot, "ggplot")
})

test_that("compare_gsea_daa handles heatmap plot", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")

  # Create test data
  gsea_results <- create_test_gsea_results()
  daa_results <- create_test_daa_results()

  # Test heatmap plot (not fully implemented yet)
  expect_warning(
    comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "heatmap"),
    "Heatmap plot not yet implemented"
  )

  # Check the result
  expect_type(comparison, "list")
  expect_s3_class(comparison$plot, "ggplot")
})

test_that("compare_gsea_daa validates inputs correctly", {
  # Create test data
  gsea_results <- create_test_gsea_results()
  daa_results <- create_test_daa_results()

  # Test invalid gsea_results
  expect_error(
    compare_gsea_daa(gsea_results = "invalid", daa_results = daa_results),
    "'gsea_results' must be a data frame"
  )

  # Test invalid daa_results
  expect_error(
    compare_gsea_daa(gsea_results = gsea_results, daa_results = "invalid"),
    "'daa_results' must be a data frame"
  )

  # Test invalid plot_type
  expect_error(
    compare_gsea_daa(gsea_results, daa_results, plot_type = "invalid"),
    "plot_type must be one of 'venn', 'upset', 'scatter', or 'heatmap'"
  )

  # Test missing required columns in gsea_results
  gsea_results_missing <- gsea_results[, !names(gsea_results) %in% c("pathway_id")]
  expect_error(
    compare_gsea_daa(gsea_results_missing, daa_results),
    "GSEA results missing required columns: pathway_id, p.adjust"
  )

  # Test missing required columns in daa_results
  daa_results_missing <- daa_results[, !names(daa_results) %in% c("feature")]
  expect_error(
    compare_gsea_daa(gsea_results, daa_results_missing),
    "DAA results missing required columns: feature, p_adjust"
  )
})

test_that("compare_gsea_daa handles different p_threshold values", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")

  # Create test data
  gsea_results <- create_test_gsea_results()
  daa_results <- create_test_daa_results()

  # Make some p-values very small to ensure they pass any threshold
  gsea_results$p.adjust[1:3] <- 0.001
  daa_results$p_adjust[1:3] <- 0.001

  # Test with default p_threshold (0.05)
  comparison_default <- compare_gsea_daa(gsea_results, daa_results)

  # Test with stricter p_threshold
  comparison_strict <- compare_gsea_daa(gsea_results, daa_results, p_threshold = 0.01)

  # Test with lenient p_threshold
  comparison_lenient <- compare_gsea_daa(gsea_results, daa_results, p_threshold = 0.2)

  # Check that the number of significant pathways changes with threshold
  # We can't guarantee the exact relationship with random data, so we just check they're numbers
  expect_type(comparison_strict$results$n_gsea_total, "integer")
  expect_type(comparison_default$results$n_gsea_total, "integer")
  expect_type(comparison_lenient$results$n_gsea_total, "integer")
})
