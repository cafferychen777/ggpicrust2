# Comprehensive tests for compare_gsea_daa function
library(testthat)

# Helper function to create realistic GSEA results
create_realistic_gsea_results <- function(n_pathways = 20, sig_proportion = 0.3) {
  set.seed(123)
  
  # Create pathway IDs with realistic KEGG format
  pathway_ids <- paste0("ko", sprintf("%05d", sample(10:999, n_pathways)))
  
  # Set some pathways to be significant
  n_significant <- round(n_pathways * sig_proportion)
  p_adjust_values <- c(
    runif(n_significant, 0.001, 0.049),  # Significant
    runif(n_pathways - n_significant, 0.051, 0.2)  # Non-significant
  )
  
  data.frame(
    pathway_id = pathway_ids,
    pathway_name = paste("Pathway", pathway_ids),
    size = sample(10:200, n_pathways, replace = TRUE),
    ES = rnorm(n_pathways, 0, 0.5),
    NES = rnorm(n_pathways, 0, 1.5),
    pvalue = runif(n_pathways, 0.001, 0.1),
    p.adjust = sample(p_adjust_values),  # Randomize order
    leading_edge = replicate(n_pathways, {
      genes <- paste0("K", sprintf("%05d", sample(1:10000, sample(5:20, 1))))
      paste(genes, collapse = ";")
    }),
    method = rep("fgsea", n_pathways),
    stringsAsFactors = FALSE
  )
}

# Helper function to create realistic DAA results
create_realistic_daa_results <- function(n_features = 20, sig_proportion = 0.4, overlap_with_gsea = TRUE) {
  set.seed(if(overlap_with_gsea) 123 else 456)  # Control overlap
  
  if (overlap_with_gsea) {
    # Create features that overlap with GSEA pathway_ids
    feature_ids <- paste0("ko", sprintf("%05d", sample(10:999, n_features)))
  } else {
    # Create features with no overlap
    feature_ids <- paste0("ko", sprintf("%05d", sample(1000:1999, n_features)))
  }
  
  # Set some features to be significant
  n_significant <- round(n_features * sig_proportion)
  p_adjust_values <- c(
    runif(n_significant, 0.001, 0.049),  # Significant
    runif(n_features - n_significant, 0.051, 0.3)  # Non-significant
  )
  
  data.frame(
    feature = feature_ids,
    method = rep("ALDEx2", n_features),
    group1 = rep("Control", n_features),
    group2 = rep("Treatment", n_features),
    p_values = runif(n_features, 0.001, 0.1),
    p_adjust = sample(p_adjust_values),  # Randomize order
    log_2_fold_change = rnorm(n_features, 0, 2),
    effect_size = rnorm(n_features, 0, 1),
    stringsAsFactors = FALSE
  )
}

test_that("compare_gsea_daa creates venn diagram with overlapping results", {
  skip_if_not_installed("ggplot2")
  
  # Create test data with known overlap
  gsea_results <- create_realistic_gsea_results(20, 0.4)  # 8 significant
  daa_results <- create_realistic_daa_results(20, 0.3, overlap_with_gsea = TRUE)  # 6 significant
  
  # Test basic venn diagram creation
  comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "venn")
  
  expect_type(comparison, "list")
  expect_true("plot" %in% names(comparison))
  expect_true("results" %in% names(comparison))
  expect_s3_class(comparison$plot, "ggplot")
  
  # Check comparison results structure
  expect_true("overlap" %in% names(comparison$results))
  expect_true("gsea_only" %in% names(comparison$results))
  expect_true("daa_only" %in% names(comparison$results))
  expect_true("n_overlap" %in% names(comparison$results))
  expect_true("n_gsea_only" %in% names(comparison$results))
  expect_true("n_daa_only" %in% names(comparison$results))
  expect_true("n_gsea_total" %in% names(comparison$results))
  expect_true("n_daa_total" %in% names(comparison$results))
  
  # Check that counts are consistent
  expect_equal(
    comparison$results$n_gsea_total,
    comparison$results$n_overlap + comparison$results$n_gsea_only
  )
  expect_equal(
    comparison$results$n_daa_total,
    comparison$results$n_overlap + comparison$results$n_daa_only
  )
})

test_that("compare_gsea_daa handles no overlapping pathways", {
  skip_if_not_installed("ggplot2")
  
  # Create test data with no overlap
  gsea_results <- create_realistic_gsea_results(10, 0.5)
  daa_results <- create_realistic_daa_results(10, 0.5, overlap_with_gsea = FALSE)
  
  comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "venn")
  
  expect_type(comparison, "list")
  expect_equal(comparison$results$n_overlap, 0)
  expect_equal(length(comparison$results$overlap), 0)
  expect_true(comparison$results$n_gsea_only >= 0)
  expect_true(comparison$results$n_daa_only >= 0)
})

test_that("compare_gsea_daa creates scatter plot with merged data", {
  skip_if_not_installed("ggplot2")
  
  # Create test data with guaranteed overlap
  gsea_results <- create_realistic_gsea_results(15, 0.6)
  # Ensure some DAA features match GSEA pathway_ids
  daa_results <- create_realistic_daa_results(15, 0.4, overlap_with_gsea = TRUE)
  # Force some overlap by using same IDs
  daa_results$feature[1:5] <- gsea_results$pathway_id[1:5]
  
  comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "scatter")
  
  expect_type(comparison, "list")
  expect_s3_class(comparison$plot, "ggplot")
  
  # Check that the plot has the expected structure
  plot_data <- ggplot2::ggplot_build(comparison$plot)
  expect_true(length(plot_data$data) > 0)
})

test_that("compare_gsea_daa handles scatter plot with no overlapping features", {
  skip_if_not_installed("ggplot2")
  
  # Create test data with no overlap
  gsea_results <- create_realistic_gsea_results(10)
  daa_results <- create_realistic_daa_results(10, overlap_with_gsea = FALSE)
  
  expect_warning(
    comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "scatter"),
    "No overlapping pathways found for scatter plot"
  )
  
  expect_s3_class(comparison$plot, "ggplot")
})

test_that("compare_gsea_daa creates upset plot", {
  skip_if_not_installed("ggplot2")
  
  gsea_results <- create_realistic_gsea_results(12, 0.5)
  daa_results <- create_realistic_daa_results(12, 0.5, overlap_with_gsea = TRUE)
  
  comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "upset")
  
  expect_type(comparison, "list")
  expect_true("plot" %in% names(comparison))
  expect_true("results" %in% names(comparison))
})

test_that("compare_gsea_daa handles heatmap plot (not implemented)", {
  skip_if_not_installed("ggplot2")
  
  gsea_results <- create_realistic_gsea_results(10)
  daa_results <- create_realistic_daa_results(10)
  
  expect_warning(
    comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "heatmap"),
    "Heatmap plot not yet implemented"
  )
  
  expect_s3_class(comparison$plot, "ggplot")
})

test_that("compare_gsea_daa respects different p_threshold values", {
  skip_if_not_installed("ggplot2")
  
  # Create test data with known p-values
  # Note: compare_gsea_daa uses strict inequality (<), not <=
  gsea_results <- data.frame(
    pathway_id = paste0("ko", sprintf("%05d", 1:10)),
    pathway_name = paste("Pathway", 1:10),
    size = sample(10:100, 10, replace = TRUE),
    ES = rnorm(10, 0, 0.5),
    NES = rnorm(10, 0, 1.5),
    pvalue = runif(10, 0.001, 0.1),
    p.adjust = c(0.005, 0.02, 0.03, 0.06, 0.08, 0.10, 0.12, 0.15, 0.18, 0.19),
    leading_edge = replicate(10, paste(sample(LETTERS, 5), collapse = ";")),
    method = rep("fgsea", 10),
    stringsAsFactors = FALSE
  )

  daa_results <- data.frame(
    feature = paste0("ko", sprintf("%05d", 1:10)),
    method = rep("ALDEx2", 10),
    group1 = rep("Control", 10),
    group2 = rep("Treatment", 10),
    p_values = runif(10, 0.001, 0.1),
    p_adjust = c(0.005, 0.015, 0.025, 0.04, 0.07, 0.09, 0.11, 0.14, 0.17, 0.19),
    log_2_fold_change = rnorm(10, 0, 2),
    stringsAsFactors = FALSE
  )

  # Test with strict threshold
  comparison_strict <- compare_gsea_daa(gsea_results, daa_results, p_threshold = 0.01)

  # Test with lenient threshold
  comparison_lenient <- compare_gsea_daa(gsea_results, daa_results, p_threshold = 0.2)

  # Strict threshold should have fewer significant pathways
  expect_true(comparison_strict$results$n_gsea_total <= comparison_lenient$results$n_gsea_total)
  expect_true(comparison_strict$results$n_daa_total <= comparison_lenient$results$n_daa_total)

  # With p_threshold = 0.01 (strict <):
  # GSEA: only 0.005 < 0.01, so 1 significant
  # DAA: only 0.005 < 0.01, so 1 significant
  expect_equal(comparison_strict$results$n_gsea_total, 1)
  expect_equal(comparison_strict$results$n_daa_total, 1)

  # With p_threshold = 0.2: all values are < 0.2
  expect_equal(comparison_lenient$results$n_gsea_total, 10)
  expect_equal(comparison_lenient$results$n_daa_total, 10)
})

test_that("compare_gsea_daa validates input parameters thoroughly", {
  gsea_results <- create_realistic_gsea_results(5)
  daa_results <- create_realistic_daa_results(5)
  
  # Test invalid gsea_results types
  expect_error(
    compare_gsea_daa(gsea_results = list(), daa_results = daa_results),
    "'gsea_results' must be a data frame"
  )
  
  expect_error(
    compare_gsea_daa(gsea_results = matrix(1:10, ncol = 2), daa_results = daa_results),
    "'gsea_results' must be a data frame"
  )
  
  # Test invalid daa_results types
  expect_error(
    compare_gsea_daa(gsea_results = gsea_results, daa_results = "invalid"),
    "'daa_results' must be a data frame"
  )
  
  # Test invalid plot_type
  expect_error(
    compare_gsea_daa(gsea_results, daa_results, plot_type = "invalid"),
    "plot_type must be one of 'venn', 'upset', 'scatter', or 'heatmap'"
  )
  
  expect_error(
    compare_gsea_daa(gsea_results, daa_results, plot_type = c("venn", "upset")),
    "plot_type must be one of 'venn', 'upset', 'scatter', or 'heatmap'"
  )
  
  # Test missing required columns in gsea_results
  gsea_missing_pathway_id <- gsea_results[, !names(gsea_results) %in% "pathway_id"]
  expect_error(
    compare_gsea_daa(gsea_missing_pathway_id, daa_results),
    "GSEA results missing required columns: pathway_id, p.adjust"
  )
  
  gsea_missing_p_adjust <- gsea_results[, !names(gsea_results) %in% "p.adjust"]
  expect_error(
    compare_gsea_daa(gsea_missing_p_adjust, daa_results),
    "GSEA results missing required columns: pathway_id, p.adjust"
  )
  
  # Test missing required columns in daa_results
  daa_missing_feature <- daa_results[, !names(daa_results) %in% "feature"]
  expect_error(
    compare_gsea_daa(gsea_results, daa_missing_feature),
    "DAA results missing required columns: feature, p_adjust"
  )
  
  daa_missing_p_adjust <- daa_results[, !names(daa_results) %in% "p_adjust"]
  expect_error(
    compare_gsea_daa(gsea_results, daa_missing_p_adjust),
    "DAA results missing required columns: feature, p_adjust"
  )
})

test_that("compare_gsea_daa handles empty results", {
  skip_if_not_installed("ggplot2")
  
  # Create empty data frames with proper structure
  empty_gsea <- data.frame(
    pathway_id = character(0),
    pathway_name = character(0),
    size = integer(0),
    ES = numeric(0),
    NES = numeric(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    leading_edge = character(0),
    method = character(0),
    stringsAsFactors = FALSE
  )
  
  empty_daa <- data.frame(
    feature = character(0),
    method = character(0),
    group1 = character(0),
    group2 = character(0),
    p_values = numeric(0),
    p_adjust = numeric(0),
    log_2_fold_change = numeric(0),
    stringsAsFactors = FALSE
  )
  
  comparison <- compare_gsea_daa(empty_gsea, empty_daa)
  
  expect_type(comparison, "list")
  expect_equal(comparison$results$n_overlap, 0)
  expect_equal(comparison$results$n_gsea_total, 0)
  expect_equal(comparison$results$n_daa_total, 0)
})

test_that("compare_gsea_daa handles data with NA values", {
  skip_if_not_installed("ggplot2")
  
  gsea_results <- create_realistic_gsea_results(10)
  daa_results <- create_realistic_daa_results(10)
  
  # Introduce NA values
  gsea_results$p.adjust[1] <- NA
  daa_results$p_adjust[2] <- NA
  gsea_results$pathway_id[3] <- NA
  daa_results$feature[4] <- NA
  
  # Should still work but might affect overlap calculations
  comparison <- compare_gsea_daa(gsea_results, daa_results)
  
  expect_type(comparison, "list")
  expect_s3_class(comparison$plot, "ggplot")
})

test_that("compare_gsea_daa set operations work correctly", {
  # Create controlled test data to verify set operations
  gsea_pathways <- c("ko00001", "ko00002", "ko00003", "ko00004", "ko00005")
  daa_features <- c("ko00003", "ko00004", "ko00005", "ko00006", "ko00007")
  
  gsea_results <- data.frame(
    pathway_id = gsea_pathways,
    pathway_name = paste("Pathway", gsea_pathways),
    size = rep(50, 5),
    ES = rep(0.5, 5),
    NES = rep(1.5, 5),
    pvalue = rep(0.01, 5),
    p.adjust = rep(0.02, 5),  # All significant at 0.05
    leading_edge = rep("K00001;K00002", 5),
    method = rep("fgsea", 5),
    stringsAsFactors = FALSE
  )
  
  daa_results <- data.frame(
    feature = daa_features,
    method = rep("ALDEx2", 5),
    group1 = rep("Control", 5),
    group2 = rep("Treatment", 5),
    p_values = rep(0.01, 5),
    p_adjust = rep(0.03, 5),  # All significant at 0.05
    log_2_fold_change = rep(1.0, 5),
    stringsAsFactors = FALSE
  )
  
  comparison <- compare_gsea_daa(gsea_results, daa_results, p_threshold = 0.05)
  
  # Expected results:
  # GSEA significant: ko00001, ko00002, ko00003, ko00004, ko00005 (5 total)
  # DAA significant: ko00003, ko00004, ko00005, ko00006, ko00007 (5 total)
  # Overlap: ko00003, ko00004, ko00005 (3)
  # GSEA only: ko00001, ko00002 (2)
  # DAA only: ko00006, ko00007 (2)
  
  expect_equal(comparison$results$n_gsea_total, 5)
  expect_equal(comparison$results$n_daa_total, 5)
  expect_equal(comparison$results$n_overlap, 3)
  expect_equal(comparison$results$n_gsea_only, 2)
  expect_equal(comparison$results$n_daa_only, 2)
  
  expect_setequal(comparison$results$overlap, c("ko00003", "ko00004", "ko00005"))
  expect_setequal(comparison$results$gsea_only, c("ko00001", "ko00002"))
  expect_setequal(comparison$results$daa_only, c("ko00006", "ko00007"))
})

test_that("compare_gsea_daa handles different column types", {
  skip_if_not_installed("ggplot2")
  
  gsea_results <- create_realistic_gsea_results(5)
  daa_results <- create_realistic_daa_results(5)
  
  # Convert some columns to factors
  gsea_results$pathway_id <- as.factor(gsea_results$pathway_id)
  daa_results$feature <- as.factor(daa_results$feature)
  
  comparison <- compare_gsea_daa(gsea_results, daa_results)
  
  expect_type(comparison, "list")
  expect_s3_class(comparison$plot, "ggplot")
})

test_that("compare_gsea_daa fallback plots work when packages missing", {
  skip_if_not_installed("ggplot2")
  
  gsea_results <- create_realistic_gsea_results(5)
  daa_results <- create_realistic_daa_results(5)
  
  # Mock package availability checks to return FALSE
  with_mock(
    `requireNamespace` = function(pkg, quietly = TRUE) {
      if (pkg %in% c("ggVennDiagram", "UpSetR")) return(FALSE)
      return(TRUE)
    },
    {
      # Test venn fallback
      expect_warning(
        comparison_venn <- compare_gsea_daa(gsea_results, daa_results, plot_type = "venn"),
        "Package 'ggVennDiagram' is required for Venn diagrams"
      )
      expect_s3_class(comparison_venn$plot, "ggplot")
      
      # Test upset fallback
      expect_warning(
        comparison_upset <- compare_gsea_daa(gsea_results, daa_results, plot_type = "upset"),
        "Package 'UpSetR' is required for UpSet plots"
      )
      expect_s3_class(comparison_upset$plot, "ggplot")
    }
  )
})

test_that("compare_gsea_daa performance with large datasets", {
  skip_if_not_installed("ggplot2")
  
  # Create larger datasets
  large_gsea <- create_realistic_gsea_results(500, 0.1)
  large_daa <- create_realistic_daa_results(500, 0.1, overlap_with_gsea = TRUE)
  
  start_time <- Sys.time()
  comparison <- compare_gsea_daa(large_gsea, large_daa)
  end_time <- Sys.time()
  
  # Should complete in reasonable time (less than 10 seconds)
  expect_true(as.numeric(end_time - start_time, units = "secs") < 10)
  
  expect_type(comparison, "list")
  expect_s3_class(comparison$plot, "ggplot")
})