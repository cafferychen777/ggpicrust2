# Advanced comprehensive tests for compare_gsea_daa function
library(testthat)

# Advanced helper functions for creating controlled test scenarios
create_controlled_gsea_results <- function(significant_pathways, non_significant_pathways, 
                                          p_threshold = 0.05, effect_sizes = NULL) {
  all_pathways <- c(significant_pathways, non_significant_pathways)
  n_total <- length(all_pathways)
  
  # Set p.adjust values based on significance
  p_adjust_values <- numeric(n_total)
  p_adjust_values[1:length(significant_pathways)] <- runif(length(significant_pathways), 0.001, p_threshold - 0.001)
  p_adjust_values[(length(significant_pathways) + 1):n_total] <- runif(length(non_significant_pathways), p_threshold + 0.001, 0.2)
  
  # Set effect sizes if provided
  if (is.null(effect_sizes)) {
    effect_sizes <- rnorm(n_total, 0, 1.5)
  }
  
  data.frame(
    pathway_id = all_pathways,
    pathway_name = paste("Pathway", all_pathways),
    size = sample(15:150, n_total, replace = TRUE),
    ES = effect_sizes,
    NES = effect_sizes * 1.2 + rnorm(n_total, 0, 0.3),
    pvalue = runif(n_total, 0.001, 0.1),
    p.adjust = p_adjust_values,
    leading_edge = replicate(n_total, {
      genes <- paste0("K", sprintf("%05d", sample(1:15000, sample(5:25, 1))))
      paste(genes, collapse = ";")
    }),
    method = rep("fgsea", n_total),
    stringsAsFactors = FALSE
  )
}

create_controlled_daa_results <- function(significant_features, non_significant_features, 
                                         p_threshold = 0.05, fold_changes = NULL) {
  all_features <- c(significant_features, non_significant_features)
  n_total <- length(all_features)
  
  # Set p_adjust values based on significance
  p_adjust_values <- numeric(n_total)
  p_adjust_values[1:length(significant_features)] <- runif(length(significant_features), 0.001, p_threshold - 0.001)
  p_adjust_values[(length(significant_features) + 1):n_total] <- runif(length(non_significant_features), p_threshold + 0.001, 0.25)
  
  # Set fold changes if provided
  if (is.null(fold_changes)) {
    fold_changes <- rnorm(n_total, 0, 2)
  }
  
  data.frame(
    feature = all_features,
    method = rep("ALDEx2", n_total),
    group1 = rep("Control", n_total),
    group2 = rep("Treatment", n_total),
    p_values = runif(n_total, 0.001, 0.1),
    p_adjust = p_adjust_values,
    log_2_fold_change = fold_changes,
    effect_size = fold_changes * 0.8 + rnorm(n_total, 0, 0.4),
    stringsAsFactors = FALSE
  )
}

test_that("compare_gsea_daa identifies significant pathway overlaps accurately", {
  skip_if_not_installed("ggplot2")
  
  # Define controlled overlap scenario
  common_pathways <- c("ko00010", "ko00020", "ko00030", "ko00040", "ko00050")
  gsea_only_pathways <- c("ko00060", "ko00070", "ko00080")
  daa_only_features <- c("ko00090", "ko00100", "ko00110", "ko00120")
  
  # All common pathways should be significant in both analyses
  gsea_results <- create_controlled_gsea_results(
    significant_pathways = c(common_pathways, gsea_only_pathways),
    non_significant_pathways = c("ko00200", "ko00210", "ko00220"),
    p_threshold = 0.05
  )
  
  daa_results <- create_controlled_daa_results(
    significant_features = c(common_pathways, daa_only_features),
    non_significant_features = c("ko00300", "ko00310", "ko00320"),
    p_threshold = 0.05
  )
  
  # Test with default threshold
  comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "venn")
  
  # Verify exact overlap identification
  expect_setequal(comparison$results$overlap, common_pathways)
  expect_setequal(comparison$results$gsea_only, gsea_only_pathways)
  expect_setequal(comparison$results$daa_only, daa_only_features)
  
  # Verify counts
  expect_equal(comparison$results$n_overlap, length(common_pathways))
  expect_equal(comparison$results$n_gsea_only, length(gsea_only_pathways))
  expect_equal(comparison$results$n_daa_only, length(daa_only_features))
  expect_equal(comparison$results$n_gsea_total, length(common_pathways) + length(gsea_only_pathways))
  expect_equal(comparison$results$n_daa_total, length(common_pathways) + length(daa_only_features))
})

test_that("compare_gsea_daa overlap analysis responds correctly to p_threshold changes", {
  skip_if_not_installed("ggplot2")
  
  # Create data with known p-value distributions
  pathways <- paste0("ko", sprintf("%05d", 1:20))
  
  # Create GSEA results with specific p.adjust values
  gsea_results <- data.frame(
    pathway_id = pathways,
    pathway_name = paste("Pathway", pathways),
    size = rep(50, 20),
    ES = rnorm(20, 0, 0.5),
    NES = rnorm(20, 0, 1.2),
    pvalue = runif(20, 0.001, 0.1),
    p.adjust = c(0.005, 0.01, 0.02, 0.03, 0.04,    # 5 very significant
                 0.06, 0.07, 0.08, 0.09, 0.095,     # 5 moderately significant
                 0.11, 0.12, 0.13, 0.14, 0.15,      # 5 not significant
                 0.16, 0.17, 0.18, 0.19, 0.20),     # 5 not significant
    leading_edge = rep("K00001;K00002", 20),
    method = rep("fgsea", 20),
    stringsAsFactors = FALSE
  )
  
  # Create DAA results with specific p_adjust values  
  daa_results <- data.frame(
    feature = pathways,
    method = rep("ALDEx2", 20),
    group1 = rep("Control", 20),
    group2 = rep("Treatment", 20),
    p_values = runif(20, 0.001, 0.1),
    p_adjust = c(0.003, 0.008, 0.015, 0.025, 0.035,  # 5 very significant
                 0.055, 0.065, 0.075, 0.085, 0.092,   # 5 moderately significant  
                 0.105, 0.115, 0.125, 0.135, 0.145,   # 5 not significant
                 0.155, 0.165, 0.175, 0.185, 0.195),  # 5 not significant
    log_2_fold_change = rnorm(20, 0, 1.5),
    stringsAsFactors = FALSE
  )
  
  # Test with very strict threshold (p < 0.01)
  comparison_strict <- compare_gsea_daa(gsea_results, daa_results, p_threshold = 0.01)
  # GSEA: pathways 1-2 significant (p.adjust <= 0.01)
  # DAA: pathways 1-2 significant (p_adjust <= 0.01) 
  expect_equal(comparison_strict$results$n_gsea_total, 2)
  expect_equal(comparison_strict$results$n_daa_total, 2) 
  expect_equal(comparison_strict$results$n_overlap, 2)
  
  # Test with standard threshold (p < 0.05)
  comparison_standard <- compare_gsea_daa(gsea_results, daa_results, p_threshold = 0.05)
  # GSEA: pathways 1-5 significant (p.adjust <= 0.05)
  # DAA: pathways 1-5 significant (p_adjust <= 0.05)
  expect_equal(comparison_standard$results$n_gsea_total, 5)
  expect_equal(comparison_standard$results$n_daa_total, 5)
  expect_equal(comparison_standard$results$n_overlap, 5)
  
  # Test with lenient threshold (p < 0.1)
  comparison_lenient <- compare_gsea_daa(gsea_results, daa_results, p_threshold = 0.1)
  # GSEA: pathways 1-10 significant (p.adjust <= 0.1)
  # DAA: pathways 1-10 significant (p_adjust <= 0.1)
  expect_equal(comparison_lenient$results$n_gsea_total, 10)
  expect_equal(comparison_lenient$results$n_daa_total, 10)
  expect_equal(comparison_lenient$results$n_overlap, 10)
  
  # Verify the progression makes sense
  expect_true(comparison_strict$results$n_overlap <= comparison_standard$results$n_overlap)
  expect_true(comparison_standard$results$n_overlap <= comparison_lenient$results$n_overlap)
})

test_that("compare_gsea_daa venn diagram visualization handles different overlap scenarios", {
  skip_if_not_installed("ggplot2")
  
  # Test scenario 1: High overlap
  high_overlap_gsea <- create_controlled_gsea_results(
    significant_pathways = c("ko00010", "ko00020", "ko00030", "ko00040", "ko00050"),
    non_significant_pathways = c("ko00100")
  )
  high_overlap_daa <- create_controlled_daa_results(
    significant_features = c("ko00010", "ko00020", "ko00030", "ko00040", "ko00060"),
    non_significant_features = c("ko00200")
  )
  
  comparison_high <- compare_gsea_daa(high_overlap_gsea, high_overlap_daa, plot_type = "venn")
  expect_equal(comparison_high$results$n_overlap, 4)  # ko00010-ko00040
  expect_s3_class(comparison_high$plot, "ggplot")
  
  # Test scenario 2: No overlap
  no_overlap_gsea <- create_controlled_gsea_results(
    significant_pathways = c("ko00010", "ko00020"),
    non_significant_pathways = c("ko00100")
  )
  no_overlap_daa <- create_controlled_daa_results(
    significant_features = c("ko00030", "ko00040"),
    non_significant_features = c("ko00200")
  )
  
  comparison_none <- compare_gsea_daa(no_overlap_gsea, no_overlap_daa, plot_type = "venn")
  expect_equal(comparison_none$results$n_overlap, 0)
  expect_equal(comparison_none$results$n_gsea_only, 2)
  expect_equal(comparison_none$results$n_daa_only, 2)
  expect_s3_class(comparison_none$plot, "ggplot")
  
  # Test scenario 3: Complete overlap
  complete_overlap_pathways <- c("ko00010", "ko00020", "ko00030")
  complete_overlap_gsea <- create_controlled_gsea_results(
    significant_pathways = complete_overlap_pathways,
    non_significant_pathways = c("ko00100")
  )
  complete_overlap_daa <- create_controlled_daa_results(
    significant_features = complete_overlap_pathways,
    non_significant_features = c("ko00200")
  )
  
  comparison_complete <- compare_gsea_daa(complete_overlap_gsea, complete_overlap_daa, plot_type = "venn")
  expect_equal(comparison_complete$results$n_overlap, 3)
  expect_equal(comparison_complete$results$n_gsea_only, 0)
  expect_equal(comparison_complete$results$n_daa_only, 0)
  expect_s3_class(comparison_complete$plot, "ggplot")
})

test_that("compare_gsea_daa upset plot creation and fallback behavior", {
  skip_if_not_installed("ggplot2")
  
  # Create test data
  gsea_results <- create_controlled_gsea_results(
    significant_pathways = c("ko00010", "ko00020", "ko00030"),
    non_significant_pathways = c("ko00040")
  )
  daa_results <- create_controlled_daa_results(
    significant_features = c("ko00020", "ko00030", "ko00050"),
    non_significant_features = c("ko00060")
  )
  
  # Test UpSet plot creation (may use fallback if UpSetR not available)
  comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "upset")
  
  expect_type(comparison, "list")
  expect_true("plot" %in% names(comparison))
  expect_true("results" %in% names(comparison))
  expect_equal(comparison$results$n_overlap, 2)  # ko00020, ko00030
  
  # Test fallback behavior when UpSetR is not available
  with_mock(
    `requireNamespace` = function(pkg, quietly = TRUE) {
      if (pkg == "UpSetR") return(FALSE)
      return(TRUE)
    },
    {
      expect_warning(
        comparison_fallback <- compare_gsea_daa(gsea_results, daa_results, plot_type = "upset"),
        "Package 'UpSetR' is required for UpSet plots"
      )
      expect_s3_class(comparison_fallback$plot, "ggplot")
    }
  )
})

test_that("compare_gsea_daa scatter plot handles correlation analysis correctly", {
  skip_if_not_installed("ggplot2")
  
  # Create overlapping data with known correlations
  overlapping_pathways <- c("ko00010", "ko00020", "ko00030", "ko00040", "ko00050")
  
  # Create correlated effect sizes
  nes_values <- c(2.1, 1.5, -1.8, -2.3, 0.8)
  log2fc_values <- c(1.8, 1.2, -1.5, -2.0, 0.6)  # Positively correlated
  
  gsea_results <- create_controlled_gsea_results(
    significant_pathways = overlapping_pathways,
    non_significant_pathways = character(0),
    effect_sizes = nes_values
  )
  
  daa_results <- create_controlled_daa_results(
    significant_features = overlapping_pathways,
    non_significant_features = character(0),
    fold_changes = log2fc_values
  )
  
  # Create scatter plot
  comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "scatter")
  
  expect_s3_class(comparison$plot, "ggplot")
  expect_type(comparison, "list")
  
  # Verify the plot contains the expected data
  plot_data <- ggplot2::ggplot_build(comparison$plot)
  expect_true(length(plot_data$data) > 0)
  
  # Test scatter plot with no overlapping pathways
  non_overlapping_gsea <- create_controlled_gsea_results(
    significant_pathways = c("ko00100", "ko00200"),
    non_significant_pathways = character(0)
  )
  non_overlapping_daa <- create_controlled_daa_results(
    significant_features = c("ko00300", "ko00400"),
    non_significant_features = character(0)
  )
  
  expect_warning(
    comparison_no_overlap <- compare_gsea_daa(non_overlapping_gsea, non_overlapping_daa, plot_type = "scatter"),
    "No overlapping pathways found for scatter plot"
  )
  
  expect_s3_class(comparison_no_overlap$plot, "ggplot")
})

test_that("compare_gsea_daa statistical comparison metrics are accurate", {
  skip_if_not_installed("ggplot2")
  
  # Create precisely controlled data
  gsea_sig <- c("A", "B", "C", "D", "E")  # 5 significant
  gsea_nonsig <- c("F", "G", "H")         # 3 non-significant
  daa_sig <- c("C", "D", "E", "I", "J")   # 5 significant  
  daa_nonsig <- c("K", "L", "M")          # 3 non-significant
  
  gsea_results <- create_controlled_gsea_results(gsea_sig, gsea_nonsig, p_threshold = 0.05)
  daa_results <- create_controlled_daa_results(daa_sig, daa_nonsig, p_threshold = 0.05)
  
  comparison <- compare_gsea_daa(gsea_results, daa_results)
  
  # Expected overlap: C, D, E (3 pathways)
  # Expected GSEA only: A, B (2 pathways) 
  # Expected DAA only: I, J (2 pathways)
  
  expect_equal(comparison$results$n_overlap, 3)
  expect_equal(comparison$results$n_gsea_only, 2) 
  expect_equal(comparison$results$n_daa_only, 2)
  expect_equal(comparison$results$n_gsea_total, 5)
  expect_equal(comparison$results$n_daa_total, 5)
  
  expect_setequal(comparison$results$overlap, c("C", "D", "E"))
  expect_setequal(comparison$results$gsea_only, c("A", "B"))
  expect_setequal(comparison$results$daa_only, c("I", "J"))
  
  # Verify mathematical consistency
  expect_equal(
    comparison$results$n_gsea_total,
    comparison$results$n_overlap + comparison$results$n_gsea_only
  )
  expect_equal(
    comparison$results$n_daa_total, 
    comparison$results$n_overlap + comparison$results$n_daa_only
  )
})

test_that("compare_gsea_daa package dependency handling works correctly", {
  skip_if_not_installed("ggplot2")
  
  gsea_results <- create_controlled_gsea_results(
    significant_pathways = c("ko00010", "ko00020"),
    non_significant_pathways = c("ko00030")
  )
  daa_results <- create_controlled_daa_results(
    significant_features = c("ko00010", "ko00040"),
    non_significant_features = c("ko00050")
  )
  
  # Test ggVennDiagram dependency handling
  with_mock(
    `requireNamespace` = function(pkg, quietly = TRUE) {
      if (pkg == "ggVennDiagram") return(FALSE)
      return(TRUE)
    },
    {
      expect_warning(
        comparison_venn <- compare_gsea_daa(gsea_results, daa_results, plot_type = "venn"),
        "Package 'ggVennDiagram' is required for Venn diagrams. Using a basic plot instead."
      )
      
      expect_s3_class(comparison_venn$plot, "ggplot")
      expect_type(comparison_venn$results, "list")
    }
  )
  
  # Test UpSetR dependency handling  
  with_mock(
    `requireNamespace` = function(pkg, quietly = TRUE) {
      if (pkg == "UpSetR") return(FALSE)
      return(TRUE)
    },
    {
      expect_warning(
        comparison_upset <- compare_gsea_daa(gsea_results, daa_results, plot_type = "upset"),
        "Package 'UpSetR' is required for UpSet plots. Using a basic plot instead."
      )
      
      expect_s3_class(comparison_upset$plot, "ggplot")
      expect_type(comparison_upset$results, "list")
    }
  )
  
  # Test that scatter plots don't require external packages
  comparison_scatter <- compare_gsea_daa(gsea_results, daa_results, plot_type = "scatter")
  expect_s3_class(comparison_scatter$plot, "ggplot")
})

test_that("compare_gsea_daa handles edge cases in comparison calculations", {
  skip_if_not_installed("ggplot2")
  
  # Test with empty significant sets
  gsea_empty_sig <- create_controlled_gsea_results(
    significant_pathways = character(0),
    non_significant_pathways = c("ko00010", "ko00020"),
    p_threshold = 0.05
  )
  daa_empty_sig <- create_controlled_daa_results(
    significant_features = character(0),
    non_significant_features = c("ko00010", "ko00020"),
    p_threshold = 0.05
  )
  
  comparison_empty <- compare_gsea_daa(gsea_empty_sig, daa_empty_sig)
  
  expect_equal(comparison_empty$results$n_overlap, 0)
  expect_equal(comparison_empty$results$n_gsea_only, 0)
  expect_equal(comparison_empty$results$n_daa_only, 0)
  expect_equal(comparison_empty$results$n_gsea_total, 0)
  expect_equal(comparison_empty$results$n_daa_total, 0)
  expect_equal(length(comparison_empty$results$overlap), 0)
  
  # Test with one empty analysis
  gsea_with_sig <- create_controlled_gsea_results(
    significant_pathways = c("ko00010", "ko00020"),
    non_significant_pathways = c("ko00030")
  )
  daa_empty <- create_controlled_daa_results(
    significant_features = character(0),
    non_significant_features = c("ko00010", "ko00020", "ko00030")
  )
  
  comparison_one_empty <- compare_gsea_daa(gsea_with_sig, daa_empty)
  
  expect_equal(comparison_one_empty$results$n_overlap, 0)
  expect_equal(comparison_one_empty$results$n_gsea_only, 2)
  expect_equal(comparison_one_empty$results$n_daa_only, 0)
  expect_equal(comparison_one_empty$results$n_gsea_total, 2)
  expect_equal(comparison_one_empty$results$n_daa_total, 0)
})

test_that("compare_gsea_daa handles missing or malformed columns gracefully", {
  # Test with missing required columns
  incomplete_gsea <- data.frame(
    pathway_name = c("Path1", "Path2"),
    size = c(50, 60),
    stringsAsFactors = FALSE
  )
  
  complete_daa <- create_controlled_daa_results(
    significant_features = c("ko00010", "ko00020"),
    non_significant_features = character(0)
  )
  
  expect_error(
    compare_gsea_daa(incomplete_gsea, complete_daa),
    "GSEA results missing required columns: pathway_id, p.adjust"
  )
  
  # Test with wrong column types
  gsea_wrong_types <- data.frame(
    pathway_id = factor(c("ko00010", "ko00020")),
    pathway_name = c("Path1", "Path2"),
    p.adjust = c("0.05", "0.1"),  # Character instead of numeric
    stringsAsFactors = FALSE
  )
  
  # This should still work as R will convert character to numeric
  result <- compare_gsea_daa(gsea_wrong_types, complete_daa)
  expect_type(result, "list")
  
  # Test with NA values in key columns
  gsea_with_nas <- create_controlled_gsea_results(
    significant_pathways = c("ko00010", "ko00020", "ko00030"),
    non_significant_pathways = character(0)
  )
  gsea_with_nas$pathway_id[2] <- NA
  gsea_with_nas$p.adjust[3] <- NA
  
  daa_with_nas <- create_controlled_daa_results(
    significant_features = c("ko00010", "ko00020", "ko00030"),
    non_significant_features = character(0)
  )
  daa_with_nas$feature[1] <- NA
  
  # Should handle NAs gracefully
  comparison_nas <- compare_gsea_daa(gsea_with_nas, daa_with_nas)
  expect_type(comparison_nas, "list")
  expect_s3_class(comparison_nas$plot, "ggplot")
})

test_that("compare_gsea_daa visualization consistency across plot types", {
  skip_if_not_installed("ggplot2")
  
  # Create consistent test data
  gsea_results <- create_controlled_gsea_results(
    significant_pathways = c("ko00010", "ko00020", "ko00030"),
    non_significant_pathways = c("ko00040", "ko00050")
  )
  daa_results <- create_controlled_daa_results(
    significant_features = c("ko00020", "ko00030", "ko00060"),
    non_significant_features = c("ko00070", "ko00080")
  )
  
  # Create all plot types
  comparison_venn <- compare_gsea_daa(gsea_results, daa_results, plot_type = "venn")
  comparison_upset <- compare_gsea_daa(gsea_results, daa_results, plot_type = "upset")
  comparison_scatter <- compare_gsea_daa(gsea_results, daa_results, plot_type = "scatter")
  
  # All should have same comparison results
  expect_identical(comparison_venn$results, comparison_upset$results)
  expect_identical(comparison_venn$results, comparison_scatter$results)
  
  # All plots should be ggplot objects
  expect_s3_class(comparison_venn$plot, "ggplot")
  expect_s3_class(comparison_upset$plot, "ggplot")
  expect_s3_class(comparison_scatter$plot, "ggplot")
  
  # Test heatmap warning
  expect_warning(
    comparison_heatmap <- compare_gsea_daa(gsea_results, daa_results, plot_type = "heatmap"),
    "Heatmap plot not yet implemented"
  )
  expect_s3_class(comparison_heatmap$plot, "ggplot")
  expect_identical(comparison_heatmap$results, comparison_venn$results)
})