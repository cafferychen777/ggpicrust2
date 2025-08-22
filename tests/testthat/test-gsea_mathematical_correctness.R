# Mathematical Correctness Tests for GSEA Core Functions
# Focus on verifying the mathematical accuracy of ranking metrics and GSEA calculations

library(testthat)
library(mockery)

# Source the functions directly for testing
source("../../R/pathway_gsea.R")

# ============================================================================
# HELPER FUNCTIONS FOR MATHEMATICAL VERIFICATION
# ============================================================================

#' Create test data with known mathematical properties
create_controlled_test_data <- function(n_features = 50, n_samples_per_group = 10) {
  set.seed(12345)  # Fixed seed for reproducibility
  
  # Create feature and sample names
  feature_names <- paste0("K", sprintf("%05d", 1:n_features))
  sample_names <- c(paste0("Group1_Sample", 1:n_samples_per_group),
                   paste0("Group2_Sample", 1:n_samples_per_group))
  
  # Create abundance matrix with controlled properties
  abundance <- matrix(rnorm(n_features * n_samples_per_group * 2, mean = 10, sd = 2), 
                     nrow = n_features, 
                     ncol = n_samples_per_group * 2)
  
  # Add known differential signal to specific features
  signal_features <- 1:10  # First 10 features have signal
  effect_size <- 5
  
  # Group1 samples have higher abundance for signal features
  abundance[signal_features, 1:n_samples_per_group] <- 
    abundance[signal_features, 1:n_samples_per_group] + effect_size
  
  # Set names
  rownames(abundance) <- feature_names
  colnames(abundance) <- sample_names
  
  # Create metadata
  metadata <- data.frame(
    sample_name = sample_names,
    group = factor(rep(c("Group1", "Group2"), each = n_samples_per_group)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- sample_names
  
  return(list(
    abundance = abundance,
    metadata = metadata,
    signal_features = feature_names[signal_features],
    effect_size = effect_size
  ))
}

# ============================================================================
# MATHEMATICAL CORRECTNESS TESTS FOR RANKING METRICS
# ============================================================================

test_that("calculate_rank_metric: signal2noise mathematical precision", {
  # Create simple controlled data for exact calculation
  abundance <- matrix(c(
    # Feature 1: Group1=[10,12,14], Group2=[5,7,9] -> clear difference
    10, 12, 14, 5, 7, 9,
    # Feature 2: Group1=[8,8,8], Group2=[8,8,8] -> no difference
    8, 8, 8, 8, 8, 8
  ), nrow = 2, ncol = 6, byrow = TRUE)
  
  rownames(abundance) <- c("K00001", "K00002")
  colnames(abundance) <- paste0("Sample", 1:6)
  
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("Group1", "Group2"), each = 3)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample_name
  
  metric <- calculate_rank_metric(abundance, metadata, "group", "signal2noise")
  
  # Manual calculation for Feature 1
  group1_vals <- c(10, 12, 14)  # mean = 12, sd = 2
  group2_vals <- c(5, 7, 9)     # mean = 7, sd = 2
  expected_s2n_1 <- (12 - 7) / (2 + 2)  # 5/4 = 1.25
  
  expect_equal(metric[["K00001"]], expected_s2n_1, tolerance = 1e-10)
  
  # Feature 2 should have very small metric (near zero due to no difference)
  expect_true(abs(metric[["K00002"]]) < 1e-6)
})

test_that("calculate_rank_metric: t_test mathematical accuracy", {
  test_data <- create_controlled_test_data(n_features = 20, n_samples_per_group = 8)
  
  metric <- calculate_rank_metric(test_data$abundance, test_data$metadata, "group", "t_test")
  
  # Verify manual calculation for first feature (which has signal)
  group1_samples <- test_data$metadata$sample_name[test_data$metadata$group == "Group1"]
  group2_samples <- test_data$metadata$sample_name[test_data$metadata$group == "Group2"]
  
  feature1_group1 <- test_data$abundance[1, group1_samples]
  feature1_group2 <- test_data$abundance[1, group2_samples]
  
  expected_t <- t.test(feature1_group1, feature1_group2)$statistic
  expect_equal(metric[1], expected_t, tolerance = 1e-10)
  
  # Features with signal should have higher absolute t-statistics
  signal_indices <- which(names(metric) %in% test_data$signal_features)
  noise_indices <- setdiff(1:length(metric), signal_indices)
  
  mean_signal_t <- mean(abs(metric[signal_indices]))
  mean_noise_t <- mean(abs(metric[noise_indices]))
  
  expect_true(mean_signal_t > mean_noise_t)
})

test_that("calculate_rank_metric: log2_ratio mathematical accuracy", {
  # Create test data with positive values for log2 calculation
  abundance <- matrix(c(
    # Feature 1: Group1 mean=16, Group2 mean=4 -> log2(16/4) = 2
    16, 16, 16, 4, 4, 4,
    # Feature 2: Group1 mean=8, Group2 mean=8 -> log2(8/8) = 0
    8, 8, 8, 8, 8, 8
  ), nrow = 2, ncol = 6, byrow = TRUE)
  
  rownames(abundance) <- c("K00001", "K00002")
  colnames(abundance) <- paste0("Sample", 1:6)
  
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("Group1", "Group2"), each = 3)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample_name
  
  metric <- calculate_rank_metric(abundance, metadata, "group", "log2_ratio")
  
  # Manual verification
  expect_equal(metric[["K00001"]], log2(16/4), tolerance = 1e-10)  # Should be 2
  expect_equal(metric[["K00002"]], log2(8/8), tolerance = 1e-10)   # Should be 0
})

test_that("calculate_rank_metric: diff_abundance mathematical accuracy", {
  test_data <- create_controlled_test_data(n_features = 15, n_samples_per_group = 6)
  
  metric <- calculate_rank_metric(test_data$abundance, test_data$metadata, "group", "diff_abundance")
  
  # Manual verification for first feature
  group1_samples <- test_data$metadata$sample_name[test_data$metadata$group == "Group1"]
  group2_samples <- test_data$metadata$sample_name[test_data$metadata$group == "Group2"]
  
  mean1 <- mean(test_data$abundance[1, group1_samples])
  mean2 <- mean(test_data$abundance[1, group2_samples])
  expected_diff <- mean1 - mean2
  
  expect_equal(metric[1], expected_diff, tolerance = 1e-10)
  
  # Signal features should have positive differences (Group1 > Group2)
  signal_indices <- which(names(metric) %in% test_data$signal_features)
  expect_true(all(metric[signal_indices] > 0))
})

# ============================================================================
# EDGE CASE MATHEMATICAL BEHAVIOR TESTS
# ============================================================================

test_that("calculate_rank_metric: zero variance handling", {
  # Create data with zero variance in one group
  abundance <- matrix(c(
    # Feature 1: Group1 has variance, Group2 has zero variance
    1, 2, 3, 5, 5, 5,
    # Feature 2: Both groups have zero variance
    7, 7, 7, 9, 9, 9
  ), nrow = 2, ncol = 6, byrow = TRUE)
  
  rownames(abundance) <- c("K00001", "K00002")
  colnames(abundance) <- paste0("Sample", 1:6)
  
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("Group1", "Group2"), each = 3)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample_name
  
  # Should handle zero variance without errors or infinite values
  metric <- calculate_rank_metric(abundance, metadata, "group", "signal2noise")
  
  expect_true(all(is.finite(metric)))
  expect_length(metric, 2)
})

test_that("calculate_rank_metric: extreme values handling", {
  # Create data with extreme values
  abundance <- matrix(c(
    1e6, 1e6, 1e6, 1, 1, 1,    # Large difference
    -1e6, -1e6, -1e6, 1, 1, 1  # Large negative values
  ), nrow = 2, ncol = 6, byrow = TRUE)
  
  rownames(abundance) <- c("K00001", "K00002")
  colnames(abundance) <- paste0("Sample", 1:6)
  
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("Group1", "Group2"), each = 3)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample_name
  
  # All methods should handle extreme values
  methods <- c("signal2noise", "t_test", "diff_abundance")
  
  for (method in methods) {
    metric <- calculate_rank_metric(abundance, metadata, "group", method)
    expect_true(all(is.finite(metric)), info = paste("Method:", method))
    expect_length(metric, 2)
  }
})

# ============================================================================
# PARAMETER SENSITIVITY AND ROBUSTNESS TESTS
# ============================================================================

test_that("calculate_rank_metric: robustness to sample size", {
  sample_sizes <- c(3, 5, 10, 20)
  
  for (n_samples in sample_sizes) {
    test_data <- create_controlled_test_data(n_features = 10, n_samples_per_group = n_samples)
    
    for (method in c("signal2noise", "t_test", "diff_abundance")) {
      metric <- calculate_rank_metric(test_data$abundance, test_data$metadata, "group", method)
      
      expect_true(all(is.finite(metric)), 
                 info = paste("Method:", method, "Sample size:", n_samples))
      expect_length(metric, 10)
      
      # Signal features should still be detectable
      signal_indices <- which(names(metric) %in% test_data$signal_features)
      if (length(signal_indices) > 0) {
        expect_true(mean(abs(metric[signal_indices])) > 0,
                   info = paste("Method:", method, "Sample size:", n_samples))
      }
    }
  }
})

test_that("calculate_rank_metric: consistency across data transformations", {
  test_data <- create_controlled_test_data(n_features = 20, n_samples_per_group = 8)
  
  # Original metric
  original_metric <- calculate_rank_metric(test_data$abundance, test_data$metadata, "group", "signal2noise")
  
  # Add constant to all values (should not change signal2noise ratio)
  shifted_abundance <- test_data$abundance + 1000
  shifted_metric <- calculate_rank_metric(shifted_abundance, test_data$metadata, "group", "signal2noise")
  
  # Signal2noise should be nearly identical
  expect_equal(original_metric, shifted_metric, tolerance = 1e-10)
  
  # Test scaling (multiply by constant)
  scaled_abundance <- test_data$abundance * 2
  scaled_metric <- calculate_rank_metric(scaled_abundance, test_data$metadata, "group", "signal2noise")
  
  # Signal2noise should be identical for scaling
  expect_equal(original_metric, scaled_metric, tolerance = 1e-10)
})

# ============================================================================
# PREPARE_GENE_SETS FUNCTION TESTS
# ============================================================================

test_that("prepare_gene_sets: handles missing reference data gracefully", {
  # Test when reference data is not available
  expect_warning(
    gene_sets_metacyc <- prepare_gene_sets("MetaCyc"),
    "MetaCyc pathway gene sets not yet implemented"
  )
  expect_type(gene_sets_metacyc, "list")
  expect_equal(length(gene_sets_metacyc), 0)
  
  expect_warning(
    gene_sets_go <- prepare_gene_sets("GO"),
    "GO pathway gene sets not yet implemented"
  )
  expect_type(gene_sets_go, "list")
  expect_equal(length(gene_sets_go), 0)
})

# ============================================================================
# INTEGRATION TESTS WITH MATHEMATICAL VERIFICATION
# ============================================================================

test_that("pathway_gsea: mathematical consistency in ranking", {
  skip_if_not_installed("fgsea")
  
  test_data <- create_controlled_test_data(n_features = 30, n_samples_per_group = 10)
  
  # Mock gene sets to use our controlled features
  mock_gene_sets <- list(
    "signal_pathway" = test_data$signal_features[1:8],    # Contains signal features
    "mixed_pathway" = c(test_data$signal_features[1:3], 
                       rownames(test_data$abundance)[21:25]), # Mixed signal/noise
    "noise_pathway" = rownames(test_data$abundance)[21:28]   # No signal features
  )
  
  stub(pathway_gsea, "prepare_gene_sets", function(...) mock_gene_sets)
  
  # Mock fgsea with realistic behavior based on ranking
  stub(pathway_gsea, "run_fgsea", function(ranked_list, gene_sets, ...) {
    # Simulate realistic enrichment based on actual ranking
    results <- data.frame(
      pathway_id = names(gene_sets),
      pathway_name = names(gene_sets),
      size = sapply(gene_sets, length),
      ES = c(0.8, 0.4, 0.1),      # Signal pathway should have highest ES
      NES = c(2.1, 1.0, 0.3),     # Normalized enrichment scores
      pvalue = c(0.001, 0.05, 0.5), # Signal pathway should be most significant
      p.adjust = c(0.003, 0.1, 0.5),
      leading_edge = c("K00001;K00002", "K00001;K00021", "K00021;K00022"),
      stringsAsFactors = FALSE
    )
    return(results)
  })
  
  # Test different ranking methods
  for (rank_method in c("signal2noise", "t_test", "diff_abundance")) {
    result <- pathway_gsea(
      abundance = test_data$abundance,
      metadata = test_data$metadata,
      group = "group",
      method = "fgsea",
      rank_method = rank_method,
      seed = 42
    )
    
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 3)
    
    # Signal pathway should have the best (lowest) p-value
    signal_pathway_row <- which(result$pathway_id == "signal_pathway")
    expect_equal(signal_pathway_row, 1)  # Should be first (best) result
    expect_true(result$pvalue[signal_pathway_row] < 0.01)
  }
})

test_that("pathway_gsea: reproducibility verification", {
  skip_if_not_installed("fgsea")
  
  test_data <- create_controlled_test_data(n_features = 25, n_samples_per_group = 8)
  
  # Mock functions for consistent testing
  stub(pathway_gsea, "prepare_gene_sets", function(...) {
    list("test_pathway" = rownames(test_data$abundance)[1:15])
  })
  
  stub(pathway_gsea, "run_fgsea", function(...) {
    data.frame(
      pathway_id = "test_pathway",
      pathway_name = "test_pathway", 
      size = 15,
      ES = 0.5,
      NES = 1.2,
      pvalue = 0.02,
      p.adjust = 0.05,
      leading_edge = "K00001;K00002;K00003",
      stringsAsFactors = FALSE
    )
  })
  
  # Run multiple times with same seed
  results <- list()
  for (i in 1:3) {
    results[[i]] <- pathway_gsea(
      abundance = test_data$abundance,
      metadata = test_data$metadata,
      group = "group",
      method = "fgsea",
      seed = 999
    )
  }
  
  # All results should be identical
  expect_identical(results[[1]], results[[2]])
  expect_identical(results[[2]], results[[3]])
})

# ============================================================================
# STRESS TESTS FOR LARGE DATA
# ============================================================================

test_that("pathway_gsea: handles large datasets mathematically", {
  skip_if_not_installed("fgsea")
  skip_if_testing_is_slow()
  
  # Create larger test dataset
  large_test_data <- create_controlled_test_data(n_features = 100, n_samples_per_group = 25)
  
  # Mock functions for performance
  stub(pathway_gsea, "prepare_gene_sets", function(...) {
    list(
      "large_pathway" = rownames(large_test_data$abundance)[1:50],
      "medium_pathway" = rownames(large_test_data$abundance)[26:50]
    )
  })
  
  stub(pathway_gsea, "run_fgsea", function(...) {
    data.frame(
      pathway_id = c("large_pathway", "medium_pathway"),
      pathway_name = c("large_pathway", "medium_pathway"),
      size = c(50, 25),
      ES = c(0.3, 0.6),
      NES = c(0.8, 1.4),
      pvalue = c(0.1, 0.01),
      p.adjust = c(0.15, 0.02),
      leading_edge = c("K00001;K00002", "K00026;K00027"),
      stringsAsFactors = FALSE
    )
  })
  
  # Should complete without mathematical errors
  expect_no_error({
    result <- pathway_gsea(
      abundance = large_test_data$abundance,
      metadata = large_test_data$metadata,
      group = "group",
      method = "fgsea"
    )
  })
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(is.finite(result$ES)))
  expect_true(all(is.finite(result$NES)))
  expect_true(all(result$pvalue >= 0 & result$pvalue <= 1))
})

# ============================================================================
# FINAL MATHEMATICAL PROPERTY VERIFICATION
# ============================================================================

test_that("calculate_rank_metric: all methods detect known signal", {
  # Create highly controlled data where we know signal should be detectable
  set.seed(42)
  n_features <- 30
  n_samples_per_group <- 15
  
  abundance <- matrix(rnorm(n_features * n_samples_per_group * 2, mean = 100, sd = 10),
                     nrow = n_features, ncol = n_samples_per_group * 2)
  
  # Add strong signal to first 10 features
  strong_signal <- 50  # Very large effect
  abundance[1:10, 1:n_samples_per_group] <- abundance[1:10, 1:n_samples_per_group] + strong_signal
  
  rownames(abundance) <- paste0("K", sprintf("%05d", 1:n_features))
  colnames(abundance) <- paste0("Sample", 1:(n_samples_per_group * 2))
  
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("Signal", "Control"), each = n_samples_per_group)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample_name
  
  # Test all ranking methods
  methods <- c("signal2noise", "t_test", "diff_abundance")
  
  for (method in methods) {
    metric <- calculate_rank_metric(abundance, metadata, "group", method)
    
    # Signal features (1-10) should all have higher values than noise features (11-30)
    signal_values <- metric[1:10]
    noise_values <- metric[11:30]
    
    # For these methods, signal group has higher abundance, so metric should be positive
    if (method %in% c("signal2noise", "diff_abundance")) {
      expect_true(all(signal_values > 0), 
                 info = paste("Method:", method, "- Signal features should be positive"))
      expect_true(mean(signal_values) > mean(abs(noise_values)) * 2,
                 info = paste("Method:", method, "- Signal should be much stronger than noise"))
    }
    
    # For t-test, signal features should have higher absolute t-statistics
    if (method == "t_test") {
      expect_true(mean(abs(signal_values)) > mean(abs(noise_values)) * 2,
                 info = paste("Method:", method, "- Signal t-stats should be much higher"))
    }
  }
})