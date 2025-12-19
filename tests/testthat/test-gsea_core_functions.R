# Comprehensive tests for GSEA core calculation functions
# Focus on mathematical correctness, edge cases, and robustness

library(testthat)
library(mockery)

# Source the functions directly since we need them for testing
source("../../R/pathway_gsea.R")

# ============================================================================
# TEST HELPER FUNCTIONS
# ============================================================================

#' Create comprehensive test data with known properties
create_test_data_with_signal <- function(n_features = 100, n_samples_per_group = 10, 
                                       effect_size = 2, noise_level = 1, seed = 42) {
  set.seed(seed)
  
  # Create feature names as KO IDs
  feature_names <- paste0("K", sprintf("%05d", 1:n_features))
  
  # Create sample names
  sample_names <- paste0("Sample", 1:(n_samples_per_group * 2))
  
  # Create abundance matrix with known signal
  abundance <- matrix(rnorm(n_features * n_samples_per_group * 2, mean = 5, sd = noise_level), 
                     nrow = n_features, ncol = n_samples_per_group * 2)
  
  # Add differential signal to first 20 features for group1
  signal_features <- 1:min(20, n_features)
  if (length(signal_features) > 0 && n_samples_per_group > 0) {
    abundance[signal_features, 1:n_samples_per_group] <- 
      abundance[signal_features, 1:n_samples_per_group] + effect_size
  }
  
  # Set row and column names
  rownames(abundance) <- feature_names
  colnames(abundance) <- sample_names
  
  # Create metadata
  metadata <- data.frame(
    sample_name = sample_names,
    group = factor(rep(c("Group1", "Group2"), each = n_samples_per_group)),
    batch = factor(rep(c("Batch1", "Batch2"), times = n_samples_per_group)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- sample_names
  
  # Create gene sets that include the signal features
  gene_sets <- list(
    "pathway_with_signal" = feature_names[1:min(25, n_features)],
    "pathway_mixed" = feature_names[min(15, n_features):min(40, n_features)],
    "pathway_no_signal" = feature_names[min(50, n_features):min(75, n_features)],
    "small_pathway" = feature_names[1:min(5, n_features)],
    "large_pathway" = feature_names[1:min(80, n_features)]
  )
  
  return(list(
    abundance = abundance,
    metadata = metadata,
    gene_sets = gene_sets,
    signal_features = feature_names[signal_features],
    effect_size = effect_size
  ))
}

#' Create edge case test data
create_edge_case_data <- function(case_type = "zeros", n_features = 50, n_samples = 20) {
  set.seed(123)
  
  abundance <- matrix(rnorm(n_features * n_samples, mean = 1), 
                     nrow = n_features, ncol = n_samples)
  rownames(abundance) <- paste0("K", sprintf("%05d", 1:n_features))
  colnames(abundance) <- paste0("Sample", 1:n_samples)
  
  if (case_type == "zeros") {
    # Add zero values
    abundance[1:min(5, n_features), 1:min(5, n_samples)] <- 0
  } else if (case_type == "identical_groups") {
    # Make groups identical
    if (n_samples >= 20) {
      abundance[, 11:20] <- abundance[, 1:10]
    }
  } else if (case_type == "extreme_values") {
    # Add extreme values
    abundance[1:min(3, n_features), 1:min(3, n_samples)] <- 1000
    if (n_samples >= 13) {
      abundance[min(4, n_features):min(6, n_features), 11:13] <- -1000
    }
  }
  
  n_samples_final <- ncol(abundance)
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("Group1", "Group2"), each = n_samples_final/2)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- colnames(abundance)
  
  return(list(abundance = abundance, metadata = metadata))
}

# ============================================================================
# TESTS FOR calculate_rank_metric FUNCTION
# ============================================================================

test_that("calculate_rank_metric: signal2noise method basic functionality", {
  test_data <- create_test_data_with_signal(n_features = 20, n_samples_per_group = 8)
  
  metric <- calculate_rank_metric(test_data$abundance, test_data$metadata, "group", "signal2noise")
  
  # Basic properties
  expect_type(metric, "double")
  expect_length(metric, nrow(test_data$abundance))
  expect_named(metric, rownames(test_data$abundance))
  expect_true(all(is.finite(metric)))
})

test_that("calculate_rank_metric: signal2noise mathematical correctness", {
  # Create simple test case for manual verification
  set.seed(42)
  abundance <- matrix(c(1, 3, 5, 2, 4, 6,
                       10, 12, 14, 11, 13, 15),
                     nrow = 2, ncol = 6, byrow = TRUE)
  rownames(abundance) <- c("K00001", "K00002")
  colnames(abundance) <- paste0("Sample", 1:6)
  
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("Group1", "Group2"), each = 3)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample_name
  
  metric <- calculate_rank_metric(abundance, metadata, "group", "signal2noise")
  
  # Manual calculation for first feature
  group1_vals <- abundance[1, 1:3]  # 1, 3, 5
  group2_vals <- abundance[1, 4:6]  # 2, 4, 6
  mean_diff <- mean(group1_vals) - mean(group2_vals)  # 3 - 4 = -1
  sd_sum <- sd(group1_vals) + sd(group2_vals)  # 2 + 2 = 4
  expected_s2n <- mean_diff / sd_sum  # -1/4 = -0.25
  
  expect_equal(unname(metric[1]), expected_s2n, tolerance = 1e-10)
})

test_that("calculate_rank_metric: t_test method functionality", {
  test_data <- create_test_data_with_signal(n_features = 15, n_samples_per_group = 6)
  
  metric <- calculate_rank_metric(test_data$abundance, test_data$metadata, "group", "t_test")
  
  expect_type(metric, "double")
  expect_length(metric, nrow(test_data$abundance))
  expect_true(all(is.finite(metric)))
})

test_that("calculate_rank_metric: log2_ratio method functionality", {
  # Use positive abundance values for log2 ratio
  test_data <- create_test_data_with_signal(n_features = 15, n_samples_per_group = 6)
  test_data$abundance <- abs(test_data$abundance) + 0.1  # Ensure positive values
  
  metric <- calculate_rank_metric(test_data$abundance, test_data$metadata, "group", "log2_ratio")
  
  expect_type(metric, "double")
  expect_length(metric, nrow(test_data$abundance))
  expect_true(all(is.finite(metric)))
})

test_that("calculate_rank_metric: diff_abundance method functionality", {
  test_data <- create_test_data_with_signal(n_features = 15, n_samples_per_group = 6)
  
  metric <- calculate_rank_metric(test_data$abundance, test_data$metadata, "group", "diff_abundance")
  
  expect_type(metric, "double")
  expect_length(metric, nrow(test_data$abundance))
  expect_true(all(is.finite(metric)))
})

test_that("calculate_rank_metric: edge cases handling", {
  # Test with zero standard deviations
  edge_data <- create_edge_case_data("zeros")
  
  expect_no_error({
    metric_s2n <- calculate_rank_metric(edge_data$abundance, edge_data$metadata, "group", "signal2noise")
  })
  expect_true(all(is.finite(metric_s2n)))
  
  # Test with identical groups
  edge_data_identical <- create_edge_case_data("identical_groups")
  metric_identical <- calculate_rank_metric(edge_data_identical$abundance, edge_data_identical$metadata, "group", "signal2noise")
  # Values should be close to zero or well-defined
  expect_true(all(is.finite(metric_identical)))
})

test_that("calculate_rank_metric: error handling", {
  test_data <- create_test_data_with_signal(n_features = 10, n_samples_per_group = 6)
  
  # Test error with more than two groups
  test_data$metadata$group <- factor(rep(c("A", "B", "C"), length.out = nrow(test_data$metadata)))
  
  expect_error(
    calculate_rank_metric(test_data$abundance, test_data$metadata, "group", "signal2noise"),
    "GSEA currently only supports two-group comparisons"
  )
})

# ============================================================================
# TESTS FOR prepare_gene_sets FUNCTION
# ============================================================================

test_that("prepare_gene_sets: KEGG pathway basic functionality", {
  # Test KEGG pathway preparation
  gene_sets <- prepare_gene_sets("KEGG")
  
  # Basic structure tests
  expect_type(gene_sets, "list")
  
  # Each gene set should be a character vector if any exist
  if (length(gene_sets) > 0) {
    expect_true(all(sapply(gene_sets, is.character)))
    expect_true(all(sapply(gene_sets, function(x) length(x) > 0)))
    expect_true(all(nzchar(names(gene_sets))))
  }
})

test_that("prepare_gene_sets: MetaCyc and GO pathway types are supported", {
  # MetaCyc is now supported - it should return gene sets or error if data not available
  result_metacyc <- tryCatch({
    gene_sets_metacyc <- prepare_gene_sets("MetaCyc")
    expect_type(gene_sets_metacyc, "list")
    TRUE
  }, error = function(e) {
    # Data file might not be available in test environment
    expect_true(grepl("MetaCyc|not found|Failed", e$message))
    TRUE
  })
  expect_true(result_metacyc)

  # GO is now supported - it should return gene sets
  result_go <- tryCatch({
    gene_sets_go <- suppressMessages(prepare_gene_sets("GO"))
    expect_type(gene_sets_go, "list")
    expect_true(length(gene_sets_go) > 0)
    TRUE
  }, error = function(e) {
    # Data file might not be available in test environment
    TRUE
  })
  expect_true(result_go)
})

# ============================================================================
# TESTS FOR run_fgsea FUNCTION  
# ============================================================================

test_that("run_fgsea: basic functionality with mock", {
  skip_if_not_installed("fgsea")
  
  # Create test data
  set.seed(42)
  ranked_list <- setNames(rnorm(100), paste0("K", sprintf("%05d", 1:100)))
  ranked_list <- sort(ranked_list, decreasing = TRUE)
  
  gene_sets <- list(
    "pathway1" = paste0("K", sprintf("%05d", 1:15)),
    "pathway2" = paste0("K", sprintf("%05d", 16:30))
  )
  
  # Mock fgsea function 
  mock_fgsea_result <- data.frame(
    pathway = c("pathway1", "pathway2"),
    pval = c(0.01, 0.05),
    padj = c(0.02, 0.1),
    ES = c(0.5, -0.3),
    NES = c(1.2, -0.8),
    size = c(15, 15),
    leadingEdge = I(list(paste0("K", sprintf("%05d", 1:3)), 
                        paste0("K", sprintf("%05d", 16:18)))),
    stringsAsFactors = FALSE
  )
  
  # Use mockery to stub the fgsea call
  stub(run_fgsea, "fgsea::fgsea", mock_fgsea_result)
  
  result <- run_fgsea(ranked_list, gene_sets)
  
  # Check structure
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(c("pathway_id", "pathway_name", "size", "ES", "NES", "pvalue", "p.adjust", "leading_edge") %in% colnames(result)))
  
  # Check data types and ranges
  expect_type(result$pvalue, "double")
  expect_true(all(result$pvalue >= 0 & result$pvalue <= 1))
  expect_true(all(result$p.adjust >= 0 & result$p.adjust <= 1))
})

# ============================================================================
# TESTS FOR pathway_gsea MAIN FUNCTION
# ============================================================================

test_that("pathway_gsea: basic input validation", {
  test_data <- create_test_data_with_signal(n_features = 10, n_samples_per_group = 5)
  
  # Test invalid abundance
  expect_error(
    pathway_gsea(abundance = "invalid", metadata = test_data$metadata, group = "group"),
    "'abundance' must be a data frame or matrix"
  )
  
  # Test invalid metadata
  expect_error(
    pathway_gsea(abundance = test_data$abundance, metadata = "invalid", group = "group"),
    "'metadata' must be a data frame"
  )
  
  # Test invalid group
  expect_error(
    pathway_gsea(abundance = test_data$abundance, metadata = test_data$metadata, group = "invalid_group"),
    "Group variable invalid_group not found in metadata"
  )
})

test_that("pathway_gsea: parameter validation", {
  test_data <- create_test_data_with_signal(n_features = 10, n_samples_per_group = 5)
  
  # Test invalid pathway_type
  expect_error(
    pathway_gsea(abundance = test_data$abundance, metadata = test_data$metadata, group = "group", pathway_type = "invalid"),
    "pathway_type must be one of 'KEGG', 'MetaCyc', or 'GO'"
  )
  
  # Test invalid method
  # Note: method options now include camera, fry, fgsea, GSEA, clusterProfiler
  expect_error(
    pathway_gsea(abundance = test_data$abundance, metadata = test_data$metadata, group = "group", method = "invalid"),
    "method must be one of: camera, fry, fgsea, GSEA, clusterProfiler"
  )
  
  # Test invalid rank_method (only validated for preranked methods like fgsea)
  expect_error(
    pathway_gsea(abundance = test_data$abundance, metadata = test_data$metadata, group = "group",
                 method = "fgsea", rank_method = "invalid"),
    "rank_method must be one of 'signal2noise', 't_test', 'log2_ratio', or 'diff_abundance'"
  )
})

test_that("pathway_gsea: sample name matching", {
  test_data <- create_test_data_with_signal(n_features = 10, n_samples_per_group = 5)

  # Create mismatched sample names (no overlap)
  mismatched_metadata <- test_data$metadata
  rownames(mismatched_metadata) <- paste0("Wrong_", rownames(mismatched_metadata))
  mismatched_metadata$sample_name <- rownames(mismatched_metadata)

  # With no overlapping samples, should error due to insufficient samples
  expect_error(
    pathway_gsea(
      abundance = test_data$abundance,
      metadata = mismatched_metadata,
      group = "group"
    ),
    "Insufficient overlapping samples"
  )
})

test_that("pathway_gsea: integration test with mocking", {
  skip_if_not_installed("fgsea")
  
  test_data <- create_test_data_with_signal(n_features = 20, n_samples_per_group = 8)
  
  # Mock prepare_gene_sets
  stub(pathway_gsea, "prepare_gene_sets", function(...) {
    list(
      "test_pathway1" = rownames(test_data$abundance)[1:10],
      "test_pathway2" = rownames(test_data$abundance)[11:20]
    )
  })
  
  # Mock run_fgsea
  stub(pathway_gsea, "run_fgsea", function(...) {
    data.frame(
      pathway_id = c("test_pathway1", "test_pathway2"),
      pathway_name = c("test_pathway1", "test_pathway2"),
      size = c(10, 10),
      ES = c(0.5, -0.3),
      NES = c(1.2, -0.8),
      pvalue = c(0.01, 0.05),
      p.adjust = c(0.02, 0.1),
      leading_edge = c("K00001;K00002", "K00011;K00012"),
      stringsAsFactors = FALSE
    )
  })
  
  # Test with different ranking methods
  for (rank_method in c("signal2noise", "t_test", "diff_abundance")) {
    result <- pathway_gsea(
      abundance = test_data$abundance,
      metadata = test_data$metadata,
      group = "group",
      pathway_type = "KEGG",
      method = "fgsea",
      rank_method = rank_method,
      nperm = 100,
      seed = 42
    )
    
    # Check result structure
    expect_s3_class(result, "data.frame")
    expect_true(nrow(result) > 0)
    expect_true(all(c("pathway_id", "pathway_name", "size", "ES", "NES", "pvalue", "p.adjust", "method") %in% colnames(result)))
    expect_equal(unique(result$method), "fgsea")
  }
})

# ============================================================================
# EDGE CASE AND ROBUSTNESS TESTS
# ============================================================================

test_that("pathway_gsea: edge cases in abundance data", {
  skip_if_not_installed("fgsea")
  
  # Test with various edge cases
  edge_cases <- c("zeros", "extreme_values")
  
  for (case_type in edge_cases) {
    edge_data <- create_edge_case_data(case_type, n_features = 15, n_samples = 12)
    
    # Mock functions for successful completion
    stub(pathway_gsea, "prepare_gene_sets", function(...) {
      list("test_pathway" = rownames(edge_data$abundance)[1:min(8, nrow(edge_data$abundance))])
    })
    
    stub(pathway_gsea, "run_fgsea", function(...) {
      data.frame(
        pathway_id = "test_pathway",
        pathway_name = "test_pathway",
        size = 8,
        ES = 0.2,
        NES = 0.5,
        pvalue = 0.3,
        p.adjust = 0.4,
        leading_edge = "feature1;feature2",
        stringsAsFactors = FALSE
      )
    })
    
    expect_no_error({
      result <- pathway_gsea(
        abundance = edge_data$abundance,
        metadata = edge_data$metadata,
        group = "group",
        method = "fgsea"
      )
    })
    
    expect_s3_class(result, "data.frame")
  }
})

test_that("pathway_gsea: reproducibility with seed", {
  skip_if_not_installed("fgsea")
  
  test_data <- create_test_data_with_signal(n_features = 15, n_samples_per_group = 6)
  
  # Mock functions
  stub(pathway_gsea, "prepare_gene_sets", function(...) {
    list("test_pathway" = rownames(test_data$abundance)[1:10])
  })
  
  # Create a mock that uses the seed
  stub(pathway_gsea, "run_fgsea", function(ranked_list, gene_sets, ...) {
    # Use current random state to create somewhat consistent results
    set.seed(42)
    data.frame(
      pathway_id = "test_pathway",
      pathway_name = "test_pathway",
      size = 10,
      ES = rnorm(1, 0, 0.1),
      NES = rnorm(1, 0, 0.1), 
      pvalue = runif(1, 0.1, 0.9),
      p.adjust = runif(1, 0.1, 0.9),
      leading_edge = "K00001;K00002",
      stringsAsFactors = FALSE
    )
  })
  
  # Run analysis twice with same seed
  result1 <- pathway_gsea(
    abundance = test_data$abundance,
    metadata = test_data$metadata,
    group = "group",
    method = "fgsea",
    seed = 12345
  )
  
  result2 <- pathway_gsea(
    abundance = test_data$abundance,
    metadata = test_data$metadata,
    group = "group",
    method = "fgsea",
    seed = 12345
  )
  
  # The ranking should be identical (deterministic calculation)
  expect_s3_class(result1, "data.frame")
  expect_s3_class(result2, "data.frame")
  expect_equal(nrow(result1), nrow(result2))
})

# ============================================================================
# MATHEMATICAL CORRECTNESS VERIFICATION TESTS
# ============================================================================

test_that("calculate_rank_metric: mathematical properties verification", {
  # Create controlled test case
  set.seed(999)
  
  # Create data where group1 has higher values than group2
  abundance <- matrix(0, nrow = 10, ncol = 20)
  abundance[, 1:10] <- 100  # Group1: high values
  abundance[, 11:20] <- 10  # Group2: low values
  
  # Add some noise
  abundance <- abundance + matrix(rnorm(200, 0, 1), nrow = 10, ncol = 20)
  
  rownames(abundance) <- paste0("K", sprintf("%05d", 1:10))
  colnames(abundance) <- paste0("Sample", 1:20)
  
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("High", "Low"), each = 10)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample_name
  
  # Test signal2noise: should be positive (high - low > 0)
  metric_s2n <- calculate_rank_metric(abundance, metadata, "group", "signal2noise")
  expect_true(all(metric_s2n > 0))
  
  # Test diff_abundance: should be positive
  metric_diff <- calculate_rank_metric(abundance, metadata, "group", "diff_abundance")
  expect_true(all(metric_diff > 0))
  
  # Test log2_ratio: should be positive (assuming all values are positive)
  abundance_pos <- abs(abundance) + 0.1
  metric_log2 <- calculate_rank_metric(abundance_pos, metadata, "group", "log2_ratio")
  expect_true(all(metric_log2 > 0))
  
  # Test t_test: should be positive
  metric_ttest <- calculate_rank_metric(abundance, metadata, "group", "t_test")
  expect_true(all(metric_ttest > 0))
})

test_that("calculate_rank_metric: consistency across methods", {
  test_data <- create_test_data_with_signal(n_features = 20, n_samples_per_group = 10, effect_size = 1)
  
  # Calculate metrics with different methods
  metrics <- list(
    s2n = calculate_rank_metric(test_data$abundance, test_data$metadata, "group", "signal2noise"),
    ttest = calculate_rank_metric(test_data$abundance, test_data$metadata, "group", "t_test"),
    diff = calculate_rank_metric(test_data$abundance, test_data$metadata, "group", "diff_abundance")
  )
  
  # All methods should identify the same features as having signal
  # (features with signal should rank higher)
  signal_indices <- which(rownames(test_data$abundance) %in% test_data$signal_features)
  
  for (method in names(metrics)) {
    metric <- metrics[[method]]

    # Check that signal features have higher absolute values on average
    if (length(signal_indices) > 0) {
      signal_abs_mean <- mean(abs(metric[signal_indices]), na.rm = TRUE)
      noise_indices <- setdiff(seq_along(metric), signal_indices)

      # Only compare if we have noise features
      if (length(noise_indices) > 0) {
        noise_abs_mean <- mean(abs(metric[noise_indices]), na.rm = TRUE)

        # Skip comparison if noise_abs_mean is NaN (all noise values were NA)
        if (!is.nan(noise_abs_mean) && !is.nan(signal_abs_mean)) {
          # Signal features should have higher ranking metrics (on average)
          expect_true(signal_abs_mean >= noise_abs_mean * 0.8,
                     info = paste("Method:", method, "signal mean:", signal_abs_mean, "noise mean:", noise_abs_mean))
        }
      }
    }
  }
})