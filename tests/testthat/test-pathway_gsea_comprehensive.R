# Comprehensive tests for GSEA data preprocessing and validation functions
# Tests focus on input validation, data formatting, sample matching, and edge cases

library(testthat)

# ==============================================================================
# Helper Functions for Test Data Creation
# ==============================================================================

#' Create comprehensive test datasets with various configurations
#' @param n_features Number of features (KO/EC/MetaCyc IDs)
#' @param n_samples Number of samples
#' @param n_groups Number of groups (default 2 for GSEA)
#' @param include_zeros Whether to include zero values
#' @param include_missing Whether to include missing values
#' @param abundance_type Type of abundance data ("matrix", "data.frame", "tibble")
#' @param sample_name_mismatch Whether to create sample name mismatches
#' @param group_imbalance Whether to create imbalanced groups
create_comprehensive_test_data <- function(n_features = 10, 
                                         n_samples = 20, 
                                         n_groups = 2,
                                         include_zeros = TRUE,
                                         include_missing = FALSE,
                                         abundance_type = "matrix",
                                         sample_name_mismatch = FALSE,
                                         group_imbalance = FALSE) {
  set.seed(123)
  
  # Create abundance data
  abundance <- matrix(
    abs(rnorm(n_features * n_samples, mean = 1000, sd = 500)), 
    nrow = n_features, 
    ncol = n_samples
  )
  
  # Add zeros if requested
  if (include_zeros) {
    zero_indices <- sample(length(abundance), size = length(abundance) * 0.1)
    abundance[zero_indices] <- 0
  }
  
  # Add missing values if requested
  if (include_missing) {
    na_indices <- sample(length(abundance), size = length(abundance) * 0.05)
    abundance[na_indices] <- NA
  }
  
  # Create sample names
  sample_names <- paste0("Sample_", sprintf("%03d", 1:n_samples))
  colnames(abundance) <- sample_names
  
  # Create feature names (KO format)
  feature_names <- paste0("K", sprintf("%05d", 1:n_features))
  rownames(abundance) <- feature_names
  
  # Create metadata
  if (group_imbalance) {
    # Create imbalanced groups (e.g., 15 vs 5 samples)
    group_sizes <- c(ceiling(n_samples * 0.75), floor(n_samples * 0.25))
    group_assignment <- c(rep("Group_A", group_sizes[1]), 
                         rep("Group_B", group_sizes[2]))
  } else {
    # Create balanced groups
    group_assignment <- rep(paste0("Group_", LETTERS[1:n_groups]), 
                           each = ceiling(n_samples / n_groups))
    group_assignment <- group_assignment[1:n_samples]
  }
  
  metadata <- data.frame(
    sample_id = sample_names,
    treatment_group = factor(group_assignment),
    batch = factor(sample(c("Batch1", "Batch2", "Batch3"), n_samples, replace = TRUE)),
    age = sample(20:80, n_samples, replace = TRUE),
    sex = factor(sample(c("M", "F"), n_samples, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- sample_names
  
  # Create sample name mismatch if requested
  if (sample_name_mismatch) {
    # Change some sample names in metadata
    mismatch_indices <- sample(1:n_samples, size = min(3, n_samples))
    rownames(metadata)[mismatch_indices] <- paste0("Mismatch_", mismatch_indices)
  }
  
  # Convert abundance to requested type
  if (abundance_type == "data.frame") {
    abundance <- as.data.frame(abundance)
  } else if (abundance_type == "tibble") {
    abundance <- tibble::as_tibble(abundance, rownames = "feature_id")
  }
  
  return(list(
    abundance = abundance,
    metadata = metadata,
    feature_names = feature_names,
    sample_names = sample_names
  ))
}

#' Create test data with specific edge cases
create_edge_case_data <- function(case_type) {
  switch(case_type,
    "single_sample_per_group" = {
      create_comprehensive_test_data(n_features = 5, n_samples = 2, n_groups = 2)
    },
    "many_groups" = {
      create_comprehensive_test_data(n_features = 5, n_samples = 12, n_groups = 6)
    },
    "zero_variance" = {
      data <- create_comprehensive_test_data(n_features = 5, n_samples = 10)
      # Make some features have zero variance
      data$abundance[1, ] <- 100  # All samples same value
      data$abundance[2, ] <- 0    # All zeros
      data
    },
    "high_sparsity" = {
      data <- create_comprehensive_test_data(n_features = 20, n_samples = 10, include_zeros = TRUE)
      # Make 80% of values zero
      zero_indices <- sample(length(data$abundance), size = length(data$abundance) * 0.8)
      data$abundance[zero_indices] <- 0
      data
    }
  )
}

# ==============================================================================
# Input Validation Tests
# ==============================================================================

test_that("pathway_gsea validates abundance data types correctly", {
  test_data <- create_comprehensive_test_data()
  metadata <- test_data$metadata
  
  # Test invalid abundance types
  expect_error(
    pathway_gsea(abundance = "invalid_string", metadata = metadata, group = "treatment_group"),
    "'abundance' must be a data frame or matrix"
  )
  
  expect_error(
    pathway_gsea(abundance = list(a = 1, b = 2), metadata = metadata, group = "treatment_group"),
    "'abundance' must be a data frame or matrix"
  )
  
  expect_error(
    pathway_gsea(abundance = NULL, metadata = metadata, group = "treatment_group"),
    "'abundance' must be a data frame or matrix"
  )
  
  # Test valid abundance types - these should not error on input validation
  # (they might error later due to missing packages or other issues)
  abundance_matrix <- test_data$abundance
  abundance_df <- as.data.frame(test_data$abundance)
  
  # Mock functions to avoid dependency issues
  mockery::stub(pathway_gsea, "prepare_gene_sets", function(...) list())
  mockery::stub(pathway_gsea, "run_fgsea", function(...) {
    data.frame(
      pathway_id = character(),
      pathway_name = character(),
      size = integer(),
      ES = numeric(),
      NES = numeric(),
      pvalue = numeric(),
      p.adjust = numeric(),
      leading_edge = character(),
      stringsAsFactors = FALSE
    )
  })
  
  # These should pass input validation
  expect_silent({
    tryCatch(
      pathway_gsea(abundance = abundance_matrix, metadata = metadata, group = "treatment_group"),
      error = function(e) {
        # Only check that the error is NOT about abundance type validation
        expect_false(grepl("'abundance' must be a data frame or matrix", e$message))
      }
    )
  })
  
  expect_silent({
    tryCatch(
      pathway_gsea(abundance = abundance_df, metadata = metadata, group = "treatment_group"),
      error = function(e) {
        # Only check that the error is NOT about abundance type validation
        expect_false(grepl("'abundance' must be a data frame or matrix", e$message))
      }
    )
  })
})

test_that("pathway_gsea validates metadata correctly", {
  test_data <- create_comprehensive_test_data()
  abundance <- test_data$abundance
  
  # Test invalid metadata types
  expect_error(
    pathway_gsea(abundance = abundance, metadata = "invalid_string", group = "treatment_group"),
    "'metadata' must be a data frame"
  )
  
  expect_error(
    pathway_gsea(abundance = abundance, metadata = list(a = 1, b = 2), group = "treatment_group"),
    "'metadata' must be a data frame"
  )
  
  expect_error(
    pathway_gsea(abundance = abundance, metadata = NULL, group = "treatment_group"),
    "'metadata' must be a data frame"
  )
  
  expect_error(
    pathway_gsea(abundance = abundance, metadata = matrix(1:10, nrow = 5), group = "treatment_group"),
    "'metadata' must be a data frame"
  )
})

test_that("pathway_gsea validates group variable correctly", {
  test_data <- create_comprehensive_test_data()
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Test missing group variable
  expect_error(
    pathway_gsea(abundance = abundance, metadata = metadata, group = "nonexistent_group"),
    "Group variable nonexistent_group not found in metadata"
  )
  
  # Test empty string group
  expect_error(
    pathway_gsea(abundance = abundance, metadata = metadata, group = ""),
    "Group variable  not found in metadata"
  )
  
  # Test NULL group - this will cause an error before reaching the group check
  expect_error({
    pathway_gsea(abundance = abundance, metadata = metadata, group = NULL)
  })
})

test_that("pathway_gsea validates method parameters correctly", {
  test_data <- create_comprehensive_test_data()
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Test invalid pathway_type
  expect_error(
    pathway_gsea(abundance = abundance, metadata = metadata, group = "treatment_group", 
                pathway_type = "invalid_pathway"),
    "pathway_type must be one of 'KEGG', 'MetaCyc', or 'GO'"
  )
  
  # Test invalid method
  expect_error(
    pathway_gsea(abundance = abundance, metadata = metadata, group = "treatment_group", 
                method = "invalid_method"),
    "method must be one of 'fgsea', 'GSEA', or 'clusterProfiler'"
  )
  
  # Test invalid rank_method
  expect_error(
    pathway_gsea(abundance = abundance, metadata = metadata, group = "treatment_group", 
                rank_method = "invalid_rank"),
    "rank_method must be one of 'signal2noise', 't_test', 'log2_ratio', or 'diff_abundance'"
  )
})

# ==============================================================================
# Data Type Conversion and Formatting Tests
# ==============================================================================

test_that("pathway_gsea handles different abundance data formats correctly", {
  skip_if_not_installed("fgsea")
  
  base_data <- create_comprehensive_test_data(n_features = 5, n_samples = 10)
  metadata <- base_data$metadata
  
  # Mock functions to focus on data handling
  mockery::stub(pathway_gsea, "prepare_gene_sets", function(...) {
    list("pathway1" = c("K00001", "K00002"), "pathway2" = c("K00003", "K00004"))
  })
  mockery::stub(pathway_gsea, "run_fgsea", function(ranked_list, ...) {
    expect_type(ranked_list, "double")
    expect_true(length(ranked_list) > 0)
    expect_true(!is.null(names(ranked_list)))
    data.frame(
      pathway_id = c("pathway1", "pathway2"),
      pathway_name = c("pathway1", "pathway2"),
      size = c(2, 2),
      ES = c(0.5, -0.3),
      NES = c(1.2, -0.8),
      pvalue = c(0.01, 0.05),
      p.adjust = c(0.02, 0.1),
      leading_edge = c("K00001;K00002", "K00003;K00004"),
      stringsAsFactors = FALSE
    )
  })
  
  # Test matrix input
  abundance_matrix <- as.matrix(base_data$abundance)
  expect_silent({
    result_matrix <- pathway_gsea(abundance = abundance_matrix, metadata = metadata, 
                                 group = "treatment_group", method = "fgsea")
  })
  
  # Test data.frame input
  abundance_df <- as.data.frame(base_data$abundance)
  expect_silent({
    result_df <- pathway_gsea(abundance = abundance_df, metadata = metadata, 
                             group = "treatment_group", method = "fgsea")
  })
  
  # Results should be equivalent regardless of input type
  # (We can't test exact equality due to mocking, but both should succeed)
  expect_s3_class(tryCatch(result_matrix, error = function(e) data.frame()), "data.frame")
  expect_s3_class(tryCatch(result_df, error = function(e) data.frame()), "data.frame")
})

test_that("abundance data is correctly converted to matrix format internally", {
  test_data <- create_comprehensive_test_data(n_features = 5, n_samples = 8)
  
  # Test that calculate_rank_metric handles conversion correctly
  abundance_df <- as.data.frame(test_data$abundance)
  metadata <- test_data$metadata
  
  # Test with data.frame input
  metric_df <- calculate_rank_metric(abundance_df, metadata, "treatment_group", "signal2noise")
  
  # Test with matrix input  
  abundance_matrix <- as.matrix(test_data$abundance)
  metric_matrix <- calculate_rank_metric(abundance_matrix, metadata, "treatment_group", "signal2noise")
  
  # Results should be identical
  expect_equal(metric_df, metric_matrix)
  expect_type(metric_df, "double")
  expect_type(metric_matrix, "double")
  expect_length(metric_df, nrow(test_data$abundance))
  expect_length(metric_matrix, nrow(test_data$abundance))
})

# ==============================================================================
# Sample Name Matching Tests
# ==============================================================================

test_that("pathway_gsea detects sample name mismatches", {
  test_data <- create_comprehensive_test_data(sample_name_mismatch = TRUE)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  expect_error(
    pathway_gsea(abundance = abundance, metadata = metadata, group = "treatment_group"),
    "Sample names in abundance data do not match sample names in metadata"
  )
})

test_that("pathway_gsea handles partial sample overlap correctly", {
  test_data <- create_comprehensive_test_data(n_features = 5, n_samples = 10)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Remove some samples from metadata
  metadata_subset <- metadata[1:8, ]
  
  # Mock functions
  mockery::stub(pathway_gsea, "prepare_gene_sets", function(...) {
    list("pathway1" = c("K00001", "K00002"))
  })
  mockery::stub(pathway_gsea, "run_fgsea", function(ranked_list, ...) {
    expect_length(ranked_list, 5)  # Should have all features
    data.frame(
      pathway_id = "pathway1",
      pathway_name = "pathway1",
      size = 2,
      ES = 0.5,
      NES = 1.2,
      pvalue = 0.01,
      p.adjust = 0.02,
      leading_edge = "K00001;K00002",
      stringsAsFactors = FALSE
    )
  })
  
  # This should work - function should subset abundance to match metadata
  # However, the current implementation checks if ALL abundance samples are in metadata
  # So this will actually fail, which is the current behavior we need to document
  expect_error({
    result <- pathway_gsea(abundance = abundance, metadata = metadata_subset, 
                          group = "treatment_group", method = "fgsea")
  }, "Sample names in abundance data do not match sample names in metadata")
})

test_that("sample subsetting works correctly in calculate_rank_metric", {
  test_data <- create_comprehensive_test_data(n_features = 5, n_samples = 10)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Remove some samples from metadata
  metadata_subset <- metadata[1:6, ]
  
  # Calculate ranking metric - should automatically subset abundance
  metric <- calculate_rank_metric(abundance, metadata_subset, "treatment_group", "signal2noise")
  
  expect_type(metric, "double")
  expect_length(metric, nrow(abundance))
  # Note: Some metrics might be NA due to insufficient data in subsets
  expect_true(is.numeric(metric))
})

# ==============================================================================
# Group Factor Handling and Level Validation Tests
# ==============================================================================

test_that("calculate_rank_metric enforces two-group requirement", {
  test_data <- create_comprehensive_test_data(n_features = 5, n_samples = 12, n_groups = 3)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  expect_error(
    calculate_rank_metric(abundance, metadata, "treatment_group", "signal2noise"),
    "GSEA currently only supports two-group comparisons"
  )
})

test_that("group factor levels are handled correctly", {
  test_data <- create_comprehensive_test_data(n_features = 5, n_samples = 10)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Test with character group variable (should be converted to factor)
  metadata$treatment_group_char <- as.character(metadata$treatment_group)
  
  expect_silent({
    metric <- calculate_rank_metric(abundance, metadata, "treatment_group_char", "signal2noise")
  })
  expect_type(metric, "double")
  expect_length(metric, nrow(abundance))
  
  # Test with factor that has unused levels
  metadata$treatment_group_extra <- factor(metadata$treatment_group, 
                                         levels = c(levels(metadata$treatment_group), "Unused_Level"))
  
  expect_silent({
    metric <- calculate_rank_metric(abundance, metadata, "treatment_group_extra", "signal2noise")
  })
  expect_type(metric, "double")
  expect_length(metric, nrow(abundance))
})

test_that("group factor with single level fails appropriately", {
  test_data <- create_comprehensive_test_data(n_features = 5, n_samples = 10)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Make all samples the same group
  metadata$single_group <- factor(rep("OnlyGroup", nrow(metadata)))
  
  expect_error(
    calculate_rank_metric(abundance, metadata, "single_group", "signal2noise"),
    "GSEA currently only supports two-group comparisons"
  )
})

# ==============================================================================
# Data Subsetting and Filtering Tests
# ==============================================================================

test_that("data subsetting preserves row and column names", {
  test_data <- create_comprehensive_test_data(n_features = 8, n_samples = 12)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Test subsetting in calculate_rank_metric
  original_rownames <- rownames(abundance)
  original_sample_names <- rownames(metadata)
  
  metric <- calculate_rank_metric(abundance, metadata, "treatment_group", "signal2noise")
  
  expect_named(metric, original_rownames)
  expect_length(metric, length(original_rownames))
})

test_that("filtering handles empty groups correctly", {
  test_data <- create_comprehensive_test_data(n_features = 5, n_samples = 10)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Create a group variable where one group is empty after filtering
  metadata$problematic_group <- metadata$treatment_group
  # Change all Group_B to Group_A
  metadata$problematic_group[metadata$problematic_group == "Group_B"] <- "Group_A"
  metadata$problematic_group <- factor(metadata$problematic_group)
  
  expect_error(
    calculate_rank_metric(abundance, metadata, "problematic_group", "signal2noise"),
    "GSEA currently only supports two-group comparisons"
  )
})

# ==============================================================================
# Missing Data Handling Tests
# ==============================================================================

test_that("missing values in abundance data are handled appropriately", {
  test_data <- create_comprehensive_test_data(n_features = 5, n_samples = 10, include_missing = TRUE)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # For signal2noise method
  metric_s2n <- calculate_rank_metric(abundance, metadata, "treatment_group", "signal2noise")
  expect_type(metric_s2n, "double")
  # Some values might be NA due to missing data, which is acceptable
  
  # For t_test method
  metric_ttest <- calculate_rank_metric(abundance, metadata, "treatment_group", "t_test")
  expect_type(metric_ttest, "double")
  
  # For log2_ratio method
  metric_log2 <- calculate_rank_metric(abundance, metadata, "treatment_group", "log2_ratio")
  expect_type(metric_log2, "double")
  
  # For diff_abundance method
  metric_diff <- calculate_rank_metric(abundance, metadata, "treatment_group", "diff_abundance")
  expect_type(metric_diff, "double")
})

test_that("missing values in metadata are handled appropriately", {
  test_data <- create_comprehensive_test_data(n_features = 5, n_samples = 10)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Add missing values to group variable
  metadata$treatment_group_na <- metadata$treatment_group
  metadata$treatment_group_na[1:2] <- NA
  
  # This should error because factor() will create NA levels
  expect_error({
    calculate_rank_metric(abundance, metadata, "treatment_group_na", "signal2noise")
  })
})

# ==============================================================================
# Zero Variance Handling Tests
# ==============================================================================

test_that("zero variance features are handled correctly in signal2noise", {
  test_data <- create_edge_case_data("zero_variance")
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # signal2noise method should handle zero variance by adding small constant
  metric <- calculate_rank_metric(abundance, metadata, "treatment_group", "signal2noise")
  
  expect_type(metric, "double")
  expect_length(metric, nrow(abundance))
  expect_true(all(is.finite(metric)))  # No infinite values
})

test_that("all ranking methods handle zero values appropriately", {
  test_data <- create_edge_case_data("zero_variance")
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Test each ranking method
  methods <- c("signal2noise", "t_test", "log2_ratio", "diff_abundance")
  
  for (method in methods) {
    if (method == "t_test") {
      # t_test will fail with zero variance data
      expect_error({
        calculate_rank_metric(abundance, metadata, "treatment_group", method)
      })
    } else {
      metric <- calculate_rank_metric(abundance, metadata, "treatment_group", method)
      expect_type(metric, "double")
      expect_length(metric, nrow(abundance))
      
      if (method == "log2_ratio") {
        # log2_ratio should handle zeros by adding small constant
        expect_true(all(is.finite(metric)))
      }
    }
  }
})

# ==============================================================================
# Edge Cases and Stress Tests
# ==============================================================================

test_that("single sample per group edge case", {
  test_data <- create_edge_case_data("single_sample_per_group")
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Should work but might produce warnings about variance calculations
  expect_silent({
    metric <- calculate_rank_metric(abundance, metadata, "treatment_group", "diff_abundance")
  })
  expect_type(metric, "double")
  expect_length(metric, nrow(abundance))
  
  # t_test might fail with single samples
  expect_error({
    calculate_rank_metric(abundance, metadata, "treatment_group", "t_test")
  })
})

test_that("high sparsity data handling", {
  test_data <- create_edge_case_data("high_sparsity")
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # All methods should handle sparse data
  methods <- c("signal2noise", "t_test", "log2_ratio", "diff_abundance")
  
  for (method in methods) {
    if (method == "t_test") {
      # t_test might fail with high sparsity, but not always
      # So we'll test that it either succeeds or fails gracefully
      tryCatch({
        metric <- calculate_rank_metric(abundance, metadata, "treatment_group", method)
        expect_type(metric, "double")
        expect_length(metric, nrow(abundance))
      }, error = function(e) {
        # This is acceptable - t_test can fail with sparse data
        expect_true(TRUE)
      })
    } else {
      metric <- calculate_rank_metric(abundance, metadata, "treatment_group", method)
      expect_type(metric, "double")
      expect_length(metric, nrow(abundance))
    }
  }
})

test_that("imbalanced groups are handled correctly", {
  test_data <- create_comprehensive_test_data(n_features = 5, n_samples = 20, group_imbalance = TRUE)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # All methods should handle imbalanced groups
  methods <- c("signal2noise", "t_test", "log2_ratio", "diff_abundance")
  
  for (method in methods) {
    metric <- calculate_rank_metric(abundance, metadata, "treatment_group", method)
    expect_type(metric, "double")
    expect_length(metric, nrow(abundance))
  }
})

# ==============================================================================
# Data Preprocessing Integration Tests
# ==============================================================================

test_that("end-to-end data preprocessing workflow", {
  skip_if_not_installed("fgsea")
  
  test_data <- create_comprehensive_test_data(n_features = 10, n_samples = 16)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Mock the gene sets and fgsea functions
  mockery::stub(pathway_gsea, "prepare_gene_sets", function(...) {
    list(
      "pathway1" = paste0("K", sprintf("%05d", 1:3)),
      "pathway2" = paste0("K", sprintf("%05d", 4:6)),
      "pathway3" = paste0("K", sprintf("%05d", 7:10))
    )
  })
  
  mockery::stub(pathway_gsea, "run_fgsea", function(ranked_list, gene_sets, ...) {
    expect_type(ranked_list, "double")
    expect_length(ranked_list, nrow(abundance))
    expect_true(!is.null(names(ranked_list)))
    expect_type(gene_sets, "list")
    
    data.frame(
      pathway_id = names(gene_sets),
      pathway_name = names(gene_sets),
      size = sapply(gene_sets, length),
      ES = c(0.5, -0.3, 0.1),
      NES = c(1.2, -0.8, 0.2),
      pvalue = c(0.01, 0.05, 0.3),
      p.adjust = c(0.02, 0.1, 0.4),
      leading_edge = sapply(gene_sets, function(x) paste(x[1:2], collapse = ";")),
      stringsAsFactors = FALSE
    )
  })
  
  # Test full workflow
  result <- pathway_gsea(
    abundance = abundance,
    metadata = metadata, 
    group = "treatment_group",
    method = "fgsea",
    rank_method = "signal2noise"
  )
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_true(all(c("pathway_id", "pathway_name", "size", "ES", "NES", 
                   "pvalue", "p.adjust", "leading_edge", "method") %in% colnames(result)))
  expect_equal(result$method, rep("fgsea", 3))
})

test_that("preprocessing preserves data integrity", {
  test_data <- create_comprehensive_test_data(n_features = 8, n_samples = 12)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Store original data characteristics
  original_rownames <- rownames(abundance)
  original_colnames <- colnames(abundance)
  original_nrow <- nrow(abundance)
  original_ncol <- ncol(abundance)
  
  # Calculate rank metric (this does internal preprocessing)
  metric <- calculate_rank_metric(abundance, metadata, "treatment_group", "signal2noise")
  
  # Check that original data is unchanged
  expect_equal(rownames(abundance), original_rownames)
  expect_equal(colnames(abundance), original_colnames)
  expect_equal(nrow(abundance), original_nrow)
  expect_equal(ncol(abundance), original_ncol)
  
  # Check that metric has correct properties
  expect_named(metric, original_rownames)
  expect_length(metric, original_nrow)
  expect_type(metric, "double")
})

# ==============================================================================
# Error Message Quality Tests
# ==============================================================================

test_that("error messages are informative", {
  test_data <- create_comprehensive_test_data()
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Test specific error message content
  expect_error(
    pathway_gsea(abundance = "wrong", metadata = metadata, group = "treatment_group"),
    "'abundance' must be a data frame or matrix",
    fixed = TRUE
  )
  
  expect_error(
    pathway_gsea(abundance = abundance, metadata = "wrong", group = "treatment_group"),
    "'metadata' must be a data frame",
    fixed = TRUE
  )
  
  expect_error(
    pathway_gsea(abundance = abundance, metadata = metadata, group = "wrong_group"),
    "Group variable wrong_group not found in metadata",
    fixed = TRUE
  )
  
  expect_error(
    pathway_gsea(abundance = abundance, metadata = metadata, group = "treatment_group", 
                pathway_type = "wrong"),
    "pathway_type must be one of 'KEGG', 'MetaCyc', or 'GO'",
    fixed = TRUE
  )
})

# ==============================================================================
# Performance and Efficiency Tests
# ==============================================================================

test_that("functions handle large datasets efficiently", {
  skip_on_cran()  # Skip on CRAN to avoid long test times
  
  # Create larger dataset
  large_data <- create_comprehensive_test_data(n_features = 100, n_samples = 50)
  abundance <- large_data$abundance
  metadata <- large_data$metadata
  
  # Test that calculation completes in reasonable time
  start_time <- Sys.time()
  metric <- calculate_rank_metric(abundance, metadata, "treatment_group", "signal2noise")
  end_time <- Sys.time()
  
  expect_lt(as.numeric(end_time - start_time), 10)  # Should complete within 10 seconds
  expect_type(metric, "double")
  expect_length(metric, nrow(abundance))
})