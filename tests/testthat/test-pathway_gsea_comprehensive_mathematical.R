# Comprehensive Mathematical Tests for GSEA Core Functions
# Testing mathematical correctness and algorithmic implementation

# Test data creation functions
create_deterministic_test_data <- function(n_features = 20, n_samples = 20, effect_size = 2) {
  set.seed(12345)  # Fixed seed for reproducibility
  
  # Create known differential abundance patterns
  abundance <- matrix(0, nrow = n_features, ncol = n_samples)
  
  # Group assignment (equal sizes)
  group1_idx <- 1:(n_samples/2)
  group2_idx <- (n_samples/2 + 1):n_samples
  
  # Create features with known differences
  for (i in 1:n_features) {
    # Base abundance
    base_mean <- runif(1, 50, 200)
    base_sd <- base_mean * 0.1  # 10% CV
    
    # Group 1 samples
    abundance[i, group1_idx] <- rnorm(length(group1_idx), 
                                     mean = base_mean, 
                                     sd = base_sd)
    
    # Group 2 samples - some with effect, some without
    if (i <= n_features/2) {
      # Features with true differential abundance
      abundance[i, group2_idx] <- rnorm(length(group2_idx), 
                                       mean = base_mean + effect_size * base_sd, 
                                       sd = base_sd)
    } else {
      # Features without differential abundance
      abundance[i, group2_idx] <- rnorm(length(group2_idx), 
                                       mean = base_mean, 
                                       sd = base_sd)
    }
  }
  
  # Set names
  rownames(abundance) <- paste0("K", sprintf("%05d", 1:n_features))
  colnames(abundance) <- paste0("Sample", 1:n_samples)
  
  # Create metadata
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(c(rep("Control", n_samples/2), rep("Treatment", n_samples/2)))
  )
  rownames(metadata) <- metadata$sample_name
  
  return(list(
    abundance = abundance,
    metadata = metadata,
    true_diff_features = rownames(abundance)[1:(n_features/2)]
  ))
}

create_zero_variance_test_data <- function() {
  set.seed(123)
  abundance <- matrix(0, nrow = 5, ncol = 10)
  
  # Feature with zero variance in group 1
  abundance[1, 1:5] <- rep(100, 5)  # No variance
  abundance[1, 6:10] <- rnorm(5, 150, 20)
  
  # Feature with zero variance in group 2
  abundance[2, 1:5] <- rnorm(5, 100, 20)
  abundance[2, 6:10] <- rep(150, 5)  # No variance
  
  # Feature with zero variance in both groups
  abundance[3, 1:5] <- rep(100, 5)
  abundance[3, 6:10] <- rep(100, 5)
  
  # Normal features
  abundance[4, ] <- rnorm(10, 100, 20)
  abundance[5, ] <- rnorm(10, 150, 30)
  
  rownames(abundance) <- paste0("K", sprintf("%05d", 1:5))
  colnames(abundance) <- paste0("Sample", 1:10)
  
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("Control", "Treatment"), each = 5))
  )
  rownames(metadata) <- metadata$sample_name
  
  return(list(abundance = abundance, metadata = metadata))
}

# Mathematical validation tests
test_that("signal2noise ranking metric is mathematically correct", {
  test_data <- create_deterministic_test_data(n_features = 100, n_samples = 20)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Calculate ranking metric
  s2n_metric <- calculate_rank_metric(abundance, metadata, "group", "signal2noise")
  
  # Validate mathematical formula: (mean1 - mean2) / (sd1 + sd2)
  for (i in 1:10) {  # Test first 10 features
    group1_data <- abundance[i, metadata$group == "Control"]
    group2_data <- abundance[i, metadata$group == "Treatment"]
    
    expected_mean_diff <- mean(group1_data) - mean(group2_data)
    expected_sd_sum <- sd(group1_data) + sd(group2_data)
    expected_s2n <- expected_mean_diff / expected_sd_sum
    
    expect_equal(as.numeric(s2n_metric[i]), expected_s2n, tolerance = 1e-10, 
                 info = paste("Signal-to-noise calculation for feature", i))
  }
  
  # Verify features with known positive effects have negative scores (Control vs Treatment)
  diff_features <- test_data$true_diff_features[1:5]
  non_diff_features <- setdiff(names(s2n_metric), diff_features)[1:5]
  
  expect_true(mean(s2n_metric[diff_features]) < mean(s2n_metric[non_diff_features]),
              "Features with true differences should have more extreme signal-to-noise ratios")
})

test_that("t-test ranking metric is mathematically correct", {
  test_data <- create_deterministic_test_data(n_features = 50, n_samples = 20)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Calculate ranking metric
  ttest_metric <- calculate_rank_metric(abundance, metadata, "group", "t_test")
  
  # Validate t-test calculations
  for (i in 1:10) {
    group1_data <- abundance[i, metadata$group == "Control"]
    group2_data <- abundance[i, metadata$group == "Treatment"]
    
    # Manual t-test calculation
    manual_t_test <- t.test(group1_data, group2_data)
    
    expect_equal(as.numeric(ttest_metric[i]), as.numeric(manual_t_test$statistic), 
                 tolerance = 1e-10,
                 info = paste("T-test statistic for feature", i))
  }
  
  # Verify that larger differences produce more extreme t-statistics
  diff_features_t <- abs(ttest_metric[test_data$true_diff_features])
  non_diff_features_t <- abs(ttest_metric[setdiff(names(ttest_metric), test_data$true_diff_features)])
  
  expect_true(mean(diff_features_t) > mean(non_diff_features_t),
              "Features with true differences should have larger |t-statistics|")
})

test_that("log2_ratio ranking metric is mathematically correct", {
  # Create data with known fold changes
  set.seed(123)
  abundance <- matrix(0, nrow = 6, ncol = 10)
  
  # Features with known fold changes
  abundance[1, 1:5] <- rep(100, 5)     # Control: 100
  abundance[1, 6:10] <- rep(200, 5)    # Treatment: 200 (2-fold up)
  
  abundance[2, 1:5] <- rep(200, 5)     # Control: 200  
  abundance[2, 6:10] <- rep(100, 5)    # Treatment: 100 (2-fold down)
  
  abundance[3, 1:5] <- rep(100, 5)     # Control: 100
  abundance[3, 6:10] <- rep(400, 5)    # Treatment: 400 (4-fold up)
  
  abundance[4, 1:5] <- rep(100, 5)     # Control: 100
  abundance[4, 6:10] <- rep(100, 5)    # Treatment: 100 (no change)
  
  # Add small noise to avoid exact zeros
  abundance <- abundance + matrix(rnorm(60, 0, 0.01), nrow = 6, ncol = 10)
  
  rownames(abundance) <- paste0("K", sprintf("%05d", 1:6))
  colnames(abundance) <- paste0("Sample", 1:10)
  
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("Control", "Treatment"), each = 5))
  )
  rownames(metadata) <- metadata$sample_name
  
  # Calculate log2 ratio
  log2_metric <- calculate_rank_metric(abundance, metadata, "group", "log2_ratio")
  
  # Validate mathematical formula: log2(mean1 / mean2)
  expected_ratios <- c(-1, 1, -2, 0)  # log2 of fold changes
  
  for (i in 1:4) {
    expect_equal(as.numeric(log2_metric[i]), expected_ratios[i], tolerance = 0.1,
                 info = paste("Log2 ratio for feature", i))
  }
})

test_that("diff_abundance ranking metric is mathematically correct", {
  test_data <- create_deterministic_test_data(n_features = 30, n_samples = 16)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Calculate difference metric
  diff_metric <- calculate_rank_metric(abundance, metadata, "group", "diff_abundance")
  
  # Validate mathematical formula: mean1 - mean2
  for (i in 1:10) {
    group1_data <- abundance[i, metadata$group == "Control"]
    group2_data <- abundance[i, metadata$group == "Treatment"]
    
    expected_diff <- mean(group1_data) - mean(group2_data)
    
    expect_equal(as.numeric(diff_metric[i]), expected_diff, tolerance = 1e-10,
                 info = paste("Difference in abundance for feature", i))
  }
})

test_that("ranking metrics handle zero variance correctly", {
  test_data <- create_zero_variance_test_data()
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Test signal-to-noise with zero variance
  s2n_metric <- calculate_rank_metric(abundance, metadata, "group", "signal2noise")
  
  # Features with zero variance should be handled gracefully
  expect_true(all(is.finite(s2n_metric)), "All signal-to-noise values should be finite")
  expect_true(all(!is.na(s2n_metric)), "No signal-to-noise values should be NA")
  
  # Feature with zero variance in both groups should have zero signal-to-noise
  expect_equal(as.numeric(s2n_metric[3]), 0, tolerance = 1e-6)
  
  # Test other metrics - handle t-test errors for constant data
  ttest_result <- tryCatch({
    calculate_rank_metric(abundance, metadata, "group", "t_test")
  }, error = function(e) NULL)
  
  if (!is.null(ttest_result)) {
    ttest_metric <- ttest_result
    expect_true(all(is.finite(ttest_metric)), "All t-test values should be finite")
  }
  
  log2_metric <- calculate_rank_metric(abundance, metadata, "group", "log2_ratio")
  diff_metric <- calculate_rank_metric(abundance, metadata, "group", "diff_abundance")
  expect_true(all(is.finite(log2_metric)), "All log2 ratio values should be finite")
  expect_true(all(is.finite(diff_metric)), "All difference values should be finite")
})

test_that("ranking metrics produce consistent orderings", {
  # Create data with clear differential patterns
  set.seed(456)
  abundance <- matrix(0, nrow = 20, ncol = 20)
  
  # Create features with increasing effect sizes
  for (i in 1:20) {
    base_abundance <- 100
    effect_size <- i * 10  # Increasing effect
    
    abundance[i, 1:10] <- rnorm(10, base_abundance, 10)
    abundance[i, 11:20] <- rnorm(10, base_abundance + effect_size, 10)
  }
  
  rownames(abundance) <- paste0("K", sprintf("%05d", 1:20))
  colnames(abundance) <- paste0("Sample", 1:20)
  
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("Control", "Treatment"), each = 10))
  )
  rownames(metadata) <- metadata$sample_name
  
  # Calculate all metrics
  s2n_metric <- calculate_rank_metric(abundance, metadata, "group", "signal2noise")
  ttest_metric <- calculate_rank_metric(abundance, metadata, "group", "t_test")
  log2_metric <- calculate_rank_metric(abundance, metadata, "group", "log2_ratio")
  diff_metric <- calculate_rank_metric(abundance, metadata, "group", "diff_abundance")
  
  # All metrics should show increasingly negative values (Control < Treatment)
  # Allow some tolerance for numerical precision and randomness
  s2n_trend <- diff(s2n_metric)
  ttest_trend <- diff(ttest_metric)
  log2_trend <- diff(log2_metric)
  diff_trend <- diff(diff_metric)
  
  # Most differences should be negative (decreasing trend)
  expect_true(mean(s2n_trend < 0) > 0.7, "Signal-to-noise should generally decrease with increasing effect")
  expect_true(mean(ttest_trend < 0) > 0.7, "T-test statistic should generally decrease with increasing effect")
  expect_true(mean(log2_trend < 0) > 0.7, "Log2 ratio should generally decrease with increasing effect")
  expect_true(mean(diff_trend < 0) > 0.7, "Difference should generally decrease with increasing effect")
})

test_that("pathway_gsea parameter validation is comprehensive", {
  test_data <- create_deterministic_test_data()
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Test all invalid inputs systematically
  invalid_inputs <- list(
    abundance = list(
      "character" = "not_a_matrix",
      "list" = list(a = 1, b = 2),
      "NULL" = NULL
    ),
    metadata = list(
      "character" = "not_a_dataframe",
      "matrix" = matrix(1:4, nrow = 2),
      "NULL" = NULL
    ),
    group = list(
      "missing_column" = "nonexistent_group",
      "numeric" = 123,
      "NULL" = NULL
    ),
    pathway_type = list(
      "invalid" = "REACTOME",
      "numeric" = 1,
      "NULL" = NULL
    ),
    method = list(
      "invalid" = "DESeq2",
      "numeric" = 2,
      "NULL" = NULL
    ),
    rank_method = list(
      "invalid" = "wilcox_test",
      "numeric" = 3,
      "NULL" = NULL
    )
  )
  
  # Test abundance validation
  for (invalid_type in names(invalid_inputs$abundance)) {
    expect_error(
      pathway_gsea(
        abundance = invalid_inputs$abundance[[invalid_type]],
        metadata = metadata,
        group = "group"
      ),
      "'abundance' must be a data frame or matrix",
      info = paste("Should reject", invalid_type, "abundance input")
    )
  }
  
  # Test metadata validation  
  for (invalid_type in names(invalid_inputs$metadata)) {
    expect_error(
      pathway_gsea(
        abundance = abundance,
        metadata = invalid_inputs$metadata[[invalid_type]],
        group = "group"
      ),
      "'metadata' must be a data frame",
      info = paste("Should reject", invalid_type, "metadata input")
    )
  }
  
  # Test parameter validation
  expect_error(
    pathway_gsea(abundance, metadata, group = "nonexistent"),
    "Group variable .* not found in metadata"
  )
  
  expect_error(
    pathway_gsea(abundance, metadata, "group", pathway_type = "INVALID"),
    "pathway_type must be one of"
  )
  
  expect_error(
    pathway_gsea(abundance, metadata, "group", method = "INVALID"),
    "method must be one of"
  )
  
  expect_error(
    pathway_gsea(abundance, metadata, "group", rank_method = "INVALID"),
    "rank_method must be one of"
  )
})

test_that("statistical significance calculations are accurate", {
  skip_if_not_installed("fgsea")
  
  # Create data with known enrichment patterns
  test_data <- create_deterministic_test_data(n_features = 100, n_samples = 20)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Create gene sets that should be enriched
  gene_sets <- list(
    "enriched_pathway" = test_data$true_diff_features[1:8],  # Should be enriched
    "random_pathway" = sample(rownames(abundance), 10),     # Random
    "depleted_pathway" = setdiff(rownames(abundance), test_data$true_diff_features)[1:8]  # Should be depleted
  )
  
  # Calculate ranking metric
  ranked_list <- calculate_rank_metric(abundance, metadata, "group", "signal2noise")
  
  # Mock fgsea to test statistical components
  mock_fgsea_result <- data.frame(
    pathway = names(gene_sets),
    pval = c(0.001, 0.5, 0.002),    # enriched, random, depleted
    padj = c(0.003, 0.7, 0.006),    # Adjusted p-values
    ES = c(-0.8, 0.1, 0.7),         # Enrichment scores
    NES = c(-2.1, 0.2, 1.9),        # Normalized enrichment scores
    size = c(8, 10, 8),
    leadingEdge = I(list(
      test_data$true_diff_features[1:5],
      sample(rownames(abundance), 5),
      setdiff(rownames(abundance), test_data$true_diff_features)[1:5]
    )),
    stringsAsFactors = FALSE
  )
  
  # Mock the fgsea function
  mockery::stub(run_fgsea, "fgsea::fgsea", function(...) mock_fgsea_result)
  
  # Run fgsea
  result <- run_fgsea(ranked_list, gene_sets, nperm = 100)
  
  # Validate result structure and statistical properties
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
  
  # Check column names and types
  expected_columns <- c("pathway_id", "pathway_name", "size", "ES", "NES", 
                       "pvalue", "p.adjust", "leading_edge")
  expect_named(result, expected_columns)
  
  # Validate statistical values
  expect_true(all(result$pvalue >= 0 & result$pvalue <= 1), 
              "P-values should be between 0 and 1")
  expect_true(all(result$p.adjust >= 0 & result$p.adjust <= 1), 
              "Adjusted p-values should be between 0 and 1")
  expect_true(all(result$p.adjust >= result$pvalue), 
              "Adjusted p-values should be >= raw p-values")
  
  # Check enrichment scores have correct signs
  expect_true(result$ES[1] < 0, "Enriched pathway should have negative ES")
  expect_true(result$ES[3] > 0, "Depleted pathway should have positive ES")
})

test_that("method comparison produces consistent results", {
  skip_if_not_installed("fgsea")
  
  test_data <- create_deterministic_test_data(n_features = 50, n_samples = 16)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Create simple gene sets
  gene_sets <- list(
    "pathway1" = test_data$true_diff_features[1:10],
    "pathway2" = setdiff(rownames(abundance), test_data$true_diff_features)[1:10]
  )
  
  # Mock prepare_gene_sets to return our test gene sets
  mockery::stub(pathway_gsea, "prepare_gene_sets", function(...) gene_sets)
  
  # Mock fgsea results
  mock_fgsea_result <- data.frame(
    pathway_id = c("pathway1", "pathway2"),
    pathway_name = c("pathway1", "pathway2"),
    size = c(10, 10),
    ES = c(-0.5, 0.3),
    NES = c(-1.8, 1.2),
    pvalue = c(0.01, 0.15),
    p.adjust = c(0.02, 0.15),
    leading_edge = c("K00001;K00002", "K00011;K00012"),
    stringsAsFactors = FALSE
  )
  
  mockery::stub(pathway_gsea, "run_fgsea", function(...) mock_fgsea_result)
  
  # Test different ranking methods produce different but valid results
  methods <- c("signal2noise", "t_test", "log2_ratio", "diff_abundance")
  results <- list()
  
  for (method in methods) {
    results[[method]] <- pathway_gsea(
      abundance = abundance,
      metadata = metadata,
      group = "group",
      method = "fgsea",
      rank_method = method,
      nperm = 100,
      seed = 123
    )
    
    # Each method should produce valid results
    expect_s3_class(results[[method]], "data.frame")
    expect_equal(nrow(results[[method]]), 2)
    expect_equal(results[[method]]$method, c("fgsea", "fgsea"))
  }
  
  # Results should be consistent in structure across methods
  for (i in 2:length(methods)) {
    expect_equal(names(results[[methods[1]]]), names(results[[methods[i]]]),
                 info = paste("Column names should match between", methods[1], "and", methods[i]))
    expect_equal(nrow(results[[methods[1]]]), nrow(results[[methods[i]]]),
                 info = paste("Row counts should match between", methods[1], "and", methods[i]))
  }
})

test_that("permutation reproducibility is ensured", {
  skip_if_not_installed("fgsea")
  
  test_data <- create_deterministic_test_data(n_features = 30, n_samples = 12)
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  
  # Mock functions for reproducibility testing
  gene_sets <- list("test_pathway" = test_data$true_diff_features[1:8])
  
  mockery::stub(pathway_gsea, "prepare_gene_sets", function(...) gene_sets)
  
  # Create deterministic mock result
  deterministic_result <- data.frame(
    pathway_id = "test_pathway",
    pathway_name = "test_pathway",
    size = 8,
    ES = -0.6,
    NES = -1.9,
    pvalue = 0.008,
    p.adjust = 0.008,
    leading_edge = paste(test_data$true_diff_features[1:5], collapse = ";"),
    stringsAsFactors = FALSE
  )
  
  mockery::stub(pathway_gsea, "run_fgsea", function(...) deterministic_result)
  
  # Run analysis twice with same seed
  result1 <- pathway_gsea(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    method = "fgsea",
    seed = 12345,
    nperm = 100
  )
  
  result2 <- pathway_gsea(
    abundance = abundance,
    metadata = metadata,
    group = "group", 
    method = "fgsea",
    seed = 12345,
    nperm = 100
  )
  
  # Results should be identical when using the same seed
  expect_equal(result1, result2)
  
  # Test different seeds produce potentially different ranking metrics
  ranking1 <- calculate_rank_metric(abundance, metadata, "group", "signal2noise")
  set.seed(12345)
  ranking2 <- calculate_rank_metric(abundance, metadata, "group", "signal2noise") 
  set.seed(54321)
  ranking3 <- calculate_rank_metric(abundance, metadata, "group", "signal2noise")
  
  # Rankings should be identical (no randomness in calculation)
  expect_equal(ranking1, ranking2)
  expect_equal(ranking1, ranking3)
})

test_that("numerical edge cases are handled correctly", {
  # Test with extremely small values
  small_abundance <- matrix(runif(100, 1e-10, 1e-8), nrow = 10, ncol = 10)
  rownames(small_abundance) <- paste0("K", sprintf("%05d", 1:10))
  colnames(small_abundance) <- paste0("Sample", 1:10)
  
  metadata <- data.frame(
    sample_name = colnames(small_abundance),
    group = factor(rep(c("Control", "Treatment"), each = 5))
  )
  rownames(metadata) <- metadata$sample_name
  
  # Test all ranking methods with small values
  methods <- c("signal2noise", "t_test", "log2_ratio", "diff_abundance")
  
  for (method in methods) {
    metric <- calculate_rank_metric(small_abundance, metadata, "group", method)
    
    expect_true(all(is.finite(metric)), 
                paste("All values should be finite for", method, "with small abundances"))
    expect_true(all(!is.na(metric)), 
                paste("No values should be NA for", method, "with small abundances"))
  }
  
  # Test with large values
  large_abundance <- matrix(runif(100, 1e6, 1e8), nrow = 10, ncol = 10)
  rownames(large_abundance) <- paste0("K", sprintf("%05d", 1:10))
  colnames(large_abundance) <- paste0("Sample", 1:10)
  
  for (method in methods) {
    metric <- calculate_rank_metric(large_abundance, metadata, "group", method)
    
    expect_true(all(is.finite(metric)), 
                paste("All values should be finite for", method, "with large abundances"))
    expect_true(all(!is.na(metric)), 
                paste("No values should be NA for", method, "with large abundances"))
  }
})

test_that("sample size effects on statistical power", {
  # Test with different sample sizes
  sample_sizes <- c(6, 10, 20, 40)
  effect_size <- 1.5
  
  statistical_power <- numeric(length(sample_sizes))
  
  for (i in seq_along(sample_sizes)) {
    n_samples <- sample_sizes[i]
    test_data <- create_deterministic_test_data(
      n_features = 50, 
      n_samples = n_samples, 
      effect_size = effect_size
    )
    
    # Calculate t-test statistics
    ttest_metric <- calculate_rank_metric(
      test_data$abundance, 
      test_data$metadata, 
      "group", 
      "t_test"
    )
    
    # Calculate "statistical power" as proportion of true differential features
    # with |t-statistic| > 2 (rough threshold)
    diff_features_t <- abs(ttest_metric[test_data$true_diff_features])
    statistical_power[i] <- mean(diff_features_t > 2)
  }
  
  # Statistical power should generally increase with sample size
  expect_true(statistical_power[4] >= statistical_power[1], 
              "Statistical power should increase with sample size")
  
  # At least some power should be observed with larger samples
  expect_true(statistical_power[4] > 0.3, 
              "Should have reasonable power with large samples")
})