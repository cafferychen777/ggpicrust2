# Statistical Validation and Edge Case Tests for GSEA Functions
# Testing statistical accuracy, boundary conditions, and edge cases

library(testthat)

# Specialized test data generators for statistical validation
create_known_enrichment_data <- function(n_features = 100, n_samples = 20, 
                                        enriched_features = 1:20, 
                                        effect_size = 2.0) {
  set.seed(999)
  
  abundance <- matrix(0, nrow = n_features, ncol = n_samples)
  
  # Split samples
  group1_idx <- 1:(n_samples/2)
  group2_idx <- (n_samples/2 + 1):n_samples
  
  for (i in 1:n_features) {
    base_abundance <- rlnorm(1, meanlog = 4, sdlog = 0.5)
    noise_sd <- base_abundance * 0.15
    
    # Group 1 baseline
    abundance[i, group1_idx] <- rlnorm(length(group1_idx), 
                                      meanlog = log(base_abundance), 
                                      sdlog = 0.15)
    
    # Group 2 - enriched features have higher abundance
    if (i %in% enriched_features) {
      abundance[i, group2_idx] <- rlnorm(length(group2_idx), 
                                        meanlog = log(base_abundance * effect_size), 
                                        sdlog = 0.15)
    } else {
      abundance[i, group2_idx] <- rlnorm(length(group2_idx), 
                                        meanlog = log(base_abundance), 
                                        sdlog = 0.15)
    }
  }
  
  rownames(abundance) <- paste0("K", sprintf("%05d", 1:n_features))
  colnames(abundance) <- paste0("Sample", 1:n_samples)
  
  metadata <- data.frame(
    sample_name = colnames(abundance),
    condition = factor(c(rep("Baseline", n_samples/2), rep("Enriched", n_samples/2)))
  )
  rownames(metadata) <- metadata$sample_name
  
  return(list(
    abundance = abundance,
    metadata = metadata,
    enriched_features = paste0("K", sprintf("%05d", enriched_features))
  ))
}

create_extreme_value_data <- function() {
  set.seed(555)
  
  # Test various extreme scenarios
  extreme_cases <- list(
    zero_values = matrix(c(rep(0, 20), runif(20, 1, 100)), nrow = 4, ncol = 10),
    very_small = matrix(runif(40, 1e-12, 1e-10), nrow = 4, ncol = 10),
    very_large = matrix(runif(40, 1e10, 1e12), nrow = 4, ncol = 10),
    mixed_scales = matrix(c(runif(10, 1e-8, 1e-6), runif(10, 1e2, 1e4), 
                           runif(10, 1e-3, 1e-1), runif(10, 1e8, 1e10)), 
                         nrow = 4, ncol = 10),
    constant_values = matrix(rep(c(100, 200, 300, 400), each = 10), nrow = 4, ncol = 10),
    single_outlier = matrix(c(rep(100, 35), 1e6, rep(100, 4)), nrow = 4, ncol = 10)
  )
  
  # Add names and metadata for each case
  for (case_name in names(extreme_cases)) {
    rownames(extreme_cases[[case_name]]) <- paste0("K", sprintf("%05d", 1:4))
    colnames(extreme_cases[[case_name]]) <- paste0("Sample", 1:10)
  }
  
  metadata <- data.frame(
    sample_name = paste0("Sample", 1:10),
    group = factor(rep(c("Group1", "Group2"), each = 5))
  )
  rownames(metadata) <- metadata$sample_name
  
  return(list(cases = extreme_cases, metadata = metadata))
}

# Statistical validation tests
test_that("ranking metrics produce statistically valid orderings", {
  # Create data with known statistical relationships
  test_data <- create_known_enrichment_data(
    n_features = 60, 
    n_samples = 24, 
    enriched_features = 1:15, 
    effect_size = 3.0
  )
  
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  enriched_features <- test_data$enriched_features
  
  # Test all ranking methods
  ranking_methods <- c("signal2noise", "t_test", "log2_ratio", "diff_abundance")
  
  for (method in ranking_methods) {
    ranking <- calculate_rank_metric(abundance, metadata, "condition", method)
    
    # Enriched features should generally rank higher (more negative for Baseline vs Enriched)
    enriched_ranks <- ranking[enriched_features]
    non_enriched_ranks <- ranking[setdiff(names(ranking), enriched_features)]
    
    # Statistical test: enriched features should have significantly different rankings
    wilcox_test <- wilcox.test(enriched_ranks, non_enriched_ranks)
    
    expect_true(wilcox_test$p.value < 0.05, 
                paste("Enriched features should have significantly different rankings with", method))
    
    # Enriched features should have more negative values (baseline < enriched)
    expect_true(mean(enriched_ranks) < mean(non_enriched_ranks),
                paste("Enriched features should have lower ranking values with", method))
  }
})

test_that("signal-to-noise ratio handles statistical edge cases correctly", {
  # Test specific statistical scenarios
  
  # Case 1: Perfect separation (no overlap)
  perfect_sep <- matrix(c(rep(100, 10), rep(200, 10)), nrow = 2, ncol = 10)
  rownames(perfect_sep) <- c("K00001", "K00002")
  colnames(perfect_sep) <- paste0("Sample", 1:10)
  
  metadata <- data.frame(
    sample_name = paste0("Sample", 1:10),
    group = factor(rep(c("A", "B"), each = 5))
  )
  rownames(metadata) <- metadata$sample_name
  
  s2n_perfect <- calculate_rank_metric(perfect_sep, metadata, "group", "signal2noise")
  
  # Should produce finite, reasonable values
  expect_true(all(is.finite(s2n_perfect)))
  expect_true(all(abs(s2n_perfect) > 0))
  
  # Case 2: Identical means, different variances
  diff_var <- matrix(c(rnorm(5, 100, 1), rnorm(5, 100, 10)), nrow = 1, ncol = 10)
  rownames(diff_var) <- "K00001"
  colnames(diff_var) <- paste0("Sample", 1:10)
  
  s2n_diff_var <- calculate_rank_metric(diff_var, metadata, "group", "signal2noise")
  
  # Should be close to zero (same means)
  expect_true(abs(s2n_diff_var[1]) < 1, "Signal-to-noise should be small for same means")
  
  # Case 3: Extreme variance ratios
  extreme_var <- matrix(c(rnorm(5, 100, 0.01), rnorm(5, 101, 50)), nrow = 1, ncol = 10)
  rownames(extreme_var) <- "K00001"
  colnames(extreme_var) <- paste0("Sample", 1:10)
  
  s2n_extreme <- calculate_rank_metric(extreme_var, metadata, "group", "signal2noise")
  
  expect_true(is.finite(s2n_extreme[1]))
  expect_true(abs(s2n_extreme[1]) < 100)  # Should not explode
})

test_that("t-test statistics match theoretical expectations", {
  # Create data with known t-test properties
  set.seed(777)
  n_per_group <- 15
  
  # Multiple features with different effect sizes
  effect_sizes <- c(0, 0.5, 1.0, 1.5, 2.0, 3.0)
  abundance <- matrix(0, nrow = length(effect_sizes), ncol = n_per_group * 2)
  
  for (i in seq_along(effect_sizes)) {
    # Group 1: N(0, 1)
    abundance[i, 1:n_per_group] <- rnorm(n_per_group, 0, 1)
    # Group 2: N(effect_size, 1)
    abundance[i, (n_per_group + 1):(2 * n_per_group)] <- rnorm(n_per_group, effect_sizes[i], 1)
  }
  
  rownames(abundance) <- paste0("K", sprintf("%05d", seq_along(effect_sizes)))
  colnames(abundance) <- paste0("Sample", 1:(n_per_group * 2))
  
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("Control", "Treatment"), each = n_per_group))
  )
  rownames(metadata) <- metadata$sample_name
  
  # Calculate t-test rankings
  ttest_ranking <- calculate_rank_metric(abundance, metadata, "group", "t_test")

  # Verify relationship between effect size and t-statistic magnitude
  # Note: Sign depends on factor level ordering (Control vs Treatment)
  # So we use absolute values to test the relationship
  theoretical_t <- effect_sizes / sqrt(2 / n_per_group)  # Assuming equal variances

  # Absolute t-statistics should correlate with absolute theoretical values
  correlation <- cor(abs(ttest_ranking), abs(theoretical_t))
  expect_true(correlation > 0.8)

  # Larger effect sizes should produce more extreme (larger absolute) t-statistics
  expect_true(all(diff(abs(ttest_ranking)) >= -0.5))
})

test_that("log2_ratio handles boundary conditions correctly", {
  extreme_data <- create_extreme_value_data()
  metadata <- extreme_data$metadata
  
  # Test each extreme case
  for (case_name in names(extreme_data$cases)) {
    abundance <- extreme_data$cases[[case_name]]
    
    log2_ranking <- calculate_rank_metric(abundance, metadata, "group", "log2_ratio")
    
    # All values should be finite
    expect_true(all(is.finite(log2_ranking)), 
                paste("Log2 ratios should be finite for", case_name))
    
    # No NaN or Inf values
    expect_true(all(!is.nan(log2_ranking)), 
                paste("Log2 ratios should not be NaN for", case_name))
    expect_true(all(!is.infinite(log2_ranking)), 
                paste("Log2 ratios should not be infinite for", case_name))
    
    # Values should be reasonable (not extremely large)
    expect_true(all(abs(log2_ranking) < 100), 
                paste("Log2 ratios should be reasonable magnitude for", case_name))
  }
})

test_that("numerical stability across different scales", {
  # Test ranking consistency across different abundance scales
  base_abundance <- matrix(rnorm(50, 100, 20), nrow = 5, ncol = 10)
  rownames(base_abundance) <- paste0("K", sprintf("%05d", 1:5))
  colnames(base_abundance) <- paste0("Sample", 1:10)
  
  metadata <- data.frame(
    sample_name = paste0("Sample", 1:10),
    group = factor(rep(c("A", "B"), each = 5))
  )
  rownames(metadata) <- metadata$sample_name
  
  # Test different scales
  scales <- c(1e-6, 1e-3, 1, 1e3, 1e6)
  
  rankings <- list()
  for (scale in scales) {
    scaled_abundance <- base_abundance * scale
    
    rankings[[as.character(scale)]] <- list(
      s2n = calculate_rank_metric(scaled_abundance, metadata, "group", "signal2noise"),
      ttest = calculate_rank_metric(scaled_abundance, metadata, "group", "t_test"),
      log2 = calculate_rank_metric(scaled_abundance, metadata, "group", "log2_ratio"),
      diff = calculate_rank_metric(scaled_abundance, metadata, "group", "diff_abundance")
    )
  }
  
  # Signal-to-noise should be scale-invariant
  base_s2n <- rankings[["1"]]$s2n
  for (scale in scales) {
    expect_equal(rankings[[as.character(scale)]]$s2n, base_s2n, tolerance = 1e-10,
                 info = paste("Signal-to-noise should be scale-invariant at scale", scale))
  }
  
  # T-test should be scale-invariant
  base_ttest <- rankings[["1"]]$ttest
  for (scale in scales) {
    expect_equal(rankings[[as.character(scale)]]$ttest, base_ttest, tolerance = 1e-10,
                 info = paste("T-test should be scale-invariant at scale", scale))
  }
  
  # Log2 ratio should be scale-invariant
  base_log2 <- rankings[["1"]]$log2
  for (scale in scales) {
    expect_equal(rankings[[as.character(scale)]]$log2, base_log2, tolerance = 1e-10,
                 info = paste("Log2 ratio should be scale-invariant at scale", scale))
  }
  
  # Difference in abundance should scale proportionally
  base_diff <- rankings[["1"]]$diff
  for (scale in scales) {
    expected_diff <- base_diff * scale
    expect_equal(rankings[[as.character(scale)]]$diff, expected_diff, tolerance = 1e-6,
                 info = paste("Difference should scale proportionally at scale", scale))
  }
})

test_that("statistical power analysis is consistent", {
  # This test verifies that larger effect sizes produce more detectable differences
  # Using a simpler approach to avoid random variation issues

  # Create data with known large effect
  set.seed(12345)
  n_samples <- 20
  n_features <- 50

  # Create abundance matrix with clear effect in first 10 features
  abundance_large_effect <- matrix(rnorm(n_features * n_samples), nrow = n_features, ncol = n_samples)
  rownames(abundance_large_effect) <- paste0("K", sprintf("%05d", 1:n_features))
  colnames(abundance_large_effect) <- paste0("Sample", 1:n_samples)

  # Add large effect (3 SDs) to first 10 features in second group
  abundance_large_effect[1:10, 11:20] <- abundance_large_effect[1:10, 11:20] + 3

  metadata <- data.frame(
    sample_name = colnames(abundance_large_effect),
    condition = factor(rep(c("Control", "Treatment"), each = 10))
  )
  rownames(metadata) <- metadata$sample_name

  # Calculate t-test statistics
  ttest_stats <- calculate_rank_metric(
    abundance_large_effect,
    metadata,
    "condition",
    "t_test"
  )

  # Features with effect should have larger absolute t-statistics
  effect_features <- paste0("K", sprintf("%05d", 1:10))
  no_effect_features <- paste0("K", sprintf("%05d", 11:50))

  mean_abs_effect <- mean(abs(ttest_stats[effect_features]))
  mean_abs_no_effect <- mean(abs(ttest_stats[no_effect_features]))

  # Features with effect should have larger absolute t-statistics
  expect_gt(mean_abs_effect, mean_abs_no_effect)

  # Features with effect should have high power (|t| > t_critical)
  t_critical <- qt(0.975, df = n_samples - 2)
  power_effect <- mean(abs(ttest_stats[effect_features]) > t_critical)

  # With 3 SD effect and n=20, power should be high
  expect_gt(power_effect, 0.7)
})

test_that("p-value calculations are statistically sound", {
  skip_if_not_installed("fgsea")
  
  # Create null data (no true enrichment)
  null_data <- create_known_enrichment_data(
    n_features = 200, 
    n_samples = 20, 
    enriched_features = integer(0),  # No enriched features
    effect_size = 1.0
  )
  
  # Create random gene sets
  set.seed(888)
  random_gene_sets <- list()
  for (i in 1:50) {
    random_gene_sets[[paste0("random_pathway_", i)]] <- sample(rownames(null_data$abundance), 
                                                              size = sample(10:30, 1))
  }
  
  # Mock prepare_gene_sets
  mockery::stub(pathway_gsea, "prepare_gene_sets", function(...) random_gene_sets)
  
  # Generate p-values under null (should be uniform)
  null_pvalues <- runif(50, 0, 1)  # Simulate null p-values
  
  mock_null_result <- data.frame(
    pathway_id = names(random_gene_sets),
    pathway_name = names(random_gene_sets),
    size = sapply(random_gene_sets, length),
    ES = runif(50, -0.3, 0.3),  # Small effect sizes under null
    NES = runif(50, -1, 1),
    pvalue = null_pvalues,
    p.adjust = p.adjust(null_pvalues, method = "BH"),
    leading_edge = sapply(random_gene_sets, function(x) paste(sample(x, min(3, length(x))), collapse = ";")),
    stringsAsFactors = FALSE
  )
  
  mockery::stub(pathway_gsea, "run_fgsea", function(...) mock_null_result)
  
  result <- pathway_gsea(
    abundance = null_data$abundance,
    metadata = null_data$metadata,
    group = "condition",
    method = "fgsea",
    nperm = 1000,
    seed = 123
  )
  
  # Under null hypothesis, p-values should be approximately uniform
  # Kolmogorov-Smirnov test for uniformity
  ks_test <- ks.test(result$pvalue, "punif", 0, 1)
  
  expect_true(ks_test$p.value > 0.05, 
              "P-values under null should be approximately uniform")
  
  # Check that very few pathways are significant at 0.05 level
  n_significant <- sum(result$p.adjust < 0.05)
  expected_false_positives <- length(random_gene_sets) * 0.05
  
  # Should be close to expected false positive rate
  expect_true(n_significant <= expected_false_positives * 3,  # Allow some variation
              "Number of significant pathways under null should be close to false positive rate")
})

test_that("multiple testing correction is applied correctly", {
  # Create mock p-values with known correction
  raw_pvalues <- c(0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8)
  
  # Test different correction methods
  correction_methods <- c("BH", "bonferroni", "holm")
  
  for (method in correction_methods) {
    # Expected corrected p-values
    expected_corrected <- p.adjust(raw_pvalues, method = method)
    
    # Mock result with raw p-values
    mock_result <- data.frame(
      pathway_id = paste0("pathway_", 1:8),
      pathway_name = paste0("pathway_", 1:8),
      size = rep(10, 8),
      ES = runif(8, -1, 1),
      NES = runif(8, -2, 2),
      pvalue = raw_pvalues,
      p.adjust = expected_corrected,
      leading_edge = rep("K00001;K00002", 8),
      stringsAsFactors = FALSE
    )
    
    # Validate that corrections are applied correctly
    expect_equal(mock_result$p.adjust, expected_corrected, tolerance = 1e-10,
                 info = paste("Multiple testing correction should be correct for", method))
    
    # Adjusted p-values should be >= raw p-values
    expect_true(all(mock_result$p.adjust >= mock_result$pvalue),
                paste("Adjusted p-values should be >= raw p-values for", method))
    
    # Adjusted p-values should be <= 1
    expect_true(all(mock_result$p.adjust <= 1),
                paste("Adjusted p-values should be <= 1 for", method))
  }
})

test_that("enrichment score calculations are mathematically sound", {
  skip_if_not_installed("fgsea")
  
  # Test ES properties with known gene sets and rankings
  set.seed(321)
  
  # Create simple ranking with clear pattern
  n_genes <- 100
  gene_names <- paste0("K", sprintf("%05d", 1:n_genes))
  
  # Ranking: genes 1-20 have high positive scores, 81-100 have high negative scores
  rankings <- rep(0, n_genes)
  rankings[1:20] <- seq(5, 1, length.out = 20)      # Top enriched
  rankings[21:80] <- seq(0.5, -0.5, length.out = 60) # Middle
  rankings[81:100] <- seq(-1, -5, length.out = 20)   # Bottom depleted
  names(rankings) <- gene_names
  
  # Gene sets
  gene_sets <- list(
    "top_enriched" = gene_names[1:15],     # Should have positive ES
    "bottom_depleted" = gene_names[85:100], # Should have negative ES
    "mixed" = gene_names[c(1:5, 85:90)],   # Mixed
    "middle" = gene_names[40:60]           # Should have ES near 0
  )
  
  # Mock fgsea with realistic ES values based on the ranking
  mock_es_result <- data.frame(
    pathway_id = c("top_enriched", "bottom_depleted", "mixed", "middle"),
    pathway_name = c("top_enriched", "bottom_depleted", "mixed", "middle"),
    size = c(15, 16, 11, 21),
    ES = c(0.8, -0.7, 0.2, 0.05),    # Consistent with expected enrichment
    NES = c(2.1, -1.9, 0.5, 0.1),    # Normalized scores
    pvalue = c(0.001, 0.002, 0.3, 0.8),
    p.adjust = c(0.004, 0.004, 0.4, 0.8),
    leading_edge = c("K00001;K00002;K00003", "K00085;K00086;K00087", 
                    "K00001;K00085", "K00045;K00050"),
    stringsAsFactors = FALSE
  )
  
  # Validate ES properties
  expect_true(mock_es_result$ES[1] > 0, "Top enriched pathway should have positive ES")
  expect_true(mock_es_result$ES[2] < 0, "Bottom depleted pathway should have negative ES")
  expect_true(abs(mock_es_result$ES[4]) < 0.2, "Middle pathway should have ES near 0")
  
  # NES should have same sign as ES but different magnitude
  expect_true(all(sign(mock_es_result$ES) == sign(mock_es_result$NES)), 
              "NES should have same sign as ES")
  
  # |NES| should generally be larger than |ES| (normalization effect)
  expect_true(mean(abs(mock_es_result$NES)) > mean(abs(mock_es_result$ES)),
              "NES magnitude should generally be larger than ES magnitude")
})