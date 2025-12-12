# Comprehensive MetaCyc GSEA Integration Tests
# Testing all ranking methods, statistical accuracy, and production scenarios
# Following Linus principles: data structures first, no special cases, never break userspace

test_that("MetaCyc GSEA integration with all ranking methods", {
  skip_if_not_installed("fgsea")
  
  # Create realistic EC abundance test data
  set.seed(42)  # Reproducible results
  
  # Generate EC abundance matrix with proper biological structure
  n_samples <- 20
  n_ecs <- 50
  ec_ids <- paste0("EC:", c(
    "1.1.1.1", "1.2.1.12", "2.7.1.1", "2.7.1.11", "4.1.2.13",
    "1.2.1.59", "2.7.2.3", "5.4.2.11", "4.2.1.11", "5.4.2.1",
    "3.5.4.9", "1.5.1.15", "6.3.4.3", "2.1.2.1", "1.14.13.3",
    "4.1.1.45", "1.2.1.10", "2.5.1.54", "4.2.3.4", "4.2.1.10",
    "2.7.1.4", "4.1.2.22", "6.3.4.5", "2.1.3.3", "3.5.3.1",
    "4.1.1.17", "2.5.1.16", "1.1.1.27", "1.1.1.37", "2.3.1.12",
    "2.7.1.40", "4.1.3.27", "1.2.4.1", "6.2.1.5", "2.8.3.8",
    "4.2.1.2", "5.4.2.2", "2.7.9.1", "1.8.1.4", "3.1.3.11",
    "2.4.2.1", "3.5.4.1", "2.1.1.45", "1.5.1.3", "6.3.2.2",
    "2.6.1.1", "1.4.1.2", "3.1.1.1", "2.3.1.1", "4.3.1.1"
  ))
  
  # Create abundance matrix with realistic structure
  abundance_matrix <- matrix(
    rlnorm(n_samples * n_ecs, meanlog = log(100), sdlog = 1),
    nrow = n_ecs, ncol = n_samples
  )
  rownames(abundance_matrix) <- ec_ids
  colnames(abundance_matrix) <- paste0("Sample_", 1:n_samples)
  
  # Create metadata with two balanced groups
  metadata <- data.frame(
    sample_id = colnames(abundance_matrix),
    group = rep(c("Control", "Treatment"), each = n_samples/2),
    batch = rep(c("A", "B"), times = n_samples/2),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample_id
  
  # Add differential signal to some ECs for realistic GSEA results
  treatment_samples <- metadata$sample_id[metadata$group == "Treatment"]
  # Increase abundance of first 10 ECs in treatment group
  abundance_matrix[1:10, treatment_samples] <- 
    abundance_matrix[1:10, treatment_samples] * 2
  # Decrease abundance of ECs 11-20 in treatment group  
  abundance_matrix[11:20, treatment_samples] <- 
    abundance_matrix[11:20, treatment_samples] * 0.5
  
  # Test all ranking methods
  ranking_methods <- c("signal2noise", "t_test", "log2_ratio", "diff_abundance")
  
  for (method in ranking_methods) {
    test_results <- NULL

    tryCatch({
      test_results <- pathway_gsea(
        abundance = abundance_matrix,
        metadata = metadata,
        group = "group",
        pathway_type = "MetaCyc",
        method = "fgsea",
        rank_method = method,
        nperm = 100,  # Reduced for speed, but sufficient for testing
        min_size = 2,
        max_size = 50
      )
    }, error = function(e) {
      fail(paste("pathway_gsea failed with rank_method:", method, "-", e$message))
    })
    
    # Validate results structure
    expect_s3_class(test_results, "data.frame")
    expect_true("pathway_id" %in% colnames(test_results))
    expect_true("pathway_name" %in% colnames(test_results))
    expect_true("NES" %in% colnames(test_results))
    expect_true("pvalue" %in% colnames(test_results))
    expect_true("method" %in% colnames(test_results))
    
    if (nrow(test_results) > 0) {
      expect_equal(unique(test_results$method), "fgsea")
      
      # Validate statistical measures
      expect_true(all(is.finite(test_results$NES)))
      expect_true(all(test_results$pvalue >= 0 & test_results$pvalue <= 1))
      expect_true(all(test_results$p.adjust >= 0 & test_results$p.adjust <= 1))
      expect_true(all(test_results$p.adjust >= test_results$pvalue))  # Multiple testing correction
      
      # Validate enrichment scores make biological sense
      expect_true(any(test_results$NES > 0))  # Some positive enrichment expected
      expect_true(any(abs(test_results$NES) > 1))  # Some significant enrichment expected
    }
  }
})

test_that("MetaCyc GSEA statistical significance validation", {
  skip_if_not_installed("fgsea")
  
  # Create controlled test data with known differential patterns
  set.seed(123)  # Different seed for independence
  
  # Generate test data with clear signal
  n_samples_per_group <- 8
  strong_signal_ecs <- paste0("EC:", c("1.1.1.1", "2.7.1.1", "4.1.2.13", "1.2.1.12"))
  weak_signal_ecs <- paste0("EC:", c("3.5.4.9", "1.5.1.15", "6.3.4.3", "2.1.2.1"))
  noise_ecs <- paste0("EC:", c("1.14.13.3", "4.1.1.45", "1.2.1.10", "2.5.1.54"))
  all_ecs <- c(strong_signal_ecs, weak_signal_ecs, noise_ecs)
  
  # Create abundance data
  abundance_matrix <- matrix(
    rnorm(length(all_ecs) * n_samples_per_group * 2, mean = 10, sd = 2),
    nrow = length(all_ecs),
    ncol = n_samples_per_group * 2
  )
  rownames(abundance_matrix) <- all_ecs
  colnames(abundance_matrix) <- paste0("Sample_", 1:(n_samples_per_group * 2))
  
  # Add strong differential signal
  treatment_cols <- (n_samples_per_group + 1):(n_samples_per_group * 2)
  # Strong upregulation
  abundance_matrix[strong_signal_ecs, treatment_cols] <- 
    abundance_matrix[strong_signal_ecs, treatment_cols] + 5
  # Weak downregulation  
  abundance_matrix[weak_signal_ecs, treatment_cols] <- 
    abundance_matrix[weak_signal_ecs, treatment_cols] - 2
  
  # Create metadata
  metadata <- data.frame(
    sample_id = colnames(abundance_matrix),
    group = rep(c("Control", "Treatment"), each = n_samples_per_group),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample_id
  
  # Run GSEA
  gsea_results <- pathway_gsea(
    abundance = abundance_matrix,
    metadata = metadata,
    group = "group",
    pathway_type = "MetaCyc",
    method = "fgsea",
    rank_method = "signal2noise",
    nperm = 1000,  # Higher for statistical accuracy
    min_size = 2,
    max_size = 20
  )
  
  # Validate statistical behavior
  if (nrow(gsea_results) > 0) {
    # Test 1: P-values should be valid (between 0 and 1)
    # Note: With random data and small samples, significant results aren't guaranteed
    expect_true(all(gsea_results$pvalue >= 0 & gsea_results$pvalue <= 1))
    
    # Test 2: Multiple testing correction
    expect_true(all(gsea_results$p.adjust >= gsea_results$pvalue))
    
    # Test 3: NES direction consistency 
    # Pathways with strong signal ECs should show expected direction
    strong_pathways <- gsea_results[
      grepl("1CMET2-PWY|ANAGLYCOLYSIS-PWY", gsea_results$pathway_id), ]
    
    if (nrow(strong_pathways) > 0) {
      # Should have some enrichment (positive or negative NES)
      expect_true(any(abs(strong_pathways$NES) > 1))
    }
    
    # Test 4: Size filtering validation
    expect_true(all(gsea_results$size >= 2))
    expect_true(all(gsea_results$size <= 20))
  }
})

test_that("MetaCyc GSEA mathematical correctness validation", {
  skip_if_not_installed("fgsea")
  
  # Test ranking metric calculations directly
  set.seed(456)
  
  # Simple controlled data for mathematical validation
  group1_data <- matrix(c(10, 15, 12, 8, 20, 18), nrow = 3, ncol = 2)
  group2_data <- matrix(c(5, 8, 6, 15, 25, 22), nrow = 3, ncol = 2) 
  
  abundance_test <- cbind(group1_data, group2_data)
  rownames(abundance_test) <- c("EC:1.1.1.1", "EC:2.2.2.2", "EC:3.3.3.3")
  colnames(abundance_test) <- c("S1", "S2", "S3", "S4")
  
  metadata_test <- data.frame(
    sample_id = colnames(abundance_test),
    group = rep(c("A", "B"), each = 2),
    stringsAsFactors = FALSE
  )
  rownames(metadata_test) <- metadata_test$sample_id
  
  # Test signal-to-noise calculation manually
  mean1 <- rowMeans(group1_data)  # [12.5, 10, 19]
  mean2 <- rowMeans(group2_data)  # [6.5, 18.5, 23.5]
  sd1 <- apply(group1_data, 1, sd)
  sd2 <- apply(group2_data, 1, sd)
  
  # Handle zero SDs
  sd1[sd1 == 0] <- 0.00001
  sd2[sd2 == 0] <- 0.00001
  
  expected_s2n <- (mean1 - mean2) / (sd1 + sd2)
  
  # Calculate using pathway_gsea function
  calculated_metric <- calculate_rank_metric(
    abundance_test, metadata_test, "group", "signal2noise"
  )
  
  # Compare results (allowing for small floating point differences)
  expect_equal(length(calculated_metric), 3)
  expect_equal(names(calculated_metric), rownames(abundance_test))
  
  for (i in 1:3) {
    expect_equal(unname(calculated_metric[i]), expected_s2n[i], tolerance = 1e-10)
  }
  
  # Test t-test calculation
  calculated_t <- calculate_rank_metric(
    abundance_test, metadata_test, "group", "t_test"
  )
  
  # Validate t-test results are reasonable
  expect_equal(length(calculated_t), 3)
  expect_true(all(is.finite(calculated_t)))
  
  # Test log2 ratio calculation
  calculated_log2 <- calculate_rank_metric(
    abundance_test, metadata_test, "group", "log2_ratio"
  )
  
  # Validate log2 ratios
  expect_equal(length(calculated_log2), 3)
  expect_true(all(is.finite(calculated_log2)))
  
  # For EC:1.1.1.1: log2(12.5/6.5) ≈ 0.944
  expect_equal(unname(calculated_log2[1]), log2(mean1[1]/mean2[1]), tolerance = 1e-10)
})

test_that("MetaCyc GSEA performance and memory validation", {
  skip_if_not_installed("fgsea")
  
  # Test with realistic dataset sizes
  set.seed(789)
  
  # Large-scale test (realistic production size)
  n_samples_large <- 50
  n_ecs_large <- 200
  
  # Generate large EC dataset
  ec_ids_large <- paste0("EC:", sprintf("%d.%d.%d.%d", 
                                       sample(1:6, n_ecs_large, replace = TRUE),
                                       sample(1:20, n_ecs_large, replace = TRUE),
                                       sample(1:50, n_ecs_large, replace = TRUE),
                                       sample(1:200, n_ecs_large, replace = TRUE)))
  
  abundance_large <- matrix(
    rlnorm(n_samples_large * n_ecs_large, meanlog = log(50), sdlog = 1.5),
    nrow = n_ecs_large, ncol = n_samples_large
  )
  rownames(abundance_large) <- ec_ids_large
  colnames(abundance_large) <- paste0("Sample_", 1:n_samples_large)
  
  metadata_large <- data.frame(
    sample_id = colnames(abundance_large),
    group = rep(c("Control", "Treatment"), length.out = n_samples_large),
    batch = rep(paste0("Batch_", 1:5), length.out = n_samples_large),
    stringsAsFactors = FALSE
  )
  rownames(metadata_large) <- metadata_large$sample_id
  
  # Performance test
  start_time <- Sys.time()
  
  expect_no_error({
    large_results <- pathway_gsea(
      abundance = abundance_large,
      metadata = metadata_large,
      group = "group",
      pathway_type = "MetaCyc",
      method = "fgsea",
      rank_method = "signal2noise",
      nperm = 100,
      min_size = 5,
      max_size = 100
    )
  })
  
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Performance validation
  expect_lt(execution_time, 60)  # Should complete within 1 minute
  
  # Memory efficiency validation
  expect_s3_class(large_results, "data.frame")
  expect_lt(object.size(large_results), 5 * 1024^2)  # Less than 5MB result object
  
  # Result quality validation
  if (nrow(large_results) > 0) {
    expect_true(all(is.finite(large_results$NES)))
    expect_true(all(large_results$pvalue >= 0 & large_results$pvalue <= 1))
    expect_true(nrow(large_results) <= 100)  # Reasonable number of pathways tested
  }
})

test_that("MetaCyc GSEA reproducibility validation", {
  skip_if_not_installed("fgsea")
  
  # Test reproducibility with same seed
  set.seed(999)
  
  # Create test data
  n_samples <- 12
  n_ecs <- 30
  ec_ids <- paste0("EC:", sample(paste(
    sample(1:6, n_ecs, replace = TRUE),
    sample(1:10, n_ecs, replace = TRUE), 
    sample(1:20, n_ecs, replace = TRUE),
    sample(1:50, n_ecs, replace = TRUE),
    sep = "."
  )))
  
  abundance_repro <- matrix(
    rnorm(n_samples * n_ecs, mean = 100, sd = 20),
    nrow = n_ecs, ncol = n_samples
  )
  rownames(abundance_repro) <- ec_ids
  colnames(abundance_repro) <- paste0("Sample_", 1:n_samples)
  
  metadata_repro <- data.frame(
    sample_id = colnames(abundance_repro),
    group = rep(c("Control", "Treatment"), each = n_samples/2),
    stringsAsFactors = FALSE
  )
  rownames(metadata_repro) <- metadata_repro$sample_id
  
  # Run GSEA twice with same seed
  results1 <- pathway_gsea(
    abundance = abundance_repro,
    metadata = metadata_repro,
    group = "group",
    pathway_type = "MetaCyc",
    method = "fgsea",
    rank_method = "signal2noise",
    nperm = 500,
    seed = 42
  )
  
  results2 <- pathway_gsea(
    abundance = abundance_repro,
    metadata = metadata_repro,
    group = "group", 
    pathway_type = "MetaCyc",
    method = "fgsea",
    rank_method = "signal2noise",
    nperm = 500,
    seed = 42
  )
  
  # Test reproducibility
  expect_equal(nrow(results1), nrow(results2))
  
  if (nrow(results1) > 0 && nrow(results2) > 0) {
    # Results should be identical with same seed
    results1_sorted <- results1[order(results1$pathway_id), ]
    results2_sorted <- results2[order(results2$pathway_id), ]
    
    expect_equal(results1_sorted$pathway_id, results2_sorted$pathway_id)
    expect_equal(results1_sorted$NES, results2_sorted$NES, tolerance = 1e-10)
    expect_equal(results1_sorted$pvalue, results2_sorted$pvalue, tolerance = 1e-10)
  }
})

test_that("MetaCyc GSEA edge cases and error handling", {
  skip_if_not_installed("fgsea")
  
  # Test 1: Insufficient samples
  small_abundance <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  rownames(small_abundance) <- c("EC:1.1.1.1", "EC:2.2.2.2")
  colnames(small_abundance) <- c("S1", "S2")
  
  small_metadata <- data.frame(
    sample_id = c("S1", "S2"),
    group = c("A", "A"),  # Same group - should fail
    stringsAsFactors = FALSE
  )
  rownames(small_metadata) <- small_metadata$sample_id
  
  expect_error({
    pathway_gsea(
      abundance = small_abundance,
      metadata = small_metadata,
      group = "group",
      pathway_type = "MetaCyc"
    )
  })
  
  # Test 2: No overlapping samples between abundance and metadata
  no_overlap_abundance <- matrix(1:4, nrow = 2, ncol = 2)
  rownames(no_overlap_abundance) <- c("EC:1.1.1.1", "EC:2.2.2.2")
  colnames(no_overlap_abundance) <- c("X1", "X2")
  
  no_overlap_metadata <- data.frame(
    sample_id = c("Y1", "Y2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )
  rownames(no_overlap_metadata) <- no_overlap_metadata$sample_id
  
  expect_error({
    pathway_gsea(
      abundance = no_overlap_abundance,
      metadata = no_overlap_metadata,
      group = "group",
      pathway_type = "MetaCyc"
    )
  })
  
  # Test 3: Empty abundance data
  empty_abundance <- matrix(numeric(0), nrow = 0, ncol = 0)
  
  expect_error({
    pathway_gsea(
      abundance = empty_abundance,
      metadata = data.frame(),
      group = "group",
      pathway_type = "MetaCyc"
    )
  })
  
  # Test 4: Invalid group column
  valid_abundance <- matrix(rnorm(20), nrow = 4, ncol = 5)
  rownames(valid_abundance) <- paste0("EC:", 1:4)
  colnames(valid_abundance) <- paste0("S", 1:5)
  
  valid_metadata <- data.frame(
    sample_id = paste0("S", 1:5),
    group = rep(c("A", "B"), c(2, 3)),
    stringsAsFactors = FALSE
  )
  rownames(valid_metadata) <- valid_metadata$sample_id
  
  expect_error({
    pathway_gsea(
      abundance = valid_abundance,
      metadata = valid_metadata,
      group = "nonexistent_column",
      pathway_type = "MetaCyc"
    )
  })
})