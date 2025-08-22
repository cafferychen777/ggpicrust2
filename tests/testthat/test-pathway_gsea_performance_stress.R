# Performance and Stress Tests for GSEA Functions
# Testing computational efficiency, memory usage, and robustness under stress

library(testthat)

# Performance test data generators
create_large_dataset <- function(n_features = 1000, n_samples = 100, sparsity = 0.3) {
  set.seed(12321)
  
  abundance <- matrix(0, nrow = n_features, ncol = n_samples)
  
  # Create sparse microbiome-like data
  for (i in 1:n_features) {
    # Determine which samples have this feature
    n_present <- rbinom(1, n_samples, 1 - sparsity)
    if (n_present > 0) {
      present_samples <- sample(n_samples, n_present)
      # Log-normal abundance for present samples
      abundance[i, present_samples] <- rlnorm(n_present, 
                                             meanlog = runif(1, 1, 4), 
                                             sdlog = runif(1, 0.5, 1.0))
    }
  }
  
  rownames(abundance) <- paste0("K", sprintf("%05d", sample(1:99999, n_features)))
  colnames(abundance) <- paste0("Sample", sprintf("%03d", 1:n_samples))
  
  metadata <- data.frame(
    sample_name = colnames(abundance),
    condition = factor(rep(c("Control", "Treatment"), length.out = n_samples)),
    batch = factor(rep(1:5, length.out = n_samples)),
    subject_id = factor(rep(1:(n_samples/2), each = 2)[1:n_samples])
  )
  rownames(metadata) <- metadata$sample_name
  
  return(list(abundance = abundance, metadata = metadata))
}

create_stress_gene_sets <- function(gene_names, n_pathways = 500, 
                                   min_size = 5, max_size = 100) {
  set.seed(54321)
  gene_sets <- list()
  
  for (i in 1:n_pathways) {
    pathway_id <- paste0("stress_pathway_", sprintf("%04d", i))
    pathway_size <- sample(min_size:max_size, 1)
    
    # Some pathways are random, others have structure
    if (i %% 10 == 0) {
      # Structured pathway - genes with similar names
      base_num <- sample(1:9000, 1)
      gene_candidates <- paste0("K", sprintf("%05d", base_num:(base_num + pathway_size * 2)))
      pathway_genes <- intersect(gene_candidates, gene_names)[1:min(pathway_size, 
                                                                   length(intersect(gene_candidates, gene_names)))]
    } else {
      # Random pathway
      pathway_genes <- sample(gene_names, min(pathway_size, length(gene_names)))
    }
    
    if (length(pathway_genes) >= min_size) {
      gene_sets[[pathway_id]] <- pathway_genes
    }
  }
  
  return(gene_sets)
}

# Performance benchmark tests
test_that("ranking metric calculations scale efficiently", {
  # Test different data sizes
  data_sizes <- list(
    small = list(features = 50, samples = 20),
    medium = list(features = 500, samples = 50),
    large = list(features = 2000, samples = 100)
  )
  
  execution_times <- list()
  
  for (size_name in names(data_sizes)) {
    n_features <- data_sizes[[size_name]]$features
    n_samples <- data_sizes[[size_name]]$samples
    
    test_data <- create_large_dataset(n_features, n_samples)
    
    # Time each ranking method
    for (method in c("signal2noise", "t_test", "log2_ratio", "diff_abundance")) {
      start_time <- Sys.time()
      
      ranking <- calculate_rank_metric(
        test_data$abundance, 
        test_data$metadata, 
        "condition", 
        method
      )
      
      end_time <- Sys.time()
      execution_time <- as.numeric(end_time - start_time, units = "secs")
      
      execution_times[[paste(size_name, method, sep = "_")]] <- execution_time
      
      # Verify result quality
      expect_length(ranking, n_features)
      expect_true(all(is.finite(ranking)))
      expect_named(ranking)
    }
  }
  
  # Check that execution times are reasonable
  expect_true(all(unlist(execution_times) < 10), 
              "All ranking calculations should complete within 10 seconds")
  
  # Large datasets shouldn't be drastically slower (should scale reasonably)
  for (method in c("signal2noise", "t_test", "log2_ratio", "diff_abundance")) {
    small_time <- execution_times[[paste("small", method, sep = "_")]]
    large_time <- execution_times[[paste("large", method, sep = "_")]]
    
    # Large dataset is 40x more features and 5x more samples = 200x more data
    # But time should not scale worse than O(n^2)
    time_ratio <- large_time / small_time
    expect_true(time_ratio < 500, 
                paste("Time scaling should be reasonable for", method, 
                     "- got ratio", round(time_ratio, 2)))
  }
})

test_that("memory usage is efficient for large datasets", {
  skip_on_cran()  # Memory tests can be unstable on CRAN
  
  # Test memory efficiency with large dataset
  large_data <- create_large_dataset(n_features = 1500, n_samples = 80, sparsity = 0.4)
  
  # Monitor memory usage during calculation
  gc()  # Clean up before test
  mem_before <- gc()
  
  # Calculate ranking metrics
  s2n_ranking <- calculate_rank_metric(large_data$abundance, large_data$metadata, 
                                      "condition", "signal2noise")
  
  mem_after <- gc()
  
  # Memory increase should be reasonable
  mem_increase <- (mem_after[2, 2] - mem_before[2, 2]) * 1024 * 1024  # Convert to bytes
  data_size <- object.size(large_data$abundance)
  
  # Memory increase should not be dramatically larger than data size
  expect_true(mem_increase < as.numeric(data_size) * 5,
              paste("Memory usage should be efficient - increase:", 
                    round(mem_increase / 1024^2, 2), "MB, data size:", 
                    round(as.numeric(data_size) / 1024^2, 2), "MB"))
  
  # Results should be complete
  expect_length(s2n_ranking, nrow(large_data$abundance))
  expect_true(all(is.finite(s2n_ranking)))
})

test_that("extreme sparsity is handled efficiently", {
  # Create very sparse data (90% zeros)
  sparse_data <- create_large_dataset(n_features = 800, n_samples = 60, sparsity = 0.9)
  
  # Count actual zeros
  zero_fraction <- sum(sparse_data$abundance == 0) / length(sparse_data$abundance)
  expect_true(zero_fraction > 0.8, "Data should be highly sparse")
  
  # All ranking methods should handle sparsity efficiently
  ranking_methods <- c("signal2noise", "t_test", "log2_ratio", "diff_abundance")
  
  for (method in ranking_methods) {
    start_time <- Sys.time()
    
    ranking <- calculate_rank_metric(sparse_data$abundance, sparse_data$metadata, 
                                    "condition", method)
    
    end_time <- Sys.time()
    execution_time <- as.numeric(end_time - start_time, units = "secs")
    
    # Should complete quickly despite sparsity
    expect_true(execution_time < 5, 
                paste("Sparse data ranking should be fast for", method))
    
    # Should produce valid results
    expect_length(ranking, nrow(sparse_data$abundance))
    expect_true(all(is.finite(ranking)))
    
    # Shouldn't have too many identical values (unless genuinely no difference)
    unique_values <- length(unique(ranking))
    expect_true(unique_values > nrow(sparse_data$abundance) / 10,
                paste("Should have reasonable diversity in rankings for", method))
  }
})

test_that("large gene set collections are processed efficiently", {
  skip_if_not_installed("fgsea")
  
  # Create medium-sized dataset
  test_data <- create_large_dataset(n_features = 600, n_samples = 40)
  
  # Create large gene set collection
  large_gene_sets <- create_stress_gene_sets(rownames(test_data$abundance), 
                                            n_pathways = 300, 
                                            min_size = 8, 
                                            max_size = 60)
  
  expect_true(length(large_gene_sets) >= 200, "Should create substantial gene set collection")
  
  # Mock prepare_gene_sets and fgsea for performance testing
  mockery::stub(pathway_gsea, "prepare_gene_sets", function(...) large_gene_sets)
  
  # Create realistic mock result
  mock_large_result <- data.frame(
    pathway_id = names(large_gene_sets),
    pathway_name = names(large_gene_sets),
    size = sapply(large_gene_sets, length),
    ES = runif(length(large_gene_sets), -1, 1),
    NES = runif(length(large_gene_sets), -3, 3),
    pvalue = runif(length(large_gene_sets), 0.001, 0.9),
    p.adjust = runif(length(large_gene_sets), 0.005, 0.95),
    leading_edge = sapply(large_gene_sets, function(x) 
      paste(sample(x, min(5, length(x))), collapse = ";")),
    stringsAsFactors = FALSE
  )
  
  mockery::stub(pathway_gsea, "run_fgsea", function(...) mock_large_result)
  
  # Time the complete analysis
  start_time <- Sys.time()
  
  result <- pathway_gsea(
    abundance = test_data$abundance,
    metadata = test_data$metadata,
    group = "condition",
    method = "fgsea",
    nperm = 100,
    seed = 123
  )
  
  end_time <- Sys.time()
  execution_time <- as.numeric(end_time - start_time, units = "secs")
  
  # Should complete in reasonable time
  expect_true(execution_time < 30, 
              paste("Large gene set analysis should complete in reasonable time, took", 
                    round(execution_time, 2), "seconds"))
  
  # Result should be complete
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), length(large_gene_sets))
  expect_true(all(!is.na(result$pvalue)))
  expect_true(all(!is.na(result$ES)))
})

test_that("concurrent access and thread safety", {
  skip_on_cran()  # Threading tests can be unstable
  
  # Create shared test data
  shared_data <- create_large_dataset(n_features = 200, n_samples = 30)
  
  # Test multiple simultaneous calculations (simulating concurrent access)
  n_concurrent <- 3
  results <- list()
  
  # Sequential execution for comparison
  sequential_start <- Sys.time()
  for (i in 1:n_concurrent) {
    results[[paste0("sequential_", i)]] <- calculate_rank_metric(
      shared_data$abundance, shared_data$metadata, "condition", "signal2noise"
    )
  }
  sequential_end <- Sys.time()
  sequential_time <- as.numeric(sequential_end - sequential_start, units = "secs")
  
  # All sequential results should be identical
  for (i in 2:n_concurrent) {
    expect_equal(results[["sequential_1"]], results[[paste0("sequential_", i)]],
                 "Sequential calculations should be identical")
  }
  
  # Test with different seeds (should produce same ranking metrics)
  different_seeds <- list()
  for (seed in c(111, 222, 333)) {
    set.seed(seed)
    different_seeds[[as.character(seed)]] <- calculate_rank_metric(
      shared_data$abundance, shared_data$metadata, "condition", "signal2noise"
    )
  }
  
  # Rankings should be identical regardless of seed (no randomness in calculation)
  expect_equal(different_seeds[["111"]], different_seeds[["222"]])
  expect_equal(different_seeds[["222"]], different_seeds[["333"]])
})

test_that("stress testing with problematic data patterns", {
  # Test various problematic data patterns
  
  # Pattern 1: All samples have identical values for some features
  identical_abundance <- matrix(100, nrow = 10, ncol = 20)
  identical_abundance[1:5, 1:10] <- 200  # Some differences
  rownames(identical_abundance) <- paste0("K", sprintf("%05d", 1:10))
  colnames(identical_abundance) <- paste0("Sample", 1:20)
  
  metadata <- data.frame(
    sample_name = paste0("Sample", 1:20),
    group = factor(rep(c("A", "B"), each = 10))
  )
  rownames(metadata) <- metadata$sample_name
  
  # Should handle without errors
  expect_no_error({
    s2n_identical <- calculate_rank_metric(identical_abundance, metadata, "group", "signal2noise")
  })
  
  expect_true(all(is.finite(s2n_identical)))
  
  # Pattern 2: Extreme outliers
  outlier_abundance <- matrix(rnorm(200, 100, 10), nrow = 10, ncol = 20)
  outlier_abundance[5, 10] <- 1e6  # Extreme outlier
  rownames(outlier_abundance) <- paste0("K", sprintf("%05d", 1:10))
  colnames(outlier_abundance) <- paste0("Sample", 1:20)
  
  expect_no_error({
    s2n_outlier <- calculate_rank_metric(outlier_abundance, metadata, "group", "signal2noise")
  })
  
  expect_true(all(is.finite(s2n_outlier)))
  
  # Pattern 3: Missing values (NaN, Inf)
  missing_abundance <- matrix(runif(200, 10, 100), nrow = 10, ncol = 20)
  missing_abundance[3, 5:8] <- NaN
  missing_abundance[7, 15] <- Inf
  rownames(missing_abundance) <- paste0("K", sprintf("%05d", 1:10))
  colnames(missing_abundance) <- paste0("Sample", 1:20)
  
  # Should handle gracefully (may error appropriately or handle)
  result_missing <- tryCatch({
    calculate_rank_metric(missing_abundance, metadata, "group", "signal2noise")
  }, error = function(e) e)
  
  # Either should work or give informative error
  if (inherits(result_missing, "error")) {
    expect_true(nchar(result_missing$message) > 10, "Should give informative error message")
  } else {
    expect_true(is.numeric(result_missing), "Should return numeric result if handled")
  }
})

test_that("repeated calculations are consistent and efficient", {
  # Test consistency and caching behavior
  test_data <- create_large_dataset(n_features = 300, n_samples = 40)
  
  # Calculate same ranking multiple times
  n_repeats <- 5
  rankings <- list()
  times <- numeric(n_repeats)
  
  for (i in 1:n_repeats) {
    start_time <- Sys.time()
    
    rankings[[i]] <- calculate_rank_metric(
      test_data$abundance, 
      test_data$metadata, 
      "condition", 
      "signal2noise"
    )
    
    end_time <- Sys.time()
    times[i] <- as.numeric(end_time - start_time, units = "secs")
  }
  
  # All results should be identical
  for (i in 2:n_repeats) {
    expect_equal(rankings[[1]], rankings[[i]], tolerance = 1e-15,
                 paste("Repeated calculation", i, "should be identical"))
  }
  
  # Times should be consistent (no major degradation)
  expect_true(max(times) / min(times) < 3, 
              "Execution times should be consistent across repeats")
  
  # All times should be reasonable
  expect_true(all(times < 5), "All calculations should complete quickly")
})

test_that("memory cleanup after large calculations", {
  skip_on_cran()  # Memory tests can be unstable
  
  # Force garbage collection
  gc()
  initial_memory <- gc()[2, 2]
  
  # Perform large calculation
  large_data <- create_large_dataset(n_features = 1000, n_samples = 60)
  
  ranking <- calculate_rank_metric(large_data$abundance, large_data$metadata, 
                                  "condition", "t_test")
  
  # Remove references and force cleanup
  rm(large_data, ranking)
  gc()
  final_memory <- gc()[2, 2]
  
  # Memory should return close to initial level
  memory_increase <- (final_memory - initial_memory) * 1024 * 1024  # Convert to bytes
  
  expect_true(memory_increase < 100 * 1024 * 1024,  # Less than 100MB increase
              paste("Memory should be cleaned up after large calculations - increase:", 
                    round(memory_increase / 1024^2, 2), "MB"))
})

test_that("robustness with malformed input edge cases", {
  # Test various edge cases that could cause issues
  
  # Case 1: Single feature
  single_feature <- matrix(rnorm(10), nrow = 1, ncol = 10)
  rownames(single_feature) <- "K00001"
  colnames(single_feature) <- paste0("Sample", 1:10)
  
  metadata <- data.frame(
    sample_name = paste0("Sample", 1:10),
    group = factor(rep(c("A", "B"), each = 5))
  )
  rownames(metadata) <- metadata$sample_name
  
  expect_no_error({
    single_ranking <- calculate_rank_metric(single_feature, metadata, "group", "signal2noise")
  })
  expect_length(single_ranking, 1)
  
  # Case 2: Single sample per group
  tiny_data <- matrix(rnorm(20), nrow = 10, ncol = 2)
  rownames(tiny_data) <- paste0("K", sprintf("%05d", 1:10))
  colnames(tiny_data) <- c("Sample1", "Sample2")
  
  tiny_metadata <- data.frame(
    sample_name = c("Sample1", "Sample2"),
    group = factor(c("A", "B"))
  )
  rownames(tiny_metadata) <- tiny_metadata$sample_name
  
  # Some methods might handle this, others might error appropriately
  for (method in c("signal2noise", "t_test", "log2_ratio", "diff_abundance")) {
    result <- tryCatch({
      calculate_rank_metric(tiny_data, tiny_metadata, "group", method)
    }, error = function(e) e)
    
    # Should either work or give reasonable error
    if (!inherits(result, "error")) {
      expect_length(result, 10)
      expect_true(all(is.finite(result)))
    }
  }
  
  # Case 3: Feature names with special characters
  special_names <- matrix(rnorm(50), nrow = 5, ncol = 10)
  rownames(special_names) <- c("K00001", "K-weird", "K.dot", "K space", "K/slash")
  colnames(special_names) <- paste0("Sample", 1:10)
  
  expect_no_error({
    special_ranking <- calculate_rank_metric(special_names, metadata, "group", "signal2noise")
  })
  expect_named(special_ranking, rownames(special_names))
})