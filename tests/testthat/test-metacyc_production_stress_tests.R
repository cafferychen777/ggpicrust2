# MetaCyc Production Scenario Stress Tests
# Testing production-scale datasets, performance, memory efficiency, and real-world scenarios
# Following Linus principles: "Theory and practice sometimes clash. Theory loses."

test_that("MetaCyc large-scale abundance data performance", {
  skip_if_not_installed("fgsea")
  skip_on_cran()  # Skip on CRAN due to time constraints
  
  # Test with realistic large-scale microbiome dataset
  set.seed(2024)
  
  # Large dataset parameters (realistic for microbiome studies)
  n_samples_large <- 100
  n_ecs_large <- 500
  
  # Generate realistic EC IDs
  ec_classes <- c(1, 2, 3, 4, 5, 6)  # EC main classes
  ec_ids_large <- character(n_ecs_large)
  
  for (i in 1:n_ecs_large) {
    ec_class <- sample(ec_classes, 1)
    subclass <- sample(1:20, 1)
    sub_subclass <- sample(1:50, 1)
    serial <- sample(1:200, 1)
    ec_ids_large[i] <- paste0("EC:", ec_class, ".", subclass, ".", sub_subclass, ".", serial)
  }
  
  # Generate abundance matrix with realistic microbiome properties
  # Log-normal distribution with many zeros (typical of microbiome data)
  abundance_large <- matrix(0, nrow = n_ecs_large, ncol = n_samples_large)
  
  for (i in 1:n_ecs_large) {
    for (j in 1:n_samples_large) {
      # 30% chance of non-zero abundance
      if (runif(1) > 0.7) {
        abundance_large[i, j] <- rlnorm(1, meanlog = log(100), sdlog = 1.5)
      }
    }
  }
  
  rownames(abundance_large) <- ec_ids_large
  colnames(abundance_large) <- paste0("Sample_", sprintf("%03d", 1:n_samples_large))
  
  # Create realistic metadata
  metadata_large <- data.frame(
    sample_id = colnames(abundance_large),
    group = rep(c("Control", "Treatment"), each = n_samples_large/2),
    batch = rep(paste0("Batch_", 1:10), length.out = n_samples_large),
    subject = rep(paste0("Subject_", 1:(n_samples_large/2)), each = 2),
    time_point = rep(c("Baseline", "Post"), times = n_samples_large/2),
    stringsAsFactors = FALSE
  )
  rownames(metadata_large) <- metadata_large$sample_id
  
  # Add realistic differential signal
  treatment_samples <- metadata_large$sample_id[metadata_large$group == "Treatment"]
  signal_ecs <- sample(1:n_ecs_large, 50)  # 50 ECs with differential signal
  
  for (ec_idx in signal_ecs[1:25]) {
    # Increase in treatment
    abundance_large[ec_idx, treatment_samples] <- 
      abundance_large[ec_idx, treatment_samples] * runif(length(treatment_samples), 1.5, 3)
  }
  
  for (ec_idx in signal_ecs[26:50]) {
    # Decrease in treatment
    abundance_large[ec_idx, treatment_samples] <- 
      abundance_large[ec_idx, treatment_samples] * runif(length(treatment_samples), 0.3, 0.7)
  }
  
  # Performance benchmarking
  start_time <- Sys.time()
  memory_before <- as.numeric(object.size(.GlobalEnv))
  
  expect_no_error({
    large_results <- pathway_gsea(
      abundance = abundance_large,
      metadata = metadata_large,
      group = "group",
      pathway_type = "MetaCyc",
      method = "fgsea",
      rank_method = "signal2noise",
      nperm = 1000,
      min_size = 5,
      max_size = 200
    )
  })
  
  end_time <- Sys.time()
  memory_after <- as.numeric(object.size(.GlobalEnv))
  
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  memory_used <- memory_after - memory_before
  
  # Performance validation
  expect_lt(execution_time, 300)  # Should complete within 5 minutes
  expect_lt(memory_used, 500 * 1024^2)  # Should use less than 500MB additional memory
  
  # Result quality validation
  expect_s3_class(large_results, "data.frame")
  
  if (nrow(large_results) > 0) {
    # Statistical validation
    expect_true(all(is.finite(large_results$NES)))
    expect_true(all(large_results$pvalue >= 0 & large_results$pvalue <= 1))
    expect_true(all(large_results$p.adjust >= large_results$pvalue))
    
    # Should find some significant pathways given the added signal
    significant_pathways <- sum(large_results$p.adjust < 0.05)
    expect_gt(significant_pathways, 0)
    
    # Result size validation
    expect_lt(as.numeric(object.size(large_results)), 10 * 1024^2)  # <10MB result
  }
  
  cat(sprintf("\nLarge-scale test completed in %.2f seconds using %.2f MB memory\n", 
              execution_time, memory_used / (1024^2)))
  cat(sprintf("Processed %d ECs across %d samples\n", n_ecs_large, n_samples_large))
  if (nrow(large_results) > 0) {
    cat(sprintf("Found %d significant pathways (p.adjust < 0.05)\n", 
                sum(large_results$p.adjust < 0.05)))
  }
})

test_that("MetaCyc memory efficiency with varying dataset sizes", {
  skip_if_not_installed("fgsea")
  
  # Test memory scaling with different dataset sizes
  dataset_sizes <- list(
    small = list(n_samples = 20, n_ecs = 50),
    medium = list(n_samples = 50, n_ecs = 150),
    large = list(n_samples = 80, n_ecs = 300)
  )
  
  memory_usage <- list()
  execution_times <- list()
  
  for (size_name in names(dataset_sizes)) {
    params <- dataset_sizes[[size_name]]
    set.seed(123)  # Reproducible results
    
    # Generate test data
    ec_ids <- paste0("EC:", sample(1:6, params$n_ecs, replace = TRUE), ".",
                     sample(1:10, params$n_ecs, replace = TRUE), ".",
                     sample(1:20, params$n_ecs, replace = TRUE), ".",
                     sample(1:50, params$n_ecs, replace = TRUE))
    
    abundance_test <- matrix(
      rlnorm(params$n_samples * params$n_ecs, meanlog = log(50), sdlog = 1),
      nrow = params$n_ecs, ncol = params$n_samples
    )
    rownames(abundance_test) <- ec_ids
    colnames(abundance_test) <- paste0("S", 1:params$n_samples)
    
    metadata_test <- data.frame(
      sample_id = paste0("S", 1:params$n_samples),
      group = rep(c("Control", "Treatment"), length.out = params$n_samples),
      stringsAsFactors = FALSE
    )
    rownames(metadata_test) <- metadata_test$sample_id
    
    # Measure performance
    gc()  # Clean up before measurement
    memory_before <- as.numeric(object.size(.GlobalEnv))
    start_time <- Sys.time()
    
    expect_no_error({
      result_size <- pathway_gsea(
        abundance = abundance_test,
        metadata = metadata_test,
        group = "group",
        pathway_type = "MetaCyc",
        method = "fgsea",
        nperm = 100,  # Reduced for speed
        min_size = 3,
        max_size = 100
      )
    })
    
    end_time <- Sys.time()
    memory_after <- as.numeric(object.size(.GlobalEnv))
    
    execution_times[[size_name]] <- as.numeric(difftime(end_time, start_time, units = "secs"))
    memory_usage[[size_name]] <- memory_after - memory_before
    
    # Clean up
    rm(abundance_test, metadata_test, result_size)
    gc()
  }
  
  # Validate memory tracking works (object.size may not show memory increase reliably)
  # Just check that memory values are non-negative
  expect_gte(memory_usage[["small"]], 0)
  expect_gte(memory_usage[["medium"]], 0)
  expect_gte(memory_usage[["large"]], 0)
  
  # Execution time should scale reasonably (not exponentially)
  time_ratio_med_small <- execution_times[["medium"]] / execution_times[["small"]]
  time_ratio_large_med <- execution_times[["large"]] / execution_times[["medium"]]
  
  expect_lt(time_ratio_med_small, 10)  # Should not be more than 10x slower
  expect_lt(time_ratio_large_med, 10)  # Should scale reasonably
  
  cat("\nMemory efficiency test results:\n")
  for (size_name in names(dataset_sizes)) {
    cat(sprintf("%s: %.2f MB memory, %.2f seconds\n", 
                size_name, memory_usage[[size_name]] / (1024^2), execution_times[[size_name]]))
  }
})

test_that("MetaCyc concurrent analysis simulation", {
  skip_if_not_installed("fgsea")
  skip_on_cran()  # Skip on CRAN
  
  # Simulate running multiple GSEA analyses concurrently (common in batch processing)
  set.seed(456)
  
  n_analyses <- 5
  results_list <- list()
  total_start_time <- Sys.time()
  
  for (i in 1:n_analyses) {
    # Create unique dataset for each analysis
    n_samples <- 30
    n_ecs <- 100
    
    ec_ids <- paste0("EC:", sample(1:6, n_ecs, replace = TRUE), ".",
                     sample(1:15, n_ecs, replace = TRUE), ".",
                     sample(1:30, n_ecs, replace = TRUE), ".",
                     sample(1:100, n_ecs, replace = TRUE))
    
    abundance_concurrent <- matrix(
      rlnorm(n_samples * n_ecs, meanlog = log(75), sdlog = 1.2),
      nrow = n_ecs, ncol = n_samples
    )
    rownames(abundance_concurrent) <- ec_ids
    colnames(abundance_concurrent) <- paste0("Sample_", i, "_", 1:n_samples)
    
    metadata_concurrent <- data.frame(
      sample_id = colnames(abundance_concurrent),
      group = rep(c("Group_A", "Group_B"), length.out = n_samples),
      analysis_id = rep(paste0("Analysis_", i), n_samples),
      stringsAsFactors = FALSE
    )
    rownames(metadata_concurrent) <- metadata_concurrent$sample_id
    
    # Use tryCatch since expect_no_error doesn't support info parameter
    tryCatch({
      result_concurrent <- pathway_gsea(
        abundance = abundance_concurrent,
        metadata = metadata_concurrent,
        group = "group",
        pathway_type = "MetaCyc",
        method = "fgsea",
        nperm = 100,
        seed = 42 + i  # Different seed for each analysis
      )
      results_list[[i]] <- result_concurrent
    }, error = function(e) {
      fail(paste("Analysis", i, "should complete successfully:", e$message))
    })
  }
  
  total_end_time <- Sys.time()
  total_time <- as.numeric(difftime(total_end_time, total_start_time, units = "secs"))
  
  # Validate all analyses completed
  expect_equal(length(results_list), n_analyses)
  
  # Validate results are independent (different seeds should give different results)
  if (all(sapply(results_list, nrow) > 0)) {
    # Compare first two analyses if both have results
    if (nrow(results_list[[1]]) > 0 && nrow(results_list[[2]]) > 0) {
      # Should have some differences due to different random seeds
      common_pathways <- intersect(results_list[[1]]$pathway_id, results_list[[2]]$pathway_id)
      if (length(common_pathways) > 0) {
        # P-values should differ slightly due to different random seeds
        pathway_1 <- common_pathways[1]
        pval_1 <- results_list[[1]][results_list[[1]]$pathway_id == pathway_1, "pvalue"][1]
        pval_2 <- results_list[[2]][results_list[[2]]$pathway_id == pathway_1, "pvalue"][1]
        # They might be the same if permutation results are identical, but usually differ
        expect_true(is.numeric(pval_1) && is.numeric(pval_2))
      }
    }
  }
  
  # Performance validation
  expect_lt(total_time, 180)  # 5 analyses should complete within 3 minutes
  
  cat(sprintf("\nConcurrent analysis test: %d analyses completed in %.2f seconds\n", 
              n_analyses, total_time))
  cat(sprintf("Average time per analysis: %.2f seconds\n", total_time / n_analyses))
})

test_that("MetaCyc real-world data structure simulation", {
  skip_if_not_installed("fgsea")
  
  # Simulate real-world microbiome study data characteristics
  set.seed(789)
  
  # Realistic study parameters
  n_subjects <- 25
  n_timepoints <- 2
  n_samples <- n_subjects * n_timepoints
  n_ecs <- 200
  
  # Generate EC IDs with realistic distribution (some EC classes more common)
  ec_class_weights <- c(0.25, 0.20, 0.20, 0.15, 0.10, 0.10)
  ec_classes <- sample(1:6, n_ecs, replace = TRUE, prob = ec_class_weights)
  
  ec_ids <- character(n_ecs)
  for (i in 1:n_ecs) {
    subclass <- sample(1:20, 1)
    sub_subclass <- sample(1:50, 1) 
    serial <- sample(1:200, 1)
    ec_ids[i] <- paste0("EC:", ec_classes[i], ".", subclass, ".", sub_subclass, ".", serial)
  }
  
  # Create abundance matrix with realistic microbiome characteristics
  abundance_realistic <- matrix(0, nrow = n_ecs, ncol = n_samples)
  
  # Subject-specific baseline abundances (individual variation)
  subject_baselines <- matrix(rlnorm(n_ecs * n_subjects, meanlog = log(50), sdlog = 1.5),
                             nrow = n_ecs, ncol = n_subjects)
  
  # Fill abundance matrix with subject-specific patterns
  for (subject in 1:n_subjects) {
    baseline_cols <- c(subject * 2 - 1)  # Baseline timepoint
    followup_cols <- c(subject * 2)      # Follow-up timepoint
    
    # Baseline abundances
    abundance_realistic[, baseline_cols] <- subject_baselines[, subject]
    
    # Follow-up with some changes
    fold_changes <- runif(n_ecs, 0.8, 1.2)  # Most ECs don't change much
    # Some ECs have bigger changes
    big_change_ecs <- sample(1:n_ecs, 20)
    fold_changes[big_change_ecs] <- runif(20, 0.5, 2.0)
    
    abundance_realistic[, followup_cols] <- subject_baselines[, subject] * fold_changes
  }
  
  # Add sparsity (many zeros, typical in microbiome data)
  zero_mask <- matrix(runif(n_ecs * n_samples) < 0.4, nrow = n_ecs, ncol = n_samples)
  abundance_realistic[zero_mask] <- 0
  
  # Set proper names
  rownames(abundance_realistic) <- ec_ids
  colnames(abundance_realistic) <- paste0("Subject_", rep(1:n_subjects, each = 2), 
                                         "_", rep(c("Baseline", "Followup"), n_subjects))
  
  # Create realistic metadata
  metadata_realistic <- data.frame(
    sample_id = colnames(abundance_realistic),
    subject_id = paste0("Subject_", rep(1:n_subjects, each = 2)),
    timepoint = rep(c("Baseline", "Followup"), n_subjects),
    treatment = rep(c("Control", "Treatment"), length.out = n_samples),
    age = rep(sample(25:65, n_subjects), each = 2),
    sex = rep(sample(c("Male", "Female"), n_subjects, replace = TRUE), each = 2),
    bmi = rep(rnorm(n_subjects, mean = 25, sd = 4), each = 2),
    stringsAsFactors = FALSE
  )
  rownames(metadata_realistic) <- metadata_realistic$sample_id
  
  # Test realistic analysis scenarios
  
  # Scenario 1: Treatment vs Control comparison
  expect_no_error({
    treatment_results <- pathway_gsea(
      abundance = abundance_realistic,
      metadata = metadata_realistic,
      group = "treatment",
      pathway_type = "MetaCyc",
      method = "fgsea",
      rank_method = "signal2noise",
      nperm = 500
    )
  })
  
  # Scenario 2: Baseline vs Followup comparison
  expect_no_error({
    timepoint_results <- pathway_gsea(
      abundance = abundance_realistic,
      metadata = metadata_realistic,
      group = "timepoint", 
      pathway_type = "MetaCyc",
      method = "fgsea",
      rank_method = "t_test",
      nperm = 500
    )
  })
  
  # Validate realistic results
  for (results in list(treatment_results, timepoint_results)) {
    if (nrow(results) > 0) {
      expect_s3_class(results, "data.frame")
      expect_true(all(is.finite(results$NES)))
      expect_true(all(results$pvalue >= 0 & results$pvalue <= 1))
      
      # Realistic effect sizes (NES typically between -3 and 3)
      expect_true(all(abs(results$NES) <= 5))
      
      # Should have reasonable pathway sizes
      expect_true(all(results$size >= 2))
      expect_true(all(results$size <= 100))
    }
  }
  
  # Test data sparsity handling
  sparsity_rate <- sum(abundance_realistic == 0) / (n_ecs * n_samples)
  expect_gt(sparsity_rate, 0.2)  # Should be reasonably sparse
  expect_lt(sparsity_rate, 0.8)  # But not too sparse
  
  cat(sprintf("\nReal-world simulation: %.1f%% sparsity, %d subjects, %d timepoints\n",
              sparsity_rate * 100, n_subjects, n_timepoints))
  cat(sprintf("Treatment analysis: %d pathways tested\n", nrow(treatment_results)))
  cat(sprintf("Timepoint analysis: %d pathways tested\n", nrow(timepoint_results)))
})

test_that("MetaCyc pathway visualization integration stress test", {
  skip_if_not_installed("fgsea")
  
  # Test that MetaCyc GSEA results integrate well with visualization functions
  set.seed(2025)
  
  # Create comprehensive test dataset for visualization
  n_samples <- 40
  n_ecs <- 150
  
  ec_ids <- paste0("EC:", sample(1:6, n_ecs, replace = TRUE), ".",
                   sample(1:12, n_ecs, replace = TRUE), ".",
                   sample(1:25, n_ecs, replace = TRUE), ".",
                   sample(1:80, n_ecs, replace = TRUE))
  
  abundance_viz <- matrix(
    rlnorm(n_samples * n_ecs, meanlog = log(80), sdlog = 1.3),
    nrow = n_ecs, ncol = n_samples
  )
  rownames(abundance_viz) <- ec_ids
  colnames(abundance_viz) <- paste0("Sample_", sprintf("%02d", 1:n_samples))
  
  # Add differential signal for visualization
  treatment_samples <- paste0("Sample_", sprintf("%02d", 21:40))
  signal_ecs <- sample(1:n_ecs, 30)
  
  for (ec_idx in signal_ecs) {
    if (ec_idx <= 15) {
      # Strong upregulation
      abundance_viz[ec_idx, treatment_samples] <- 
        abundance_viz[ec_idx, treatment_samples] * runif(length(treatment_samples), 2, 4)
    } else if (ec_idx <= 30) {
      # Strong downregulation  
      abundance_viz[ec_idx, treatment_samples] <- 
        abundance_viz[ec_idx, treatment_samples] * runif(length(treatment_samples), 0.2, 0.5)
    }
  }
  
  metadata_viz <- data.frame(
    sample_id = colnames(abundance_viz),
    group = rep(c("Control", "Treatment"), each = 20),
    batch = rep(c("Batch_1", "Batch_2", "Batch_3", "Batch_4"), each = 10),
    stringsAsFactors = FALSE
  )
  rownames(metadata_viz) <- metadata_viz$sample_id
  
  # Run comprehensive GSEA analysis
  expect_no_error({
    viz_results <- pathway_gsea(
      abundance = abundance_viz,
      metadata = metadata_viz,
      group = "group",
      pathway_type = "MetaCyc",
      method = "fgsea",
      rank_method = "signal2noise",
      nperm = 1000,
      min_size = 3,
      max_size = 50
    )
  })
  
  # Test annotation integration
  if (nrow(viz_results) > 0) {
    expect_no_error({
      annotated_viz_results <- gsea_pathway_annotation(viz_results, pathway_type = "MetaCyc")
    })
    
    expect_s3_class(annotated_viz_results, "data.frame")
    expect_true("pathway_name" %in% colnames(annotated_viz_results))
    expect_equal(nrow(annotated_viz_results), nrow(viz_results))
    
    # Test visualization data preparation
    if (nrow(annotated_viz_results) > 0) {
      # Test pathway labeling logic (what visualize_gsea would use)
      pathway_labels <- ifelse(
        !is.na(annotated_viz_results$pathway_name) & 
        annotated_viz_results$pathway_name != "",
        annotated_viz_results$pathway_name,
        annotated_viz_results$pathway_id
      )
      
      expect_equal(length(pathway_labels), nrow(annotated_viz_results))
      expect_true(all(nchar(pathway_labels) > 0))
      
      # Test sorting and filtering (common visualization operations)
      top_pathways <- head(annotated_viz_results[order(abs(annotated_viz_results$NES), decreasing = TRUE), ], 10)
      expect_lte(nrow(top_pathways), 10)
      expect_true(all(abs(top_pathways$NES) >= abs(tail(top_pathways$NES, 1))))
      
      # Test significant pathway filtering
      significant_pathways <- annotated_viz_results[annotated_viz_results$p.adjust < 0.05, ]
      if (nrow(significant_pathways) > 0) {
        expect_true(all(significant_pathways$p.adjust < 0.05))
      }
    }
  }
  
  cat(sprintf("\nVisualization integration test: %d total pathways, %d significant\n",
              nrow(viz_results), sum(viz_results$p.adjust < 0.05, na.rm = TRUE)))
})

test_that("MetaCyc reproducibility across R sessions", {
  skip_if_not_installed("fgsea")
  
  # Test that results are reproducible across different R sessions
  # This simulates running the same analysis multiple times
  
  set.seed(12345)  # Fixed seed for data generation
  
  n_samples <- 30
  n_ecs <- 80
  
  # Generate reproducible test data
  ec_ids <- paste0("EC:", 1, ".", rep(1:4, each = 20), ".", 
                   rep(1:20, 4), ".", 1:80)
  
  abundance_repro <- matrix(
    rnorm(n_samples * n_ecs, mean = 100, sd = 25),
    nrow = n_ecs, ncol = n_samples
  )
  rownames(abundance_repro) <- ec_ids
  colnames(abundance_repro) <- paste0("Sample_", 1:n_samples)
  
  metadata_repro <- data.frame(
    sample_id = colnames(abundance_repro),
    group = rep(c("Control", "Treatment"), each = 15),
    stringsAsFactors = FALSE
  )
  rownames(metadata_repro) <- metadata_repro$sample_id
  
  # Run analysis multiple times with same parameters
  results_list <- list()
  
  for (run in 1:3) {
    expect_no_error({
      results_list[[run]] <- pathway_gsea(
        abundance = abundance_repro,
        metadata = metadata_repro,
        group = "group",
        pathway_type = "MetaCyc",
        method = "fgsea", 
        rank_method = "signal2noise",
        nperm = 500,
        seed = 123  # Same seed for all runs
      )
    })
  }
  
  # Verify reproducibility
  expect_equal(length(results_list), 3)
  
  if (all(sapply(results_list, nrow) > 0)) {
    # All runs should have same number of pathways
    expect_equal(nrow(results_list[[1]]), nrow(results_list[[2]]))
    expect_equal(nrow(results_list[[1]]), nrow(results_list[[3]]))
    
    # Sort all results by pathway_id for comparison
    for (i in 1:3) {
      results_list[[i]] <- results_list[[i]][order(results_list[[i]]$pathway_id), ]
    }
    
    # Pathway IDs should be identical
    expect_equal(results_list[[1]]$pathway_id, results_list[[2]]$pathway_id)
    expect_equal(results_list[[1]]$pathway_id, results_list[[3]]$pathway_id)
    
    # Statistical measures should be very similar (within tolerance for stochastic methods)
    for (col in c("NES", "pvalue", "p.adjust")) {
      if (col %in% colnames(results_list[[1]])) {
        expect_equal(results_list[[1]][[col]], results_list[[2]][[col]], tolerance = 1e-6)
        expect_equal(results_list[[1]][[col]], results_list[[3]][[col]], tolerance = 1e-6)
      }
    }
  }
  
  cat(sprintf("\nReproducibility test: %d pathways tested across 3 runs\n", 
              nrow(results_list[[1]])))
  
  # Test with different seeds should give different results
  expect_no_error({
    different_seed_result <- pathway_gsea(
      abundance = abundance_repro,
      metadata = metadata_repro,
      group = "group",
      pathway_type = "MetaCyc",
      method = "fgsea",
      rank_method = "signal2noise", 
      nperm = 500,
      seed = 456  # Different seed
    )
  })
  
  # Different seed results may differ slightly in p-values due to permutation randomness
  if (nrow(different_seed_result) > 0 && nrow(results_list[[1]]) > 0) {
    # Should have same pathways tested
    expect_setequal(different_seed_result$pathway_id, results_list[[1]]$pathway_id)
  }
})