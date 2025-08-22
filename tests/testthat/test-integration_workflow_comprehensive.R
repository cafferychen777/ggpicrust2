# Comprehensive integration workflow tests for end-to-end GSEA + DAA analysis
library(testthat)

# Comprehensive workflow test helpers
create_realistic_microbiome_data <- function(n_samples = 24, n_ko_features = 100, 
                                           effect_proportion = 0.25, effect_magnitude = 1.5) {
  set.seed(54321)
  
  # Generate realistic KO identifiers
  ko_ids <- paste0("K", sprintf("%05d", sample(1:25000, n_ko_features)))
  sample_ids <- paste0("Sample_", 1:n_samples)
  
  # Create experimental design with multiple groups
  groups <- rep(c("Control", "Treatment"), each = n_samples / 2)
  batch <- rep(c("Batch1", "Batch2"), times = n_samples / 2)
  
  # Generate abundance data with realistic microbiome characteristics
  abundance_matrix <- matrix(0, nrow = n_ko_features, ncol = n_samples)
  
  for (i in 1:n_ko_features) {
    # Base abundance following negative binomial distribution (common in microbiome)
    base_mean <- rgamma(1, shape = 2, rate = 0.1)  # Varies per feature
    base_dispersion <- 0.5
    
    # Generate baseline abundances
    baseline <- rnbinom(n_samples, mu = base_mean, size = 1/base_dispersion)
    
    # Add group effects for some features
    if (runif(1) < effect_proportion) {
      group_effect <- rnorm(1, 0, effect_magnitude)
      for (j in 1:n_samples) {
        if (groups[j] == "Treatment") {
          new_mean <- base_mean * exp(group_effect)
          baseline[j] <- rnbinom(1, mu = new_mean, size = 1/base_dispersion)
        }
      }
    }
    
    # Add small batch effects
    batch_effect <- rnorm(1, 0, 0.3)
    for (j in 1:n_samples) {
      if (batch[j] == "Batch2") {
        batch_multiplier <- exp(batch_effect * 0.2)
        baseline[j] <- round(baseline[j] * batch_multiplier)
      }
    }
    
    abundance_matrix[i, ] <- pmax(0, baseline)  # Ensure non-negative
  }
  
  rownames(abundance_matrix) <- ko_ids
  colnames(abundance_matrix) <- sample_ids
  
  # Create abundance data frame in ggpicrust2 format
  abundance_df <- data.frame(
    `#NAME` = ko_ids,
    abundance_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Create comprehensive metadata
  metadata <- data.frame(
    sample_id = sample_ids,
    Group = factor(groups),
    Batch = factor(batch),
    Age = sample(25:75, n_samples, replace = TRUE),
    BMI = rnorm(n_samples, 25, 4),
    Sex = factor(sample(c("Male", "Female"), n_samples, replace = TRUE)),
    Collection_Date = seq(as.Date("2023-01-01"), by = "week", length.out = n_samples),
    stringsAsFactors = FALSE
  )
  
  return(list(
    abundance = abundance_df,
    metadata = metadata,
    true_effects = abundance_matrix,
    design_groups = groups
  ))
}

create_mock_pipeline_results <- function(ko_features, n_pathways = 30, overlap_rate = 0.6) {
  # Simulate realistic pipeline results with controlled overlap
  
  # Create DAA results with some significant features
  n_daa_sig <- round(length(ko_features) * 0.3)
  daa_significant <- sample(ko_features, n_daa_sig)
  
  daa_results <- data.frame(
    feature = ko_features,
    method = rep("LinDA", length(ko_features)),
    group1 = rep("Control", length(ko_features)),
    group2 = rep("Treatment", length(ko_features)),
    p_values = ifelse(ko_features %in% daa_significant,
                     runif(length(ko_features), 0.001, 0.045),
                     runif(length(ko_features), 0.06, 0.3)),
    p_adjust = ifelse(ko_features %in% daa_significant,
                     runif(length(ko_features), 0.001, 0.049),
                     runif(length(ko_features), 0.055, 0.4)),
    log_2_fold_change = ifelse(ko_features %in% daa_significant,
                              rnorm(length(ko_features), 0, 2),
                              rnorm(length(ko_features), 0, 0.8)),
    effect_size = rnorm(length(ko_features), 0, 1),
    stringsAsFactors = FALSE
  )
  
  # Create GSEA pathway results
  pathway_ids <- paste0("ko", sprintf("%05d", sample(100:2000, n_pathways)))
  n_gsea_sig <- round(n_pathways * 0.35)
  
  # Control overlap between DAA and GSEA significant results
  n_overlap_sig <- round(min(n_daa_sig, n_gsea_sig) * overlap_rate)
  gsea_significant <- c(
    sample(daa_significant, n_overlap_sig),  # Overlapping significant
    sample(setdiff(pathway_ids, daa_significant), n_gsea_sig - n_overlap_sig)  # GSEA-only significant
  )
  
  gsea_results <- data.frame(
    pathway_id = pathway_ids,
    pathway_name = paste("KEGG Pathway", pathway_ids),
    pathway_class = sample(c("Metabolism", "Genetic Information Processing", 
                           "Environmental Information Processing", "Cellular Processes"),
                          n_pathways, replace = TRUE),
    size = sample(15:200, n_pathways, replace = TRUE),
    ES = ifelse(pathway_ids %in% gsea_significant,
               rnorm(n_pathways, 0, 1.5),
               rnorm(n_pathways, 0, 0.6)),
    NES = ifelse(pathway_ids %in% gsea_significant,
                rnorm(n_pathways, 0, 2),
                rnorm(n_pathways, 0, 0.8)),
    pvalue = ifelse(pathway_ids %in% gsea_significant,
                   runif(n_pathways, 0.001, 0.045),
                   runif(n_pathways, 0.06, 0.25)),
    p.adjust = ifelse(pathway_ids %in% gsea_significant,
                     runif(n_pathways, 0.001, 0.049),
                     runif(n_pathways, 0.055, 0.3)),
    leading_edge = replicate(n_pathways, {
      genes <- sample(ko_features, sample(5:30, 1))
      paste(genes, collapse = ";")
    }),
    method = rep("fgsea", n_pathways),
    stringsAsFactors = FALSE
  )
  
  return(list(
    daa_results = list(list(
      plot = ggplot2::ggplot() + ggplot2::geom_blank() + ggplot2::ggtitle("DAA Results"),
      results = daa_results
    )),
    gsea_results = gsea_results,
    expected_overlap = intersect(daa_significant, gsea_significant),
    daa_significant = daa_significant,
    gsea_significant = gsea_significant
  ))
}

test_that("End-to-end GSEA + DAA integration workflow executes completely", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  # Create realistic test data
  microbiome_data <- create_realistic_microbiome_data(n_samples = 20, n_ko_features = 80)
  mock_results <- create_mock_pipeline_results(microbiome_data$abundance$`#NAME`, n_pathways = 25)
  
  # Mock all pipeline components with realistic behavior
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    return(mock_results$daa_results)
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    # Verify parameters passed through correctly
    args <- list(...)
    expect_true("abundance" %in% names(args))
    expect_true("metadata" %in% names(args))
    expect_true("group" %in% names(args))
    
    return(mock_results$gsea_results)
  })
  
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(gsea_results, ...) {
    # Add pathway class annotations
    if (!"pathway_class" %in% colnames(gsea_results)) {
      gsea_results$pathway_class <- sample(c("Metabolism", "Genetic Information Processing"),
                                          nrow(gsea_results), replace = TRUE)
    }
    return(gsea_results)
  })
  
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(...) {
    ggplot2::ggplot() + 
      ggplot2::geom_col(ggplot2::aes(x = 1:10, y = rnorm(10))) +
      ggplot2::ggtitle("GSEA Visualization")
  })
  
  mockery::stub(ggpicrust2_extended, "compare_gsea_daa", function(gsea_results, daa_results, ...) {
    # Calculate actual comparison based on mock data
    gsea_sig <- gsea_results$pathway_id[gsea_results$p.adjust < 0.05]
    daa_sig <- daa_results$feature[daa_results$p_adjust < 0.05]
    
    overlap <- intersect(gsea_sig, daa_sig)
    gsea_only <- setdiff(gsea_sig, daa_sig)
    daa_only <- setdiff(daa_sig, gsea_sig)
    
    list(
      plot = ggplot2::ggplot() + ggplot2::geom_blank() + ggplot2::ggtitle("DAA vs GSEA Comparison"),
      results = list(
        overlap = overlap,
        gsea_only = gsea_only,
        daa_only = daa_only,
        n_overlap = length(overlap),
        n_gsea_only = length(gsea_only),
        n_daa_only = length(daa_only),
        n_gsea_total = length(gsea_sig),
        n_daa_total = length(daa_sig)
      )
    )
  })
  
  # Execute complete workflow
  workflow_results <- ggpicrust2_extended(
    data = microbiome_data$abundance,
    metadata = microbiome_data$metadata,
    group = "Group",
    pathway = "KO",
    daa_method = "LinDA",
    ko_to_kegg = TRUE,
    run_gsea = TRUE,
    gsea_params = list(
      method = "fgsea",
      rank_method = "signal2noise",
      nperm = 1000,
      min_size = 10,
      max_size = 500
    )
  )
  
  # Verify complete workflow results structure
  expect_type(workflow_results, "list")
  expect_setequal(names(workflow_results), c("daa_results", "gsea_results", "gsea_plot", "comparison"))
  
  # Verify DAA results
  expect_type(workflow_results$daa_results, "list")
  expect_s3_class(workflow_results$daa_results[[1]]$plot, "ggplot")
  expect_true(is.data.frame(workflow_results$daa_results[[1]]$results))
  
  # Verify GSEA results
  expect_true(is.data.frame(workflow_results$gsea_results))
  expect_true(all(c("pathway_id", "NES", "p.adjust") %in% colnames(workflow_results$gsea_results)))
  
  # Verify GSEA visualization
  expect_s3_class(workflow_results$gsea_plot, "ggplot")
  
  # Verify comparison results
  expect_type(workflow_results$comparison, "list")
  expect_true("plot" %in% names(workflow_results$comparison))
  expect_true("results" %in% names(workflow_results$comparison))
  expect_s3_class(workflow_results$comparison$plot, "ggplot")
})

test_that("Integration workflow data flow between functions is consistent", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  # Create controlled test data to track data flow
  test_data <- create_realistic_microbiome_data(n_samples = 16, n_ko_features = 50)
  
  # Track data passed between functions
  captured_data <- list()
  
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    captured_data$ggpicrust2_args <<- list(...)
    create_mock_pipeline_results(test_data$abundance$`#NAME`)$daa_results
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    captured_data$pathway_gsea_args <<- list(...)
    create_mock_pipeline_results(test_data$abundance$`#NAME`)$gsea_results
  })
  
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(gsea_results, pathway_type) {
    captured_data$annotation_args <<- list(gsea_results = gsea_results, pathway_type = pathway_type)
    gsea_results$pathway_class <- "Metabolism"
    return(gsea_results)
  })
  
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(gsea_results, ...) {
    captured_data$visualize_args <<- list(gsea_results = gsea_results, ...)
    ggplot2::ggplot() + ggplot2::geom_blank()
  })
  
  mockery::stub(ggpicrust2_extended, "compare_gsea_daa", function(gsea_results, daa_results, ...) {
    captured_data$comparison_args <<- list(gsea_results = gsea_results, daa_results = daa_results, ...)
    list(plot = ggplot2::ggplot(), results = list(n_overlap = 5))
  })
  
  # Execute workflow with specific parameters
  workflow_results <- ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Group",
    pathway = "KO",
    daa_method = "ALDEx2",
    run_gsea = TRUE,
    gsea_params = list(
      method = "fgsea",
      rank_method = "log2FC",
      nperm = 500
    )
  )
  
  # Verify data flow consistency
  
  # Check ggpicrust2 received original parameters
  expect_equal(captured_data$ggpicrust2_args$group, "Group")
  expect_equal(captured_data$ggpicrust2_args$pathway, "KO")
  expect_equal(captured_data$ggpicrust2_args$daa_method, "ALDEx2")
  
  # Check pathway_gsea received processed abundance data and merged parameters
  expect_true("abundance" %in% names(captured_data$pathway_gsea_args))
  expect_equal(captured_data$pathway_gsea_args$group, "Group")
  expect_equal(captured_data$pathway_gsea_args$method, "fgsea")
  expect_equal(captured_data$pathway_gsea_args$rank_method, "log2FC")
  expect_equal(captured_data$pathway_gsea_args$nperm, 500)
  
  # Check annotation received GSEA results
  expect_true(is.data.frame(captured_data$annotation_args$gsea_results))
  expect_equal(captured_data$annotation_args$pathway_type, "KEGG")
  
  # Check visualization received annotated GSEA results
  expect_true(is.data.frame(captured_data$visualize_args$gsea_results))
  
  # Check comparison received both result sets
  expect_true(is.data.frame(captured_data$comparison_args$gsea_results))
  expect_true(is.data.frame(captured_data$comparison_args$daa_results))
})

test_that("Integration workflow result format consistency validation", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  # Create test data
  test_data <- create_realistic_microbiome_data(n_samples = 18, n_ko_features = 60)
  
  # Create mock results with specific formats to validate consistency
  mock_daa_results <- data.frame(
    feature = paste0("K", sprintf("%05d", 1:30)),
    method = rep("LinDA", 30),
    p_adjust = runif(30, 0.001, 0.2),
    log_2_fold_change = rnorm(30, 0, 1.5),
    stringsAsFactors = FALSE
  )
  
  mock_gsea_results <- data.frame(
    pathway_id = paste0("ko", sprintf("%05d", 1:25)),
    pathway_name = paste("Pathway", 1:25),
    NES = rnorm(25, 0, 2),
    p.adjust = runif(25, 0.001, 0.15),
    stringsAsFactors = FALSE
  )
  
  # Mock functions to return controlled formats
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    list(list(
      plot = ggplot2::ggplot() + ggplot2::geom_blank(),
      results = mock_daa_results
    ))
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) mock_gsea_results)
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(gsea_results, ...) {
    gsea_results$pathway_class <- "Metabolism"
    return(gsea_results)
  })
  
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(...) {
    ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(x = 1:5, y = 1:5))
  })
  
  mockery::stub(ggpicrust2_extended, "compare_gsea_daa", function(...) {
    list(
      plot = ggplot2::ggplot() + ggplot2::geom_bar(ggplot2::aes(x = c("A", "B"), y = c(1, 2)), stat = "identity"),
      results = list(
        overlap = c("K00001", "K00002"),
        n_overlap = 2,
        n_gsea_total = 10,
        n_daa_total = 12
      )
    )
  })
  
  # Execute workflow
  results <- ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Group",
    pathway = "KO",
    run_gsea = TRUE
  )
  
  # Validate result structure consistency
  expect_type(results, "list")
  expect_equal(length(results), 4)
  expect_named(results, c("daa_results", "gsea_results", "gsea_plot", "comparison"))
  
  # Validate DAA results format
  expect_type(results$daa_results, "list")
  expect_s3_class(results$daa_results[[1]]$plot, "ggplot")
  expect_true(is.data.frame(results$daa_results[[1]]$results))
  expect_true(all(c("feature", "p_adjust", "log_2_fold_change") %in% colnames(results$daa_results[[1]]$results)))
  
  # Validate GSEA results format
  expect_true(is.data.frame(results$gsea_results))
  expect_true(all(c("pathway_id", "pathway_name", "NES", "p.adjust") %in% colnames(results$gsea_results)))
  expect_true("pathway_class" %in% colnames(results$gsea_results))  # Should be added by annotation
  
  # Validate visualization format
  expect_s3_class(results$gsea_plot, "ggplot")
  
  # Validate comparison results format
  expect_type(results$comparison, "list")
  expect_s3_class(results$comparison$plot, "ggplot")
  expect_type(results$comparison$results, "list")
  expect_true(all(c("overlap", "n_overlap", "n_gsea_total", "n_daa_total") %in% names(results$comparison$results)))
})

test_that("Integration workflow error propagation handling works correctly", {
  skip_if_not_installed("dplyr")
  
  test_data <- create_realistic_microbiome_data(n_samples = 12, n_ko_features = 30)
  
  # Test error in DAA analysis propagates immediately
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    stop("DAA analysis failed: insufficient sample size")
  })
  
  expect_error(
    ggpicrust2_extended(
      data = test_data$abundance,
      metadata = test_data$metadata,
      group = "Group",
      pathway = "KO",
      run_gsea = FALSE
    ),
    "DAA analysis failed: insufficient sample size"
  )
  
  # Test error in GSEA analysis when run_gsea = TRUE
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    list(list(plot = ggplot2::ggplot(), results = data.frame(feature = character(0))))
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    stop("GSEA failed: no valid gene sets found")
  })
  
  expect_error(
    ggpicrust2_extended(
      data = test_data$abundance,
      metadata = test_data$metadata,
      group = "Group",
      pathway = "KO",
      run_gsea = TRUE
    ),
    "GSEA failed: no valid gene sets found"
  )
  
  # Test error in downstream analysis (annotation/visualization)
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    list(list(plot = ggplot2::ggplot(), results = data.frame(feature = paste0("K", 1:10))))
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    data.frame(pathway_id = paste0("ko", 1:5), p.adjust = runif(5, 0.01, 0.1))
  })
  
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(...) {
    stop("Annotation database unavailable")
  })
  
  expect_error(
    ggpicrust2_extended(
      data = test_data$abundance,
      metadata = test_data$metadata,
      group = "Group",
      pathway = "KO",
      run_gsea = TRUE
    ),
    "Annotation database unavailable"
  )
})

test_that("Integration workflow handles different experimental designs correctly", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  # Test with two-group design
  two_group_data <- create_realistic_microbiome_data(n_samples = 20, n_ko_features = 40)
  
  # Test with multi-group design
  multi_group_metadata <- data.frame(
    sample_id = paste0("Sample_", 1:24),
    Treatment = factor(rep(c("Control", "Low_Dose", "High_Dose"), each = 8)),
    stringsAsFactors = FALSE
  )
  
  mock_pipeline <- create_mock_pipeline_results(two_group_data$abundance$`#NAME`)
  
  # Mock functions for multi-group test
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) mock_pipeline$daa_results)
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    args <- list(...)
    expect_equal(args$group, "Treatment")  # Verify group parameter passed correctly
    mock_pipeline$gsea_results
  })
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(gsea_results, ...) gsea_results)
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(...) ggplot2::ggplot())
  mockery::stub(ggpicrust2_extended, "compare_gsea_daa", function(...) {
    list(plot = ggplot2::ggplot(), results = list(n_overlap = 3))
  })
  
  # Test two-group design
  results_two_group <- ggpicrust2_extended(
    data = two_group_data$abundance,
    metadata = two_group_data$metadata,
    group = "Group",
    pathway = "KO",
    run_gsea = TRUE
  )
  
  expect_type(results_two_group, "list")
  expect_equal(length(results_two_group), 4)
  
  # Test multi-group design
  results_multi_group <- ggpicrust2_extended(
    data = two_group_data$abundance,  # Reuse abundance data
    metadata = multi_group_metadata,
    group = "Treatment",
    pathway = "KO",
    run_gsea = TRUE
  )
  
  expect_type(results_multi_group, "list")
  expect_equal(length(results_multi_group), 4)
})

test_that("Integration workflow performance with larger datasets", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  # Create larger test dataset
  large_data <- create_realistic_microbiome_data(n_samples = 50, n_ko_features = 200)
  large_mock <- create_mock_pipeline_results(large_data$abundance$`#NAME`, n_pathways = 100)
  
  # Mock all functions with performance tracking
  function_call_times <- list()
  
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    start_time <- Sys.time()
    Sys.sleep(0.01)  # Simulate some computation time
    function_call_times$ggpicrust2 <<- Sys.time() - start_time
    large_mock$daa_results
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    start_time <- Sys.time()
    Sys.sleep(0.02)  # Simulate GSEA computation time
    function_call_times$pathway_gsea <<- Sys.time() - start_time
    large_mock$gsea_results
  })
  
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(gsea_results, ...) {
    start_time <- Sys.time()
    gsea_results$pathway_class <- "Metabolism"
    function_call_times$annotation <<- Sys.time() - start_time
    gsea_results
  })
  
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(...) {
    start_time <- Sys.time()
    plot <- ggplot2::ggplot() + ggplot2::geom_blank()
    function_call_times$visualization <<- Sys.time() - start_time
    plot
  })
  
  mockery::stub(ggpicrust2_extended, "compare_gsea_daa", function(...) {
    start_time <- Sys.time()
    result <- list(plot = ggplot2::ggplot(), results = list(n_overlap = 15))
    function_call_times$comparison <<- Sys.time() - start_time
    result
  })
  
  # Execute workflow and measure total time
  total_start <- Sys.time()
  results <- ggpicrust2_extended(
    data = large_data$abundance,
    metadata = large_data$metadata,
    group = "Group",
    pathway = "KO",
    run_gsea = TRUE,
    gsea_params = list(nperm = 100)  # Reduce permutations for speed
  )
  total_time <- Sys.time() - total_start
  
  # Verify results structure
  expect_type(results, "list")
  expect_equal(length(results), 4)
  
  # Verify reasonable performance (should complete in under 10 seconds with mocking)
  expect_true(as.numeric(total_time, units = "secs") < 10)
  
  # Verify all functions were called
  expect_true("ggpicrust2" %in% names(function_call_times))
  expect_true("pathway_gsea" %in% names(function_call_times))
  expect_true("annotation" %in% names(function_call_times))
  expect_true("visualization" %in% names(function_call_times))
  expect_true("comparison" %in% names(function_call_times))
})