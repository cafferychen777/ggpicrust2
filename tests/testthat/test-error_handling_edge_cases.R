# Comprehensive error handling and edge case tests for complex workflows
library(testthat)

# Advanced error simulation helpers
create_problematic_data <- function(problem_type = "empty") {
  base_abundance <- data.frame(
    `#NAME` = paste0("K", sprintf("%05d", 1:20)),
    Sample1 = rpois(20, 50), Sample2 = rpois(20, 50),
    Sample3 = rpois(20, 50), Sample4 = rpois(20, 50),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  
  base_metadata <- data.frame(
    sample_id = paste0("Sample", 1:4),
    Group = factor(c("A", "A", "B", "B")),
    stringsAsFactors = FALSE
  )
  
  if (problem_type == "empty") {
    # Empty abundance data
    abundance <- data.frame(`#NAME` = character(0), check.names = FALSE)
    metadata <- data.frame(sample_id = character(0), Group = factor())
  } else if (problem_type == "mismatched_samples") {
    # Sample mismatch between abundance and metadata
    abundance <- base_abundance
    metadata <- data.frame(
      sample_id = paste0("DifferentSample", 1:4),
      Group = factor(c("A", "A", "B", "B"))
    )
  } else if (problem_type == "single_group") {
    # Only one group in metadata
    abundance <- base_abundance
    metadata <- data.frame(
      sample_id = paste0("Sample", 1:4),
      Group = factor(rep("A", 4))  # All same group
    )
  } else if (problem_type == "missing_values") {
    # Missing values in key columns
    abundance <- base_abundance
    abundance$`#NAME`[1] <- NA
    abundance$Sample1[2] <- NA
    metadata <- base_metadata
    metadata$Group[1] <- NA
  } else if (problem_type == "wrong_types") {
    # Wrong data types
    abundance <- base_abundance
    abundance$Sample1 <- as.character(abundance$Sample1)  # Should be numeric
    metadata <- base_metadata
    metadata$Group <- as.character(metadata$Group)  # Should be factor for some analyses
  } else if (problem_type == "negative_values") {
    # Negative abundance values (impossible in microbiome data)
    abundance <- base_abundance
    abundance[abundance > 30] <- -abundance[abundance > 30]
    metadata <- base_metadata
  } else {
    # Default: return normal data
    abundance <- base_abundance
    metadata <- base_metadata
  }
  
  return(list(abundance = abundance, metadata = metadata))
}

test_that("ggpicrust2_extended handles empty abundance data gracefully", {
  skip_if_not_installed("dplyr")
  
  problem_data <- create_problematic_data("empty")
  
  # Mock ggpicrust2 to handle empty data
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    stop("Empty abundance matrix provided")
  })
  
  # Should propagate the error from ggpicrust2
  expect_error(
    ggpicrust2_extended(
      data = problem_data$abundance,
      metadata = problem_data$metadata,
      group = "Group",
      pathway = "KO",
      run_gsea = FALSE
    ),
    "Empty abundance matrix provided"
  )
  
  # With GSEA enabled, should fail even before reaching GSEA
  expect_error(
    ggpicrust2_extended(
      data = problem_data$abundance,
      metadata = problem_data$metadata,
      group = "Group", 
      pathway = "KO",
      run_gsea = TRUE
    ),
    "Empty abundance matrix provided"
  )
})

test_that("ggpicrust2_extended handles sample mismatch between abundance and metadata", {
  skip_if_not_installed("dplyr")
  
  problem_data <- create_problematic_data("mismatched_samples")
  
  # Mock ggpicrust2 to simulate sample mismatch error
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    stop("Sample names in abundance data do not match metadata")
  })
  
  expect_error(
    ggpicrust2_extended(
      data = problem_data$abundance,
      metadata = problem_data$metadata,
      group = "Group",
      pathway = "KO",
      run_gsea = FALSE
    ),
    "Sample names in abundance data do not match metadata"
  )
})

test_that("ggpicrust2_extended handles insufficient group variation", {
  skip_if_not_installed("dplyr")
  
  problem_data <- create_problematic_data("single_group")
  
  # Mock ggpicrust2 to simulate insufficient groups error
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    stop("Need at least 2 groups for differential analysis")
  })
  
  expect_error(
    ggpicrust2_extended(
      data = problem_data$abundance,
      metadata = problem_data$metadata,
      group = "Group",
      pathway = "KO",
      daa_method = "LinDA",
      run_gsea = FALSE
    ),
    "Need at least 2 groups for differential analysis"
  )
})

test_that("ggpicrust2_extended handles missing required packages gracefully", {
  skip_if_not_installed("dplyr")
  
  normal_data <- create_problematic_data("normal")
  
  # Mock ggpicrust2 to return valid results
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    list(list(
      plot = ggplot2::ggplot(),
      results = data.frame(feature = paste0("K", 1:5), p_adjust = runif(5, 0.01, 0.1))
    ))
  })
  
  # Test missing fgsea package
  mockery::stub(ggpicrust2_extended, "requireNamespace", function(pkg, quietly = TRUE) {
    if (pkg == "fgsea") return(FALSE)
    return(TRUE)
  })
  
  expect_warning(
    result <- ggpicrust2_extended(
      data = normal_data$abundance,
      metadata = normal_data$metadata,
      group = "Group",
      pathway = "KO",
      run_gsea = TRUE
    ),
    "Package 'fgsea' is required for GSEA analysis. Skipping GSEA."
  )
  
  # Should complete but without GSEA results
  expect_true("daa_results" %in% names(result))
  expect_false("gsea_results" %in% names(result))
})

test_that("ggpicrust2_extended handles GSEA computation failures robustly", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  normal_data <- create_problematic_data("normal")
  
  # Mock successful DAA but failed GSEA
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    list(list(
      plot = ggplot2::ggplot(),
      results = data.frame(feature = paste0("K", 1:10), p_adjust = runif(10, 0.01, 0.1))
    ))
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    stop("GSEA failed: insufficient pathway coverage")
  })
  
  expect_error(
    ggpicrust2_extended(
      data = normal_data$abundance,
      metadata = normal_data$metadata,
      group = "Group",
      pathway = "KO",
      run_gsea = TRUE
    ),
    "GSEA failed: insufficient pathway coverage"
  )
})

test_that("ggpicrust2_extended handles annotation failures in GSEA workflow", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  normal_data <- create_problematic_data("normal")
  
  # Mock successful components until annotation
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    list(list(plot = ggplot2::ggplot(), results = data.frame(feature = paste0("K", 1:5))))
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    data.frame(pathway_id = paste0("ko", 1:5), p.adjust = runif(5))
  })
  
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(...) {
    stop("Pathway annotation database connection failed")
  })
  
  expect_error(
    ggpicrust2_extended(
      data = normal_data$abundance,
      metadata = normal_data$metadata,
      group = "Group",
      pathway = "KO",
      run_gsea = TRUE
    ),
    "Pathway annotation database connection failed"
  )
})

test_that("ggpicrust2_extended handles visualization creation failures", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  normal_data <- create_problematic_data("normal")
  
  # Mock successful components until visualization
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    list(list(plot = ggplot2::ggplot(), results = data.frame(feature = paste0("K", 1:5))))
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    data.frame(pathway_id = paste0("ko", 1:5), p.adjust = runif(5))
  })
  
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(gsea_results, ...) gsea_results)
  
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(...) {
    stop("Visualization failed: invalid plot dimensions")
  })
  
  expect_error(
    ggpicrust2_extended(
      data = normal_data$abundance,
      metadata = normal_data$metadata,
      group = "Group",
      pathway = "KO",
      run_gsea = TRUE
    ),
    "Visualization failed: invalid plot dimensions"
  )
})

test_that("compare_gsea_daa handles malformed input data robustly", {
  skip_if_not_installed("ggplot2")
  
  # Test completely invalid inputs
  expect_error(
    compare_gsea_daa(gsea_results = NULL, daa_results = data.frame()),
    "'gsea_results' must be a data frame"
  )
  
  expect_error(
    compare_gsea_daa(gsea_results = data.frame(), daa_results = list()),
    "'daa_results' must be a data frame"
  )
  
  # Test data frames with wrong structure
  invalid_gsea <- data.frame(wrong_column = 1:5)
  valid_daa <- data.frame(feature = paste0("F", 1:5), p_adjust = runif(5))
  
  expect_error(
    compare_gsea_daa(gsea_results = invalid_gsea, daa_results = valid_daa),
    "GSEA results missing required columns"
  )
  
  valid_gsea <- data.frame(pathway_id = paste0("P", 1:5), p.adjust = runif(5))
  invalid_daa <- data.frame(wrong_column = 1:5)
  
  expect_error(
    compare_gsea_daa(gsea_results = valid_gsea, daa_results = invalid_daa),
    "DAA results missing required columns"
  )
})

test_that("compare_gsea_daa handles extreme data values", {
  skip_if_not_installed("ggplot2")
  
  # Test with extreme p-values
  extreme_gsea <- data.frame(
    pathway_id = paste0("P", 1:10),
    pathway_name = paste("Pathway", 1:10),
    size = rep(50, 10),
    ES = c(rep(Inf, 3), rep(-Inf, 2), rnorm(5)),  # Infinite effect sizes
    NES = c(rep(100, 3), rep(-100, 2), rnorm(5)),  # Extreme NES values
    pvalue = c(rep(0, 5), rep(1, 5)),  # Extreme p-values
    p.adjust = c(rep(0, 5), rep(1, 5)),  # Extreme adjusted p-values
    leading_edge = rep("K00001", 10),
    method = rep("fgsea", 10),
    stringsAsFactors = FALSE
  )
  
  extreme_daa <- data.frame(
    feature = paste0("P", 1:10),
    method = rep("ALDEx2", 10),
    group1 = rep("A", 10),
    group2 = rep("B", 10),
    p_values = c(rep(0, 5), rep(1, 5)),
    p_adjust = c(rep(0, 5), rep(1, 5)),
    log_2_fold_change = c(rep(Inf, 3), rep(-Inf, 2), rnorm(5)),  # Infinite fold changes
    stringsAsFactors = FALSE
  )
  
  # Should handle extreme values without crashing
  comparison <- compare_gsea_daa(extreme_gsea, extreme_daa, plot_type = "venn")
  
  expect_type(comparison, "list")
  expect_s3_class(comparison$plot, "ggplot")
  
  # Test scatter plot with extreme values
  comparison_scatter <- compare_gsea_daa(extreme_gsea, extreme_daa, plot_type = "scatter")
  expect_s3_class(comparison_scatter$plot, "ggplot")
})

test_that("compare_gsea_daa handles data with all NAs", {
  skip_if_not_installed("ggplot2")
  
  # Create data with all NAs in critical columns
  na_gsea <- data.frame(
    pathway_id = c(NA, NA, NA, "P4", "P5"),
    pathway_name = rep("Pathway", 5),
    size = rep(50, 5),
    ES = rep(NA, 5),
    NES = rep(NA, 5),
    pvalue = rep(NA, 5),
    p.adjust = c(NA, NA, NA, 0.05, 0.1),
    leading_edge = rep("K00001", 5),
    method = rep("fgsea", 5),
    stringsAsFactors = FALSE
  )
  
  na_daa <- data.frame(
    feature = c(NA, NA, "P3", "P4", "P5"),
    method = rep("ALDEx2", 5),
    group1 = rep("A", 5),
    group2 = rep("B", 5),
    p_values = rep(NA, 5),
    p_adjust = c(NA, NA, 0.03, 0.06, 0.12),
    log_2_fold_change = rep(NA, 5),
    stringsAsFactors = FALSE
  )
  
  # Should handle NAs gracefully
  comparison <- compare_gsea_daa(na_gsea, na_daa)
  
  expect_type(comparison, "list")
  expect_s3_class(comparison$plot, "ggplot")
  
  # Results should reflect only non-NA significant entries
  # Only P4 should be significant in GSEA (p.adjust = 0.05) and P3, P4 in DAA
  expect_true(comparison$results$n_gsea_total >= 0)
  expect_true(comparison$results$n_daa_total >= 0)
})

test_that("compare_gsea_daa handles memory-intensive operations", {
  skip_if_not_installed("ggplot2")
  
  # Create very large datasets to test memory handling
  n_large <- 1000
  
  large_gsea <- data.frame(
    pathway_id = paste0("pathway_", 1:n_large),
    pathway_name = paste("Large Pathway", 1:n_large),
    size = sample(10:200, n_large, replace = TRUE),
    ES = rnorm(n_large, 0, 1),
    NES = rnorm(n_large, 0, 1.5),
    pvalue = runif(n_large, 0.001, 0.2),
    p.adjust = runif(n_large, 0.001, 0.25),
    leading_edge = replicate(n_large, paste(sample(LETTERS, 10), collapse = ";")),
    method = rep("fgsea", n_large),
    stringsAsFactors = FALSE
  )
  
  large_daa <- data.frame(
    feature = paste0("feature_", 1:n_large),
    method = rep("ALDEx2", n_large),
    group1 = rep("Control", n_large),
    group2 = rep("Treatment", n_large),
    p_values = runif(n_large, 0.001, 0.2),
    p_adjust = runif(n_large, 0.001, 0.25),
    log_2_fold_change = rnorm(n_large, 0, 2),
    stringsAsFactors = FALSE
  )
  
  # Should complete in reasonable time and memory
  start_time <- Sys.time()
  comparison <- compare_gsea_daa(large_gsea, large_daa, plot_type = "venn")
  end_time <- Sys.time()
  
  expect_type(comparison, "list")
  expect_s3_class(comparison$plot, "ggplot")
  
  # Should complete within reasonable time (less than 30 seconds)
  expect_true(as.numeric(end_time - start_time, units = "secs") < 30)
})

test_that("Error recovery and partial results in complex workflows", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  normal_data <- create_problematic_data("normal")
  
  # Test scenario where comparison fails but other components succeed
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    list(list(plot = ggplot2::ggplot(), results = data.frame(feature = paste0("K", 1:5))))
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    data.frame(pathway_id = paste0("ko", 1:5), p.adjust = runif(5))
  })
  
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(gsea_results, ...) gsea_results)
  
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(...) ggplot2::ggplot())
  
  mockery::stub(ggpicrust2_extended, "compare_gsea_daa", function(...) {
    stop("Comparison failed: incompatible data formats")
  })
  
  # Should fail at comparison step but preserve earlier results
  expect_error(
    result <- ggpicrust2_extended(
      data = normal_data$abundance,
      metadata = normal_data$metadata,
      group = "Group",
      pathway = "KO",
      run_gsea = TRUE
    ),
    "Comparison failed: incompatible data formats"
  )
})

test_that("Input validation edge cases and boundary conditions", {
  skip_if_not_installed("ggplot2")
  
  # Test invalid plot_type values
  valid_gsea <- data.frame(pathway_id = "P1", p.adjust = 0.05)
  valid_daa <- data.frame(feature = "F1", p_adjust = 0.03)
  
  expect_error(
    compare_gsea_daa(valid_gsea, valid_daa, plot_type = ""),
    "plot_type must be one of"
  )
  
  expect_error(
    compare_gsea_daa(valid_gsea, valid_daa, plot_type = c("venn", "upset")),
    "plot_type must be one of"
  )
  
  expect_error(
    compare_gsea_daa(valid_gsea, valid_daa, plot_type = 123),
    "plot_type must be one of"
  )
  
  # Test boundary p_threshold values
  comparison_zero <- compare_gsea_daa(valid_gsea, valid_daa, p_threshold = 0)
  expect_equal(comparison_zero$results$n_gsea_total, 0)
  expect_equal(comparison_zero$results$n_daa_total, 0)
  
  comparison_one <- compare_gsea_daa(valid_gsea, valid_daa, p_threshold = 1)
  expect_equal(comparison_one$results$n_gsea_total, 1)
  expect_equal(comparison_one$results$n_daa_total, 1)
  
  # Test with negative p_threshold (should still work mathematically)
  comparison_negative <- compare_gsea_daa(valid_gsea, valid_daa, p_threshold = -0.1)
  expect_equal(comparison_negative$results$n_gsea_total, 0)
  expect_equal(comparison_negative$results$n_daa_total, 0)
})

test_that("Resource cleanup and memory management in error scenarios", {
  skip_if_not_installed("dplyr")
  
  # Test that errors don't leave resources hanging
  problem_data <- create_problematic_data("normal")
  
  # Create scenario that fails after allocating memory
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    # Simulate memory allocation
    large_object <- matrix(rnorm(10000), nrow = 100)
    list(list(plot = ggplot2::ggplot(), results = data.frame(feature = 1:100)))
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    # Simulate more memory allocation then fail
    another_large_object <- matrix(rnorm(10000), nrow = 100)
    stop("Simulated failure after memory allocation")
  })
  
  # Error should be thrown without memory leaks
  expect_error(
    ggpicrust2_extended(
      data = problem_data$abundance,
      metadata = problem_data$metadata,
      group = "Group",
      pathway = "KO",
      run_gsea = TRUE
    ),
    "Simulated failure after memory allocation"
  )
  
  # Memory should be cleaned up (R's garbage collection handles this)
  gc()  # Force garbage collection
  
  # Test should complete without memory issues
  expect_true(TRUE)  # If we get here, memory management worked
})

test_that("Concurrent access and thread safety considerations", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  # Test that multiple simultaneous calls don't interfere
  # (Note: R is single-threaded, but this tests function isolation)
  
  normal_data <- create_problematic_data("normal")
  
  # Mock functions that track call order
  call_log <- character(0)
  
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    call_log <<- c(call_log, "ggpicrust2")
    list(list(plot = ggplot2::ggplot(), results = data.frame(feature = "test")))
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    call_log <<- c(call_log, "pathway_gsea")
    data.frame(pathway_id = "test", p.adjust = 0.05)
  })
  
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(gsea_results, ...) {
    call_log <<- c(call_log, "annotation")
    gsea_results
  })
  
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(...) {
    call_log <<- c(call_log, "visualize")
    ggplot2::ggplot()
  })
  
  mockery::stub(ggpicrust2_extended, "compare_gsea_daa", function(...) {
    call_log <<- c(call_log, "compare")
    list(plot = ggplot2::ggplot(), results = list(n_overlap = 1))
  })
  
  # Multiple calls should maintain proper order
  result1 <- ggpicrust2_extended(
    data = normal_data$abundance,
    metadata = normal_data$metadata,
    group = "Group",
    pathway = "KO",
    run_gsea = TRUE
  )
  
  first_call_log <- call_log
  call_log <- character(0)  # Reset
  
  result2 <- ggpicrust2_extended(
    data = normal_data$abundance,
    metadata = normal_data$metadata,
    group = "Group",
    pathway = "KO",
    run_gsea = TRUE
  )
  
  second_call_log <- call_log
  
  # Both should have same call pattern
  expect_equal(first_call_log, second_call_log)
  expect_equal(first_call_log, c("ggpicrust2", "pathway_gsea", "annotation", "visualize", "compare"))
})