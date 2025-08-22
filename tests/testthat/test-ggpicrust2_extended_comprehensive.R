# Comprehensive tests for ggpicrust2_extended function integration workflows
library(testthat)

# Advanced test helpers for creating realistic data
create_comprehensive_test_data <- function(n_samples = 20, n_features = 50, groups = c("Control", "Treatment")) {
  set.seed(12345)
  
  # Create abundance matrix with realistic KO identifiers
  ko_ids <- paste0("K", sprintf("%05d", sample(1:20000, n_features)))
  sample_ids <- paste0("Sample", 1:n_samples)
  
  # Generate abundance data with group effects
  group_assignment <- rep(groups, length.out = n_samples)
  abundance_matrix <- matrix(0, nrow = n_features, ncol = n_samples)
  
  for (i in 1:n_features) {
    # Base abundance
    base_abundance <- rpois(n_samples, lambda = 50)
    
    # Add group effect for some features
    if (runif(1) < 0.3) {  # 30% of features have group differences
      effect_size <- rnorm(1, mean = 0, sd = 1.5)
      for (j in 1:n_samples) {
        if (group_assignment[j] == "Treatment") {
          base_abundance[j] <- rpois(1, lambda = pmax(1, base_abundance[j] * exp(effect_size)))
        }
      }
    }
    abundance_matrix[i, ] <- base_abundance
  }
  
  rownames(abundance_matrix) <- ko_ids
  colnames(abundance_matrix) <- sample_ids
  
  # Convert to data frame format expected by ggpicrust2
  abundance_df <- data.frame(
    `#NAME` = rownames(abundance_matrix),
    abundance_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Create metadata
  metadata <- data.frame(
    sample_id = sample_ids,
    Group = factor(group_assignment),
    Age = sample(20:80, n_samples, replace = TRUE),
    BMI = rnorm(n_samples, mean = 25, sd = 5),
    stringsAsFactors = FALSE
  )
  
  return(list(
    abundance = abundance_df,
    metadata = metadata,
    true_group_assignment = group_assignment
  ))
}

create_mock_ggpicrust2_results <- function(features) {
  n_features <- length(features)
  
  list(
    list(
      plot = ggplot2::ggplot() + ggplot2::geom_blank(),
      results = data.frame(
        feature = features,
        method = rep("LinDA", n_features),
        group1 = rep("Control", n_features),
        group2 = rep("Treatment", n_features),
        p_values = runif(n_features, 0.001, 0.2),
        p_adjust = runif(n_features, 0.001, 0.25),
        log_2_fold_change = rnorm(n_features, 0, 2),
        effect_size = rnorm(n_features, 0, 1),
        stringsAsFactors = FALSE
      )
    )
  )
}

create_mock_gsea_results <- function(pathway_ids) {
  n_pathways <- length(pathway_ids)
  
  data.frame(
    pathway_id = pathway_ids,
    pathway_name = paste("Pathway", pathway_ids),
    size = sample(10:200, n_pathways, replace = TRUE),
    ES = rnorm(n_pathways, 0, 0.6),
    NES = rnorm(n_pathways, 0, 1.8),
    pvalue = runif(n_pathways, 0.001, 0.15),
    p.adjust = runif(n_pathways, 0.001, 0.2),
    leading_edge = replicate(n_pathways, {
      genes <- paste0("K", sprintf("%05d", sample(1:20000, sample(3:15, 1))))
      paste(genes, collapse = ";")
    }),
    method = rep("fgsea", n_pathways),
    stringsAsFactors = FALSE
  )
}

test_that("ggpicrust2_extended integrates with standard ggpicrust2 workflow correctly", {
  skip_if_not_installed("dplyr")
  
  # Create comprehensive test data
  test_data <- create_comprehensive_test_data(n_samples = 16, n_features = 30)
  
  # Mock ggpicrust2 with realistic results
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    create_mock_ggpicrust2_results(test_data$abundance$`#NAME`)
  })
  
  # Test basic integration without GSEA
  result <- ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Group",
    pathway = "KO",
    daa_method = "LinDA",
    run_gsea = FALSE
  )
  
  # Verify integration structure
  expect_type(result, "list")
  expect_true("daa_results" %in% names(result))
  expect_false("gsea_results" %in% names(result))
  expect_false("gsea_plot" %in% names(result))
  expect_false("comparison" %in% names(result))
  
  # Verify DAA results structure matches ggpicrust2 output
  expect_s3_class(result$daa_results[[1]]$plot, "ggplot")
  expect_true(is.data.frame(result$daa_results[[1]]$results))
  expect_true(all(c("feature", "p_adjust", "log_2_fold_change") %in% colnames(result$daa_results[[1]]$results)))
})

test_that("ggpicrust2_extended handles GSEA parameter passing correctly", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  # Create test data
  test_data <- create_comprehensive_test_data(n_samples = 20, n_features = 40)
  
  # Mock all dependencies
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    create_mock_ggpicrust2_results(test_data$abundance$`#NAME`)
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    args <- list(...)
    
    # Verify that GSEA parameters are passed correctly
    expect_true("abundance" %in% names(args))
    expect_true("metadata" %in% names(args))
    expect_true("group" %in% names(args))
    expect_true("pathway_type" %in% names(args))
    expect_equal(args$method, "fgsea")
    expect_equal(args$rank_method, "log2FC")
    expect_equal(args$nperm, 5000)
    expect_equal(args$min_size, 5)
    expect_equal(args$max_size, 800)
    
    # Return mock results
    create_mock_gsea_results(sample(test_data$abundance$`#NAME`, 25))
  })
  
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(gsea_results, ...) {
    gsea_results$pathway_class <- sample(c("Metabolism", "Genetic Information Processing", "Environmental Information Processing"), 
                                        nrow(gsea_results), replace = TRUE)
    return(gsea_results)
  })
  
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(...) {
    ggplot2::ggplot() + ggplot2::geom_blank()
  })
  
  mockery::stub(ggpicrust2_extended, "compare_gsea_daa", function(...) {
    list(
      plot = ggplot2::ggplot() + ggplot2::geom_blank(),
      results = list(
        overlap = c("K00001", "K00002"),
        gsea_only = c("K00003", "K00004"),
        daa_only = c("K00005", "K00006"),
        n_overlap = 2,
        n_gsea_only = 2,
        n_daa_only = 2,
        n_gsea_total = 4,
        n_daa_total = 4
      )
    )
  })
  
  # Test with custom GSEA parameters
  result <- ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Group",
    pathway = "KO",
    daa_method = "LinDA",
    run_gsea = TRUE,
    gsea_params = list(
      method = "fgsea",
      rank_method = "log2FC",
      nperm = 5000,
      min_size = 5,
      max_size = 800,
      p.adjust = "bonferroni"
    )
  )
  
  # Verify complete workflow results
  expect_type(result, "list")
  expect_true(all(c("daa_results", "gsea_results", "gsea_plot", "comparison") %in% names(result)))
})

test_that("ggpicrust2_extended handles run_gsea flag correctly in different scenarios", {
  skip_if_not_installed("dplyr")
  
  test_data <- create_comprehensive_test_data(n_samples = 12, n_features = 20)
  
  # Mock ggpicrust2
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    create_mock_ggpicrust2_results(test_data$abundance$`#NAME`)
  })
  
  # Test run_gsea = FALSE (default)
  result_no_gsea <- ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Group",
    pathway = "KO",
    daa_method = "ALDEx2"
  )
  
  expect_false("gsea_results" %in% names(result_no_gsea))
  expect_false("gsea_plot" %in% names(result_no_gsea))
  expect_false("comparison" %in% names(result_no_gsea))
  
  # Test run_gsea = FALSE explicitly
  result_explicit_no_gsea <- ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Group",
    pathway = "KO",
    daa_method = "ALDEx2",
    run_gsea = FALSE
  )
  
  expect_identical(names(result_no_gsea), names(result_explicit_no_gsea))
})

test_that("ggpicrust2_extended creates consistent combined result structure", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  test_data <- create_comprehensive_test_data(n_samples = 18, n_features = 35)
  
  # Mock all functions with detailed tracking
  gsea_call_count <- 0
  annotation_call_count <- 0
  visualization_call_count <- 0
  comparison_call_count <- 0
  
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    create_mock_ggpicrust2_results(test_data$abundance$`#NAME`)
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    gsea_call_count <<- gsea_call_count + 1
    create_mock_gsea_results(sample(test_data$abundance$`#NAME`, 20))
  })
  
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(gsea_results, ...) {
    annotation_call_count <<- annotation_call_count + 1
    gsea_results$pathway_class <- sample(c("Metabolism", "Genetic Information Processing"), 
                                        nrow(gsea_results), replace = TRUE)
    return(gsea_results)
  })
  
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(...) {
    visualization_call_count <<- visualization_call_count + 1
    ggplot2::ggplot() + ggplot2::geom_blank() + ggplot2::ggtitle("GSEA Results")
  })
  
  mockery::stub(ggpicrust2_extended, "compare_gsea_daa", function(...) {
    comparison_call_count <<- comparison_call_count + 1
    list(
      plot = ggplot2::ggplot() + ggplot2::geom_blank() + ggplot2::ggtitle("Comparison"),
      results = list(
        overlap = c("K00001", "K00002", "K00003"),
        gsea_only = c("K00004", "K00005"),
        daa_only = c("K00006", "K00007"),
        n_overlap = 3,
        n_gsea_only = 2,
        n_daa_only = 2,
        n_gsea_total = 5,
        n_daa_total = 5
      )
    )
  })
  
  # Run full workflow
  result <- ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Group",
    pathway = "KO",
    daa_method = "LinDA",
    run_gsea = TRUE,
    gsea_params = list(method = "fgsea", nperm = 1000)
  )
  
  # Verify all functions were called exactly once
  expect_equal(gsea_call_count, 1)
  expect_equal(annotation_call_count, 1)
  expect_equal(visualization_call_count, 1)
  expect_equal(comparison_call_count, 1)
  
  # Verify result structure
  expect_type(result, "list")
  expect_setequal(names(result), c("daa_results", "gsea_results", "gsea_plot", "comparison"))
  
  # Verify each component structure
  expect_type(result$daa_results, "list")
  expect_true(is.data.frame(result$gsea_results))
  expect_s3_class(result$gsea_plot, "ggplot")
  expect_type(result$comparison, "list")
  expect_true("plot" %in% names(result$comparison))
  expect_true("results" %in% names(result$comparison))
})

test_that("ggpicrust2_extended handles workflow error propagation correctly", {
  skip_if_not_installed("dplyr")
  
  test_data <- create_comprehensive_test_data()
  
  # Test error in ggpicrust2 (should propagate immediately)
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    stop("ggpicrust2 failed with invalid input")
  })
  
  expect_error(
    ggpicrust2_extended(
      data = test_data$abundance,
      metadata = test_data$metadata,
      group = "Group",
      pathway = "KO",
      daa_method = "LinDA"
    ),
    "ggpicrust2 failed with invalid input"
  )
  
  # Test error in GSEA pathway (should handle gracefully)
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    create_mock_ggpicrust2_results(test_data$abundance$`#NAME`)
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    stop("GSEA computation failed")
  })
  
  # Error in GSEA should propagate when run_gsea = TRUE
  expect_error(
    ggpicrust2_extended(
      data = test_data$abundance,
      metadata = test_data$metadata,
      group = "Group",
      pathway = "KO",
      daa_method = "LinDA",
      run_gsea = TRUE
    ),
    "GSEA computation failed"
  )
})

test_that("ggpicrust2_extended validates input parameters for integrated workflow", {
  # Test missing abundance data when GSEA is requested
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    list(list(plot = ggplot2::ggplot(), results = data.frame()))
  })
  
  # Mock requireNamespace to return TRUE for all packages
  mockery::stub(ggpicrust2_extended, "requireNamespace", function(...) TRUE)
  
  expect_error(
    ggpicrust2_extended(
      metadata = data.frame(sample_id = 1:5, group = rep(c("A", "B"), length.out = 5)),
      group = "group",
      pathway = "KO",
      run_gsea = TRUE
    ),
    "No abundance data provided for GSEA analysis"
  )
  
  # Test invalid gsea_params
  test_data <- create_comprehensive_test_data(n_samples = 10, n_features = 15)
  
  # Should error when gsea_params is not a list
  expect_error(
    ggpicrust2_extended(
      data = test_data$abundance,
      metadata = test_data$metadata,
      group = "Group",
      pathway = "KO",
      run_gsea = TRUE,
      gsea_params = "invalid_params"
    )
  )
})

test_that("ggpicrust2_extended handles different pathway types correctly", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  test_data <- create_comprehensive_test_data(n_samples = 14, n_features = 25)
  
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    create_mock_ggpicrust2_results(test_data$abundance$`#NAME`)
  })
  
  # Track pathway_type parameter
  captured_pathway_types <- character()
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    args <- list(...)
    captured_pathway_types <<- c(captured_pathway_types, args$pathway_type)
    create_mock_gsea_results(sample(test_data$abundance$`#NAME`, 15))
  })
  
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(gsea_results, ...) gsea_results)
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(...) ggplot2::ggplot())
  mockery::stub(ggpicrust2_extended, "compare_gsea_daa", function(...) list(plot = ggplot2::ggplot(), results = list()))
  
  # Test KO pathway
  ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Group",
    pathway = "KO",
    run_gsea = TRUE
  )
  
  # Test KEGG pathway
  ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Group",
    pathway = "KEGG",
    run_gsea = TRUE
  )
  
  # Test MetaCyc pathway
  ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Group",
    pathway = "MetaCyc",
    run_gsea = TRUE
  )
  
  # Test EC pathway
  ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Group",
    pathway = "EC",
    run_gsea = TRUE
  )
  
  # Verify pathway types were set correctly
  expect_equal(captured_pathway_types, c("KEGG", "KEGG", "MetaCyc", "EC"))
})

test_that("ggpicrust2_extended handles data format conversion correctly", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  # Create test data in different formats
  test_data <- create_comprehensive_test_data(n_samples = 16, n_features = 30)
  
  # Convert to matrix format (should be handled)
  abundance_matrix <- as.matrix(test_data$abundance[, -1])
  rownames(abundance_matrix) <- test_data$abundance$`#NAME`
  
  captured_abundance_data <- NULL
  
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    create_mock_ggpicrust2_results(rownames(abundance_matrix))
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    args <- list(...)
    captured_abundance_data <<- args$abundance
    create_mock_gsea_results(sample(rownames(abundance_matrix), 15))
  })
  
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(gsea_results, ...) gsea_results)
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(...) ggplot2::ggplot())
  mockery::stub(ggpicrust2_extended, "compare_gsea_daa", function(...) list(plot = ggplot2::ggplot(), results = list()))
  
  # Test with data frame format
  result_df <- ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Group",
    pathway = "KO",
    run_gsea = TRUE
  )
  
  # Verify abundance data was processed correctly
  expect_true(is.matrix(captured_abundance_data) || is.data.frame(captured_abundance_data))
  expect_true(nrow(captured_abundance_data) > 0)
  expect_true(ncol(captured_abundance_data) > 0)
})

test_that("ggpicrust2_extended provides informative progress messages", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")
  
  test_data <- create_comprehensive_test_data(n_samples = 12, n_features = 20)
  
  # Mock all functions
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    create_mock_ggpicrust2_results(test_data$abundance$`#NAME`)
  })
  
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    create_mock_gsea_results(sample(test_data$abundance$`#NAME`, 10))
  })
  
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(gsea_results, ...) gsea_results)
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(...) ggplot2::ggplot())
  mockery::stub(ggpicrust2_extended, "compare_gsea_daa", function(...) list(plot = ggplot2::ggplot(), results = list()))
  
  # Capture messages
  expect_message(
    expect_message(
      expect_message(
        expect_message(
          ggpicrust2_extended(
            data = test_data$abundance,
            metadata = test_data$metadata,
            group = "Group",
            pathway = "KO",
            run_gsea = TRUE
          ),
          "Performing Gene Set Enrichment Analysis"
        ),
        "Annotating GSEA results"
      ),
      "Creating GSEA visualization"
    ),
    "Comparing GSEA and DAA results"
  )
})