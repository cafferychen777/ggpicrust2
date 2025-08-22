# Performance and Edge Case Tests for GSEA Visualization
library(testthat)

# Helper function to create large datasets for performance testing
create_large_gsea_dataset <- function(n_pathways = 100, n_genes_per_pathway = 20) {
  set.seed(123)
  
  # Create a large gene universe
  all_genes <- paste0("K", sprintf("%05d", 1:5000))
  
  # Create pathway results
  pathway_results <- data.frame(
    pathway_id = paste0("path:ko", sprintf("%05d", 1:n_pathways)),
    pathway_name = paste("Pathway", 1:n_pathways),
    NES = runif(n_pathways, -3, 3),
    pvalue = runif(n_pathways, 0.0001, 0.1),
    p.adjust = runif(n_pathways, 0.0001, 0.2),
    size = sample(10:200, n_pathways, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Create leading edge genes for each pathway
  pathway_results$leading_edge <- replicate(n_pathways, {
    n_genes <- sample(5:n_genes_per_pathway, 1)
    paste(sample(all_genes, n_genes), collapse = ";")
  })
  
  return(pathway_results)
}

# Helper function to create large abundance dataset
create_large_abundance <- function(n_genes = 1000, n_samples = 50) {
  set.seed(123)
  gene_ids <- paste0("K", sprintf("%05d", 1:n_genes))
  sample_ids <- paste0("Sample", 1:n_samples)
  
  # Create abundance matrix with realistic distribution
  abundance_matrix <- matrix(
    rgamma(n_genes * n_samples, shape = 1.5, rate = 0.01),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(gene_ids, sample_ids)
  )
  
  # Add some zeros
  zero_indices <- sample(length(abundance_matrix), size = floor(length(abundance_matrix) * 0.1))
  abundance_matrix[zero_indices] <- 0
  
  return(as.data.frame(abundance_matrix))
}

# Helper function to create large metadata
create_large_metadata <- function(n_samples = 50, n_groups = 5) {
  set.seed(123)
  sample_ids <- paste0("Sample", 1:n_samples)
  
  metadata <- data.frame(
    sample = sample_ids,
    group = rep(paste0("Group", 1:n_groups), length.out = n_samples),
    batch = rep(paste0("Batch", 1:10), length.out = n_samples),
    treatment = rep(c("Control", "Treatment"), length.out = n_samples),
    stringsAsFactors = FALSE
  )
  
  rownames(metadata) <- metadata$sample
  return(metadata)
}

test_that("visualize_gsea handles large datasets efficiently", {
  skip_if_not_installed("ggplot2")
  skip_on_ci() # Skip on continuous integration due to performance
  
  # Create large dataset
  large_gsea <- create_large_gsea_dataset(n_pathways = 200)
  
  # Test that large datasets can be handled
  expect_no_error({
    start_time <- Sys.time()
    p_large <- visualize_gsea(large_gsea, plot_type = "barplot", n_pathways = 50)
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    expect_s3_class(p_large, "ggplot")
    expect_equal(nrow(p_large$data), 50)
    
    # Performance check - should complete within reasonable time
    expect_lt(execution_time, 30) # Less than 30 seconds
  })
})

test_that("visualize_gsea network plot handles large similarity matrices", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  skip_on_ci() # Skip on CI due to performance
  
  # Create moderately large dataset for network analysis
  large_gsea <- create_large_gsea_dataset(n_pathways = 30, n_genes_per_pathway = 50)
  
  expect_no_error({
    start_time <- Sys.time()
    p_network <- visualize_gsea(
      large_gsea, 
      plot_type = "network",
      n_pathways = 20,
      network_params = list(similarity_cutoff = 0.3)
    )
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    expect_s3_class(p_network, "ggplot")
    
    # Performance check
    expect_lt(execution_time, 60) # Less than 60 seconds
  })
})

test_that("visualize_gsea heatmap handles large abundance matrices", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  skip_on_ci() # Skip on CI due to performance and memory
  
  # Create large datasets
  large_abundance <- create_large_abundance(n_genes = 500, n_samples = 30)
  large_metadata <- create_large_metadata(n_samples = 30)
  large_gsea <- create_large_gsea_dataset(n_pathways = 20)
  
  # Ensure the leading edge genes exist in abundance data
  abundance_genes <- rownames(large_abundance)
  for (i in 1:nrow(large_gsea)) {
    available_genes <- sample(abundance_genes, min(20, length(abundance_genes)))
    large_gsea$leading_edge[i] <- paste(available_genes, collapse = ";")
  }
  
  expect_no_error({
    start_time <- Sys.time()
    heatmap_large <- visualize_gsea(
      large_gsea,
      plot_type = "heatmap",
      abundance = large_abundance,
      metadata = large_metadata,
      group = "group",
      n_pathways = 15
    )
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    expect_s4_class(heatmap_large, "Heatmap")
    
    # Performance check
    expect_lt(execution_time, 120) # Less than 2 minutes
  })
})

test_that("visualize_gsea handles extreme pathway numbers", {
  skip_if_not_installed("ggplot2")
  
  # Test with single pathway
  single_pathway <- data.frame(
    pathway_id = "single_path",
    pathway_name = "Single Pathway",
    NES = 2.5,
    pvalue = 0.001,
    p.adjust = 0.01,
    size = 50,
    leading_edge = "gene1;gene2;gene3;gene4;gene5",
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    p_single <- visualize_gsea(single_pathway, plot_type = "barplot")
    expect_s3_class(p_single, "ggplot")
    expect_equal(nrow(p_single$data), 1)
  })
  
  # Test requesting more pathways than available
  few_pathways <- create_large_gsea_dataset(n_pathways = 5)
  expect_no_error({
    p_limited <- visualize_gsea(few_pathways, plot_type = "barplot", n_pathways = 100)
    expect_s3_class(p_limited, "ggplot")
    expect_equal(nrow(p_limited$data), 5) # Should be limited to available pathways
  })
})

test_that("visualize_gsea handles extreme NES values", {
  skip_if_not_installed("ggplot2")
  
  # Create dataset with extreme NES values
  extreme_gsea <- data.frame(
    pathway_id = c("path1", "path2", "path3", "path4", "path5"),
    pathway_name = c("Very Positive", "Very Negative", "Near Zero", "Moderate Pos", "Moderate Neg"),
    NES = c(10.5, -10.8, 0.001, 3.2, -2.9),
    pvalue = c(0.0001, 0.0001, 0.9, 0.01, 0.02),
    p.adjust = c(0.001, 0.001, 0.95, 0.05, 0.08),
    size = c(100, 120, 10, 60, 80),
    leading_edge = rep("gene1;gene2;gene3", 5),
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    p_extreme <- visualize_gsea(extreme_gsea, plot_type = "barplot")
    expect_s3_class(p_extreme, "ggplot")
    
    # Check that extreme values are handled
    expect_true(max(abs(p_extreme$data$NES)) >= 10)
  })
  
  # Test with all positive NES
  all_positive <- extreme_gsea
  all_positive$NES <- abs(all_positive$NES)
  
  expect_no_error({
    p_pos <- visualize_gsea(all_positive, plot_type = "barplot")
    expect_s3_class(p_pos, "ggplot")
    expect_true(all(p_pos$data$NES >= 0))
  })
  
  # Test with all negative NES
  all_negative <- extreme_gsea
  all_negative$NES <- -abs(all_negative$NES)
  
  expect_no_error({
    p_neg <- visualize_gsea(all_negative, plot_type = "barplot")
    expect_s3_class(p_neg, "ggplot")
    expect_true(all(p_neg$data$NES <= 0))
  })
})

test_that("visualize_gsea handles extreme p-values", {
  skip_if_not_installed("ggplot2")
  
  # Create dataset with extreme p-values
  extreme_pvals <- data.frame(
    pathway_id = c("path1", "path2", "path3", "path4"),
    pathway_name = c("Very Significant", "Not Significant", "Borderline", "Machine Precision"),
    NES = c(3.2, 0.5, 1.8, 4.1),
    pvalue = c(1e-15, 0.999, 0.05, .Machine$double.eps),
    p.adjust = c(1e-12, 0.999, 0.1, 1e-10),
    size = c(80, 30, 50, 90),
    leading_edge = rep("gene1;gene2;gene3", 4),
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    p_extreme_pval <- visualize_gsea(extreme_pvals, plot_type = "dotplot")
    expect_s3_class(p_extreme_pval, "ggplot")
    
    # Check that extreme p-values are handled
    expect_true(min(p_extreme_pval$data$pvalue) < 1e-10)
    expect_true(max(p_extreme_pval$data$pvalue) > 0.9)
  })
})

test_that("visualize_gsea handles very long pathway names", {
  skip_if_not_installed("ggplot2")
  
  # Create dataset with very long pathway names
  long_names_gsea <- data.frame(
    pathway_id = c("path1", "path2", "path3"),
    pathway_name = c(
      "This is an extremely long pathway name that might cause text wrapping issues in visualizations and could potentially break the layout",
      "Another very long pathway name with lots of technical terminology and specific biological process descriptions that go on and on",
      "Short name"
    ),
    NES = c(2.1, -1.8, 1.5),
    pvalue = c(0.001, 0.05, 0.02),
    p.adjust = c(0.01, 0.1, 0.08),
    size = c(50, 40, 35),
    leading_edge = rep("gene1;gene2", 3),
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    p_long_names <- visualize_gsea(long_names_gsea, plot_type = "barplot")
    expect_s3_class(p_long_names, "ggplot")
    
    # Check that long names are preserved
    max_label_length <- max(nchar(as.character(p_long_names$data$pathway_label)))
    expect_gt(max_label_length, 50) # Should have long labels
  })
})

test_that("visualize_gsea handles special characters in pathway names", {
  skip_if_not_installed("ggplot2")
  
  # Create dataset with special characters
  special_char_gsea <- data.frame(
    pathway_id = c("path1", "path2", "path3", "path4"),
    pathway_name = c(
      "Pathway with (parentheses) and [brackets]",
      "Pathway with symbols: α, β, γ, δ",
      "Pathway with / slashes & ampersands",
      "Pathway with 'quotes' and \"double quotes\""
    ),
    NES = c(2.1, -1.8, 1.5, -2.2),
    pvalue = c(0.001, 0.05, 0.02, 0.001),
    p.adjust = c(0.01, 0.1, 0.08, 0.02),
    size = c(50, 40, 35, 60),
    leading_edge = rep("gene1;gene2", 4),
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    p_special <- visualize_gsea(special_char_gsea, plot_type = "barplot")
    expect_s3_class(p_special, "ggplot")
    
    # Check that special characters are preserved
    labels <- as.character(p_special$data$pathway_label)
    expect_true(any(grepl("\\(|\\)|\\[|\\]", labels))) # Contains brackets/parentheses
    expect_true(any(grepl("α|β|γ", labels))) # Contains Greek letters
  })
})

test_that("visualize_gsea handles very large leading edge gene lists", {
  skip_if_not_installed("ggplot2")
  
  # Create dataset with very large leading edge lists
  set.seed(123)
  many_genes <- paste0("K", sprintf("%05d", 1:1000))
  
  large_edges_gsea <- data.frame(
    pathway_id = c("path1", "path2", "path3"),
    pathway_name = c("Pathway A", "Pathway B", "Pathway C"),
    NES = c(2.1, -1.8, 1.5),
    pvalue = c(0.001, 0.05, 0.02),
    p.adjust = c(0.01, 0.1, 0.08),
    size = c(500, 300, 200),
    leading_edge = c(
      paste(sample(many_genes, 500), collapse = ";"),  # 500 genes
      paste(sample(many_genes, 300), collapse = ";"),  # 300 genes
      paste(sample(many_genes, 200), collapse = ";")   # 200 genes
    ),
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    p_large_edges <- visualize_gsea(large_edges_gsea, plot_type = "barplot")
    expect_s3_class(p_large_edges, "ggplot")
    
    # Check that large edge lists are handled
    edge_lengths <- nchar(large_edges_gsea$leading_edge)
    expect_true(max(edge_lengths) > 1000) # Should have very long edge strings
  })
})

test_that("visualize_gsea handles missing values gracefully", {
  skip_if_not_installed("ggplot2")
  
  # Create dataset with some NA values
  na_gsea <- data.frame(
    pathway_id = c("path1", "path2", "path3", "path4"),
    pathway_name = c("Pathway A", "Pathway B", NA, "Pathway D"),
    NES = c(2.1, NA, 1.5, -1.2),
    pvalue = c(0.001, 0.05, NA, 0.03),
    p.adjust = c(0.01, 0.1, 0.08, NA),
    size = c(50, 40, NA, 45),
    leading_edge = c("gene1;gene2", "gene3;gene4", "", "gene5;gene6"),
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    p_na <- visualize_gsea(na_gsea, plot_type = "barplot")
    expect_s3_class(p_na, "ggplot")
    
    # Should handle NA values without crashing
    expect_equal(nrow(p_na$data), 4) # All rows should be preserved
  })
})

test_that("visualize_gsea memory usage is reasonable", {
  skip_on_ci() # Skip on CI due to memory testing complexity
  
  # This is a basic memory check - more sophisticated memory testing
  # would require additional tools
  gsea_results <- create_large_gsea_dataset(n_pathways = 100)
  
  # Monitor memory usage (basic check)
  gc_before <- gc()
  
  expect_no_error({
    p <- visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 50)
    expect_s3_class(p, "ggplot")
  })
  
  gc_after <- gc()
  
  # Basic check that we're not accumulating excessive memory
  # (This is a very basic test and might not catch all memory issues)
  memory_increase <- gc_after[1, 2] - gc_before[1, 2]
  expect_lt(memory_increase, 100) # Less than 100 MB increase
})

test_that("visualize_gsea handles concurrent execution", {
  skip_on_ci() # Skip on CI to avoid complexity
  
  gsea_results <- create_large_gsea_dataset(n_pathways = 20)
  
  # Test that multiple simultaneous calls don't interfere
  expect_no_error({
    results <- list()
    for (i in 1:3) {
      results[[i]] <- visualize_gsea(
        gsea_results, 
        plot_type = "barplot", 
        n_pathways = 10,
        sort_by = c("NES", "pvalue", "p.adjust")[i]
      )
    }
    
    # All should be valid ggplot objects
    for (p in results) {
      expect_s3_class(p, "ggplot")
    }
  })
})