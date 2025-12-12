# Detailed Tests for Heatmap Visualization Functions
library(testthat)

# Helper function to create realistic abundance data with known structure
create_structured_abundance <- function(n_genes = 50, n_samples = 12, add_batch_effects = TRUE) {
  set.seed(42)
  gene_ids <- paste0("K", sprintf("%05d", 1:n_genes))
  sample_ids <- paste0("Sample", 1:n_samples)
  
  # Create base abundance matrix
  abundance_matrix <- matrix(
    rgamma(n_genes * n_samples, shape = 2, rate = 0.02),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(gene_ids, sample_ids)
  )
  
  # Add group-specific patterns (simulate differential abundance)
  group1_samples <- 1:(n_samples/2)
  group2_samples <- (n_samples/2 + 1):n_samples
  
  # Make some genes higher in group1
  high_in_group1 <- 1:(n_genes/4)
  abundance_matrix[high_in_group1, group1_samples] <- 
    abundance_matrix[high_in_group1, group1_samples] * 2
  
  # Make some genes higher in group2
  high_in_group2 <- (n_genes/2 + 1):(3*n_genes/4)
  abundance_matrix[high_in_group2, group2_samples] <- 
    abundance_matrix[high_in_group2, group2_samples] * 2
  
  # Add batch effects if requested
  if (add_batch_effects) {
    batch1_samples <- c(1, 3, 5, 7, 9, 11)
    batch2_samples <- c(2, 4, 6, 8, 10, 12)
    
    # Add subtle batch effect to some genes
    batch_genes <- (n_genes/4 + 1):(n_genes/2)
    abundance_matrix[batch_genes, batch1_samples] <- 
      abundance_matrix[batch_genes, batch1_samples] * 1.2
  }
  
  # Add some zeros
  zero_prop <- 0.05
  zero_indices <- sample(length(abundance_matrix), size = floor(length(abundance_matrix) * zero_prop))
  abundance_matrix[zero_indices] <- 0
  
  return(as.data.frame(abundance_matrix))
}

# Helper function to create structured metadata
create_structured_metadata <- function(n_samples = 12) {
  set.seed(42)
  sample_ids <- paste0("Sample", 1:n_samples)
  
  metadata <- data.frame(
    sample = sample_ids,
    group = rep(c("Control", "Treatment"), each = n_samples/2),
    batch = rep(c("Batch1", "Batch2"), length.out = n_samples),
    subject = paste0("Subject", rep(1:(n_samples/2), each = 2)),
    timepoint = rep(c("T1", "T2"), n_samples/2),
    stringsAsFactors = FALSE
  )
  
  rownames(metadata) <- metadata$sample
  return(metadata)
}

# Helper function to create GSEA results with known gene sets
create_heatmap_test_gsea <- function(abundance_genes, n_pathways = 8) {
  set.seed(42)
  
  # Create pathway results that use genes from our abundance data
  pathway_genes <- list()
  for (i in 1:n_pathways) {
    # Select 5-15 genes for each pathway from our abundance genes
    n_genes_in_pathway <- sample(5:15, 1)
    pathway_genes[[i]] <- sample(abundance_genes, n_genes_in_pathway)
  }
  
  data.frame(
    pathway_id = paste0("path:ko", sprintf("%05d", 1:n_pathways)),
    pathway_name = paste("Pathway", letters[1:n_pathways]),
    NES = runif(n_pathways, -2.5, 2.5),
    pvalue = runif(n_pathways, 0.001, 0.1),
    p.adjust = runif(n_pathways, 0.001, 0.2),
    size = sapply(pathway_genes, length),
    leading_edge = sapply(pathway_genes, paste, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

test_that("create_heatmap_plot handles basic heatmap creation", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  abundance <- create_structured_abundance()
  metadata <- create_structured_metadata()
  gsea_results <- create_heatmap_test_gsea(rownames(abundance))
  
  expect_no_error({
    heatmap <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group",
      n_pathways = 5
    )
    expect_s4_class(heatmap, "Heatmap")
  })
})

test_that("create_heatmap_plot handles different clustering options", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  abundance <- create_structured_abundance()
  metadata <- create_structured_metadata()
  gsea_results <- create_heatmap_test_gsea(rownames(abundance))
  
  # Test no clustering
  expect_no_error({
    heatmap_no_cluster <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group",
      heatmap_params = list(
        cluster_rows = FALSE,
        cluster_columns = FALSE
      )
    )
    expect_s4_class(heatmap_no_cluster, "Heatmap")
  })
  
  # Test row clustering only
  expect_no_error({
    heatmap_row_cluster <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group",
      heatmap_params = list(
        cluster_rows = TRUE,
        cluster_columns = FALSE
      )
    )
    expect_s4_class(heatmap_row_cluster, "Heatmap")
  })
  
  # Test column clustering only
  expect_no_error({
    heatmap_col_cluster <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group",
      heatmap_params = list(
        cluster_rows = FALSE,
        cluster_columns = TRUE
      )
    )
    expect_s4_class(heatmap_col_cluster, "Heatmap")
  })
})

test_that("create_heatmap_plot handles different grouping variables", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  abundance <- create_structured_abundance()
  metadata <- create_structured_metadata()
  gsea_results <- create_heatmap_test_gsea(rownames(abundance))
  
  # Test different grouping variables
  group_vars <- c("group", "batch", "subject", "timepoint")
  
  for (group_var in group_vars) {
    # expect_no_error doesn't accept info parameter
    result <- tryCatch({
      heatmap <- visualize_gsea(
        gsea_results,
        plot_type = "heatmap",
        abundance = abundance,
        metadata = metadata,
        group = group_var,
        n_pathways = 5
      )
      expect_s4_class(heatmap, "Heatmap")
      TRUE
    }, error = function(e) {
      fail(paste("Group variable:", group_var, "Error:", e$message))
      FALSE
    })
    expect_true(result)
  }
})

test_that("create_heatmap_plot handles custom annotation colors", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  abundance <- create_structured_abundance()
  metadata <- create_structured_metadata()
  gsea_results <- create_heatmap_test_gsea(rownames(abundance))
  
  # Test custom annotation colors
  custom_colors <- list(
    Group = c("Control" = "#FF6B6B", "Treatment" = "#4ECDC4")
  )
  
  expect_no_error({
    heatmap_custom <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group",
      heatmap_params = list(
        annotation_colors = custom_colors
      )
    )
    expect_s4_class(heatmap_custom, "Heatmap")
  })
})

test_that("create_heatmap_plot handles show_rownames parameter", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  abundance <- create_structured_abundance()
  metadata <- create_structured_metadata()
  gsea_results <- create_heatmap_test_gsea(rownames(abundance))
  
  # Test with rownames hidden
  expect_no_error({
    heatmap_no_names <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group",
      heatmap_params = list(show_rownames = FALSE)
    )
    expect_s4_class(heatmap_no_names, "Heatmap")
  })
  
  # Test with rownames shown
  expect_no_error({
    heatmap_with_names <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group",
      heatmap_params = list(show_rownames = TRUE)
    )
    expect_s4_class(heatmap_with_names, "Heatmap")
  })
})

test_that("create_heatmap_plot handles missing genes in abundance data", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  abundance <- create_structured_abundance()
  metadata <- create_structured_metadata()
  
  # Create GSEA results with some genes not in abundance data
  gsea_results <- data.frame(
    pathway_id = c("path1", "path2", "path3"),
    pathway_name = c("Pathway A", "Pathway B", "Pathway C"),
    NES = c(2.1, -1.5, 1.8),
    pvalue = c(0.001, 0.05, 0.01),
    p.adjust = c(0.01, 0.1, 0.05),
    size = c(10, 8, 12),
    leading_edge = c(
      paste(c(rownames(abundance)[1:5], "MISSING1", "MISSING2"), collapse = ";"),
      paste(rownames(abundance)[6:10], collapse = ";"),
      paste(c(rownames(abundance)[11:15], "MISSING3"), collapse = ";")
    ),
    stringsAsFactors = FALSE
  )
  
  # Should handle missing genes gracefully
  expect_no_error({
    heatmap <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group"
    )
    expect_s4_class(heatmap, "Heatmap")
  })
})

test_that("create_heatmap_plot handles empty leading edges", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  abundance <- create_structured_abundance()
  metadata <- create_structured_metadata()
  
  # Create GSEA results with empty leading edges
  gsea_results <- data.frame(
    pathway_id = c("path1", "path2", "path3"),
    pathway_name = c("Pathway A", "Pathway B", "Pathway C"),
    NES = c(2.1, -1.5, 1.8),
    pvalue = c(0.001, 0.05, 0.01),
    p.adjust = c(0.01, 0.1, 0.05),
    size = c(10, 8, 12),
    leading_edge = c(
      paste(rownames(abundance)[1:5], collapse = ";"),
      "", # Empty leading edge
      paste(rownames(abundance)[6:10], collapse = ";")
    ),
    stringsAsFactors = FALSE
  )
  
  # Should handle empty leading edges gracefully
  expect_no_error({
    heatmap <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group"
    )
    expect_s4_class(heatmap, "Heatmap")
  })
})

test_that("create_heatmap_plot handles single pathway", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  abundance <- create_structured_abundance()
  metadata <- create_structured_metadata()
  
  # Create single pathway GSEA results
  gsea_single <- data.frame(
    pathway_id = "path1",
    pathway_name = "Single Pathway",
    NES = 2.1,
    pvalue = 0.001,
    p.adjust = 0.01,
    size = 10,
    leading_edge = paste(rownames(abundance)[1:10], collapse = ";"),
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    heatmap <- visualize_gsea(
      gsea_single,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group"
    )
    expect_s4_class(heatmap, "Heatmap")
  })
})

test_that("create_heatmap_plot handles abundance-metadata sample mismatch", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  abundance <- create_structured_abundance()
  metadata <- create_structured_metadata()
  gsea_results <- create_heatmap_test_gsea(rownames(abundance))
  
  # Create metadata with different sample names
  metadata_mismatch <- metadata
  rownames(metadata_mismatch) <- paste0("DifferentSample", 1:nrow(metadata_mismatch))
  metadata_mismatch$sample <- rownames(metadata_mismatch)
  
  # Should handle sample mismatch (though results might be limited)
  expect_no_error({
    heatmap <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata_mismatch,
      group = "group"
    )
    expect_s4_class(heatmap, "Heatmap")
  })
})

test_that("create_heatmap_plot validates required parameters", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  abundance <- create_structured_abundance()
  metadata <- create_structured_metadata()
  gsea_results <- create_heatmap_test_gsea(rownames(abundance))
  
  # Test missing abundance parameter
  expect_error(
    visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      metadata = metadata,
      group = "group"
    ),
    "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required"
  )
  
  # Test missing metadata parameter
  expect_error(
    visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      group = "group"
    ),
    "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required"
  )
  
  # Test missing group parameter
  expect_error(
    visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata
    ),
    "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required"
  )
})

test_that("create_heatmap_plot handles different pathway label columns", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  abundance <- create_structured_abundance()
  metadata <- create_structured_metadata()
  gsea_results <- create_heatmap_test_gsea(rownames(abundance))
  
  # Test with pathway_name
  expect_no_error({
    heatmap_names <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group",
      pathway_label_column = "pathway_name"
    )
    expect_s4_class(heatmap_names, "Heatmap")
  })
  
  # Test with pathway_id
  expect_no_error({
    heatmap_ids <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group",
      pathway_label_column = "pathway_id"
    )
    expect_s4_class(heatmap_ids, "Heatmap")
  })
})

test_that("create_heatmap_plot handles scaling edge cases", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  metadata <- create_structured_metadata()
  
  # Create abundance data with extreme values
  abundance_extreme <- create_structured_abundance()
  abundance_extreme[1, ] <- 1000  # Very high values
  abundance_extreme[2, ] <- 0.001 # Very low values
  
  gsea_results <- create_heatmap_test_gsea(rownames(abundance_extreme))
  
  expect_no_error({
    heatmap <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance_extreme,
      metadata = metadata,
      group = "group"
    )
    expect_s4_class(heatmap, "Heatmap")
  })
  
  # Create abundance data with all zeros for some genes
  abundance_zeros <- create_structured_abundance()
  abundance_zeros[1:5, ] <- 0
  
  gsea_results_zeros <- data.frame(
    pathway_id = "path1",
    pathway_name = "Zero Pathway",
    NES = 1.5,
    pvalue = 0.01,
    p.adjust = 0.05,
    size = 5,
    leading_edge = paste(rownames(abundance_zeros)[1:5], collapse = ";"),
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    heatmap_zeros <- visualize_gsea(
      gsea_results_zeros,
      plot_type = "heatmap",
      abundance = abundance_zeros,
      metadata = metadata,
      group = "group"
    )
    expect_s4_class(heatmap_zeros, "Heatmap")
  })
})

test_that("create_heatmap_plot handles large datasets", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  # Create larger dataset
  abundance_large <- create_structured_abundance(n_genes = 200, n_samples = 24)
  metadata_large <- create_structured_metadata(n_samples = 24)
  gsea_results_large <- create_heatmap_test_gsea(rownames(abundance_large), n_pathways = 15)
  
  expect_no_error({
    heatmap_large <- visualize_gsea(
      gsea_results_large,
      plot_type = "heatmap",
      abundance = abundance_large,
      metadata = metadata_large,
      group = "group",
      n_pathways = 10  # Limit pathways for performance
    )
    expect_s4_class(heatmap_large, "Heatmap")
  })
})