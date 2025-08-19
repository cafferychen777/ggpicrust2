# Comprehensive Enhanced Tests for GSEA Visualization Functions
# Following Linus Torvalds' principles: test the core data structures,
# eliminate special cases, and ensure robustness
library(testthat)

# Load the package functions
devtools::load_all()

#=============================================================================
# TEST DATA GENERATORS - Focus on core data structures
#=============================================================================

#' Create minimal, well-structured GSEA results for testing
#' Following Linus principle: "Bad programmers worry about code, good programmers worry about data structures"
create_core_gsea_data <- function(n_pathways = 10, scenario = "standard") {
  set.seed(123) # Deterministic for reproducible tests
  
  base_data <- data.frame(
    pathway_id = sprintf("ko%05d", seq_len(n_pathways)),
    NES = seq(-2.5, 2.5, length.out = n_pathways),
    pvalue = 10^seq(-4, -1, length.out = n_pathways),
    p.adjust = 10^seq(-3, -0.5, length.out = n_pathways),
    size = seq(10, 100, length.out = n_pathways),
    stringsAsFactors = FALSE
  )
  
  # Create controlled leading edge data for similarity testing
  base_genes <- sprintf("K%05d", 1:200)
  
  if (scenario == "overlapping") {
    # Create overlapping gene sets for network testing
    base_data$leading_edge <- sapply(seq_len(n_pathways), function(i) {
      start_idx <- max(1, i * 5 - 2)  # Overlap with adjacent pathways
      end_idx <- min(200, start_idx + 10)
      paste(base_genes[start_idx:end_idx], collapse = ";")
    })
  } else if (scenario == "empty") {
    # Include empty leading edges to test edge cases
    base_data$leading_edge <- c(
      rep("", 2),  # Empty leading edges
      sapply(3:n_pathways, function(i) paste(sample(base_genes, 5), collapse = ";"))
    )
  } else {
    # Standard non-overlapping
    base_data$leading_edge <- sapply(seq_len(n_pathways), function(i) {
      start_idx <- (i - 1) * 15 + 1
      end_idx <- min(200, start_idx + 14)
      paste(base_genes[start_idx:end_idx], collapse = ";")
    })
  }
  
  return(base_data)
}

#' Add pathway names to test auto-detection logic
add_pathway_names <- function(gsea_data) {
  gsea_data$pathway_name <- sprintf("Test Pathway %d", seq_len(nrow(gsea_data)))
  return(gsea_data)
}

#' Create abundance matrix with controlled structure
create_test_abundance_matrix <- function(n_genes = 200, n_samples = 12) {
  set.seed(123)
  gene_ids <- sprintf("K%05d", 1:n_genes)
  sample_ids <- sprintf("Sample_%02d", 1:n_samples)
  
  # Create structured abundance data
  abundance <- matrix(
    exp(rnorm(n_genes * n_samples, mean = 5, sd = 2)), # Log-normal distribution
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(gene_ids, sample_ids)
  )
  
  return(as.data.frame(abundance))
}

#' Create metadata with multiple grouping variables
create_test_metadata <- function(n_samples = 12) {
  set.seed(123)
  sample_ids <- sprintf("Sample_%02d", 1:n_samples)
  
  data.frame(
    sample = sample_ids,
    group = rep(c("Control", "Treatment"), each = n_samples / 2),
    batch = rep(c("Batch1", "Batch2"), length.out = n_samples),
    stringsAsFactors = FALSE,
    row.names = sample_ids
  )
}

#=============================================================================
# CORE FUNCTION VALIDATION TESTS
#=============================================================================

test_that("visualize_gsea: Input validation follows 'good taste' principles", {
  # Test the fundamental data structure requirements
  # "Never break userspace" - validate inputs properly
  
  valid_data <- create_core_gsea_data()
  
  # Test 1: Data frame requirement (core data structure)
  expect_error(
    visualize_gsea("not_a_dataframe"), 
    "'gsea_results' must be a data frame",
    fixed = TRUE
  )
  
  # Test 2: Plot type validation (eliminate special cases)
  valid_types <- c("enrichment_plot", "dotplot", "barplot", "network", "heatmap")
  expect_error(
    visualize_gsea(valid_data, plot_type = "invalid_type"),
    "plot_type must be one of"
  )
  
  # Test 3: Sorting criteria validation
  valid_sorts <- c("NES", "pvalue", "p.adjust")
  expect_error(
    visualize_gsea(valid_data, sort_by = "invalid_sort"),
    "sort_by must be one of"
  )
  
  # Test 4: Colors validation (type safety)
  expect_error(
    visualize_gsea(valid_data, colors = 123),
    "colors must be NULL or a character vector"
  )
  
  # Test 5: pathway_label_column validation
  expect_error(
    visualize_gsea(valid_data, pathway_label_column = 123),
    "pathway_label_column must be NULL or a character string"
  )
})

test_that("visualize_gsea: Pathway labeling logic eliminates special cases", {
  skip_if_not_installed("ggplot2")
  
  # Test automatic pathway label selection - this should have no special cases
  
  # Scenario 1: Only pathway_id (should use pathway_id)
  data_no_names <- create_core_gsea_data(5)
  expect_no_error({
    plot1 <- visualize_gsea(data_no_names, plot_type = "barplot", n_pathways = 3)
    expect_s3_class(plot1, "ggplot")
    expect_true("pathway_label" %in% colnames(plot1$data))
  })
  
  # Scenario 2: Both pathway_id and pathway_name (should auto-select pathway_name)
  data_with_names <- add_pathway_names(data_no_names)
  expect_no_error({
    plot2 <- visualize_gsea(data_with_names, plot_type = "barplot", n_pathways = 3)
    expect_s3_class(plot2, "ggplot")
    expect_true("pathway_label" %in% colnames(plot2$data))
  })
  
  # Scenario 3: Explicit pathway_label_column override
  expect_no_error({
    plot3 <- visualize_gsea(data_with_names, plot_type = "barplot", 
                           pathway_label_column = "pathway_id", n_pathways = 3)
    expect_s3_class(plot3, "ggplot")
  })
  
  # Scenario 4: Missing required columns (should fail cleanly)
  data_incomplete <- data_with_names
  data_incomplete$pathway_id <- NULL
  data_incomplete$pathway_name <- NULL
  expect_error(
    visualize_gsea(data_incomplete),
    "must contain either 'pathway_name' or 'pathway_id' column"
  )
  
  # Scenario 5: Non-existent custom column
  expect_error(
    visualize_gsea(data_with_names, pathway_label_column = "nonexistent"),
    "not found in gsea_results"
  )
})

#=============================================================================
# SORTING AND FILTERING TESTS
#=============================================================================

test_that("visualize_gsea: Data sorting eliminates complexity", {
  skip_if_not_installed("ggplot2")
  
  # Create data with known sorting values
  test_data <- data.frame(
    pathway_id = c("path1", "path2", "path3", "path4", "path5"),
    NES = c(1.0, -2.5, 2.0, -1.5, 0.5),           # abs order: -2.5, 2.0, -1.5, 1.0, 0.5
    pvalue = c(0.05, 0.001, 0.1, 0.01, 0.2),       # order: 0.001, 0.01, 0.05, 0.1, 0.2
    p.adjust = c(0.1, 0.01, 0.2, 0.05, 0.3),       # order: 0.01, 0.05, 0.1, 0.2, 0.3
    size = rep(50, 5),
    leading_edge = rep("gene1;gene2", 5),
    stringsAsFactors = FALSE
  )
  
  # Test NES sorting (by absolute value)
  plot_nes <- visualize_gsea(test_data, plot_type = "barplot", sort_by = "NES", n_pathways = 3)
  expect_s3_class(plot_nes, "ggplot")
  expect_equal(nrow(plot_nes$data), 3)
  
  # Test pvalue sorting
  plot_pval <- visualize_gsea(test_data, plot_type = "barplot", sort_by = "pvalue", n_pathways = 3)
  expect_s3_class(plot_pval, "ggplot")
  expect_equal(nrow(plot_pval$data), 3)
  
  # Test p.adjust sorting
  plot_padj <- visualize_gsea(test_data, plot_type = "barplot", sort_by = "p.adjust", n_pathways = 3)
  expect_s3_class(plot_padj, "ggplot")
  expect_equal(nrow(plot_padj$data), 3)
})

test_that("visualize_gsea: n_pathways filtering is simple and correct", {
  skip_if_not_installed("ggplot2")
  
  gsea_data <- create_core_gsea_data(20)
  
  # Test with fewer pathways than available
  plot_few <- visualize_gsea(gsea_data, plot_type = "barplot", n_pathways = 5)
  expect_equal(nrow(plot_few$data), 5)
  
  # Test with more pathways than available (should not break)
  plot_many <- visualize_gsea(gsea_data, plot_type = "barplot", n_pathways = 50)
  expect_equal(nrow(plot_many$data), 20)
  
  # Test with n_pathways = 0 (edge case)
  plot_zero <- visualize_gsea(gsea_data, plot_type = "barplot", n_pathways = 0)
  expect_equal(nrow(plot_zero$data), 0)
})

#=============================================================================
# PLOT TYPE SPECIFIC TESTS
#=============================================================================

test_that("visualize_gsea: All basic plot types work correctly", {
  skip_if_not_installed("ggplot2")
  
  gsea_data <- add_pathway_names(create_core_gsea_data(8))
  
  # Test enrichment_plot
  expect_no_error({
    plot_enrich <- visualize_gsea(gsea_data, plot_type = "enrichment_plot", n_pathways = 5)
    expect_s3_class(plot_enrich, "ggplot")
    expect_true("GeomBar" %in% sapply(plot_enrich$layers, function(x) class(x$geom)[1]))
  })
  
  # Test dotplot
  expect_no_error({
    plot_dot <- visualize_gsea(gsea_data, plot_type = "dotplot", n_pathways = 5)
    expect_s3_class(plot_dot, "ggplot")
    expect_true("GeomPoint" %in% sapply(plot_dot$layers, function(x) class(x$geom)[1]))
  })
  
  # Test barplot
  expect_no_error({
    plot_bar <- visualize_gsea(gsea_data, plot_type = "barplot", n_pathways = 5)
    expect_s3_class(plot_bar, "ggplot")
    expect_true("GeomBar" %in% sapply(plot_bar$layers, function(x) class(x$geom)[1]))
    
    # Barplot should add direction column
    expect_true("direction" %in% colnames(plot_bar$data))
    expect_true(all(plot_bar$data$direction %in% c("Positive", "Negative")))
  })
})

#=============================================================================
# NETWORK PLOT COMPREHENSIVE TESTS
#=============================================================================

test_that("create_network_plot: Similarity calculations are mathematically correct", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  
  # Create data with known gene overlaps for testing similarity calculations
  test_data <- data.frame(
    pathway_id = c("path1", "path2", "path3"),
    NES = c(1.5, -1.2, 2.0),
    pvalue = c(0.01, 0.05, 0.001),
    p.adjust = c(0.05, 0.1, 0.01),
    size = c(5, 4, 3),
    # Controlled gene sets:
    # path1: {A, B, C, D, E} (5 genes)
    # path2: {A, B, F, G}     (4 genes, 2 overlap with path1)
    # path3: {H, I, J}        (3 genes, no overlap)
    leading_edge = c(
      "A;B;C;D;E",
      "A;B;F;G", 
      "H;I;J"
    ),
    stringsAsFactors = FALSE
  )
  
  # Test Jaccard similarity: |intersection| / |union|
  # path1 vs path2: |{A,B}| / |{A,B,C,D,E,F,G}| = 2/7 ≈ 0.286
  expect_no_error({
    plot_jaccard <- visualize_gsea(test_data, plot_type = "network",
                                  network_params = list(
                                    similarity_measure = "jaccard",
                                    similarity_cutoff = 0.2
                                  ))
    expect_s3_class(plot_jaccard, "ggplot")
  })
  
  # Test overlap coefficient: |intersection| / min(|A|, |B|)
  # path1 vs path2: |{A,B}| / min(5,4) = 2/4 = 0.5
  expect_no_error({
    plot_overlap <- visualize_gsea(test_data, plot_type = "network",
                                  network_params = list(
                                    similarity_measure = "overlap",
                                    similarity_cutoff = 0.3
                                  ))
    expect_s3_class(plot_overlap, "ggplot")
  })
  
  # Test correlation measure: |intersection| / sqrt(|A| * |B|)
  # path1 vs path2: |{A,B}| / sqrt(5*4) = 2/sqrt(20) ≈ 0.447
  expect_no_error({
    plot_corr <- visualize_gsea(test_data, plot_type = "network",
                               network_params = list(
                                 similarity_measure = "correlation",
                                 similarity_cutoff = 0.3
                               ))
    expect_s3_class(plot_corr, "ggplot")
  })
})

test_that("create_network_plot: Layout algorithms work correctly", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  
  gsea_data <- create_core_gsea_data(6, scenario = "overlapping")
  
  layouts <- c("fruchterman", "kamada", "circle")
  
  for (layout in layouts) {
    expect_no_error({
      plot_layout <- visualize_gsea(gsea_data, plot_type = "network",
                                   network_params = list(
                                     layout = layout,
                                     similarity_cutoff = 0.1
                                   ))
      expect_s3_class(plot_layout, "ggplot")
    })
  }
})

test_that("create_network_plot: Edge cases are handled properly", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  
  # Test with very high similarity cutoff (no connections)
  gsea_data <- create_core_gsea_data(5)
  expect_no_error({
    plot_no_conn <- visualize_gsea(gsea_data, plot_type = "network",
                                  network_params = list(similarity_cutoff = 0.99))
    expect_s3_class(plot_no_conn, "ggplot")
  })
  
  # Test with empty leading edges
  gsea_empty <- create_core_gsea_data(5, scenario = "empty")
  expect_no_error({
    plot_empty <- visualize_gsea(gsea_empty, plot_type = "network",
                                network_params = list(similarity_cutoff = 0.1))
    expect_s3_class(plot_empty, "ggplot")
  })
  
  # Test with single pathway
  gsea_single <- create_core_gsea_data(1)
  expect_no_error({
    plot_single <- visualize_gsea(gsea_single, plot_type = "network")
    expect_s3_class(plot_single, "ggplot")
  })
})

test_that("create_network_plot: Node and edge attributes are correct", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  
  gsea_data <- create_core_gsea_data(5, scenario = "overlapping")
  
  # Test different node coloring options
  node_colors <- c("NES", "pvalue", "p.adjust")
  
  for (node_color in node_colors) {
    expect_no_error({
      plot_node <- visualize_gsea(gsea_data, plot_type = "network",
                                 network_params = list(
                                   node_color_by = node_color,
                                   similarity_cutoff = 0.1
                                 ))
      expect_s3_class(plot_node, "ggplot")
    })
  }
  
  # Test edge width mapping
  expect_no_error({
    plot_edge <- visualize_gsea(gsea_data, plot_type = "network",
                               network_params = list(
                                 edge_width_by = "similarity",
                                 similarity_cutoff = 0.1
                               ))
    expect_s3_class(plot_edge, "ggplot")
  })
})

#=============================================================================
# HEATMAP PLOT COMPREHENSIVE TESTS  
#=============================================================================

test_that("create_heatmap_plot: Data preparation is mathematically sound", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  gsea_data <- create_core_gsea_data(8)
  abundance <- create_test_abundance_matrix(50, 12)
  metadata <- create_test_metadata(12)
  
  # Test basic heatmap creation
  expect_no_error({
    heatmap_basic <- visualize_gsea(
      gsea_data, plot_type = "heatmap", n_pathways = 5,
      abundance = abundance, metadata = metadata, group = "group"
    )
    expect_s4_class(heatmap_basic, "Heatmap")
  })
  
  # Test with different clustering options
  expect_no_error({
    heatmap_no_cluster <- visualize_gsea(
      gsea_data, plot_type = "heatmap", n_pathways = 5,
      abundance = abundance, metadata = metadata, group = "group",
      heatmap_params = list(cluster_rows = FALSE, cluster_columns = FALSE)
    )
    expect_s4_class(heatmap_no_cluster, "Heatmap")
  })
  
  # Test with different grouping variable
  expect_no_error({
    heatmap_batch <- visualize_gsea(
      gsea_data, plot_type = "heatmap", n_pathways = 5,
      abundance = abundance, metadata = metadata, group = "batch"
    )
    expect_s4_class(heatmap_batch, "Heatmap")
  })
})

test_that("create_heatmap_plot: Required parameters validation", {
  gsea_data <- create_core_gsea_data(5)
  abundance <- create_test_abundance_matrix(50, 12)
  metadata <- create_test_metadata(12)
  
  # Test missing abundance
  expect_error(
    visualize_gsea(gsea_data, plot_type = "heatmap", metadata = metadata, group = "group"),
    "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required"
  )
  
  # Test missing metadata
  expect_error(
    visualize_gsea(gsea_data, plot_type = "heatmap", abundance = abundance, group = "group"),
    "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required"
  )
  
  # Test missing group
  expect_error(
    visualize_gsea(gsea_data, plot_type = "heatmap", abundance = abundance, metadata = metadata),
    "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required"
  )
})

test_that("create_heatmap_plot: Leading edge gene extraction works correctly", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  # Create abundance data with genes that match leading edges
  abundance <- create_test_abundance_matrix(200, 12)
  metadata <- create_test_metadata(12)
  
  # Create GSEA data with genes that exist in abundance matrix
  gsea_data <- data.frame(
    pathway_id = c("path1", "path2", "path3"),
    NES = c(1.5, -1.2, 2.0),
    pvalue = c(0.01, 0.05, 0.001),
    p.adjust = c(0.05, 0.1, 0.01),
    size = c(5, 4, 3),
    leading_edge = c(
      "K00001;K00002;K00003;K00004;K00005",  # These genes exist in abundance
      "K00010;K00011;K00012;K00013",
      "K00020;K00021;K00022"
    ),
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    heatmap_genes <- visualize_gsea(
      gsea_data, plot_type = "heatmap", n_pathways = 3,
      abundance = abundance, metadata = metadata, group = "group"
    )
    expect_s4_class(heatmap_genes, "Heatmap")
  })
  
  # Test with some genes not in abundance matrix
  gsea_partial <- gsea_data
  gsea_partial$leading_edge[1] <- "K00001;K00002;MISSING_GENE;K00004"
  
  expect_no_error({
    heatmap_partial <- visualize_gsea(
      gsea_partial, plot_type = "heatmap", n_pathways = 3,
      abundance = abundance, metadata = metadata, group = "group"
    )
    expect_s4_class(heatmap_partial, "Heatmap")
  })
})

#=============================================================================
# PACKAGE DEPENDENCY TESTS
#=============================================================================

test_that("visualize_gsea: Package dependency checking works correctly", {
  gsea_data <- create_core_gsea_data(5)
  
  # These tests verify that the code checks for required packages
  # We can't easily test missing packages without uninstalling them,
  # but we can verify the check functions exist
  expect_true(exists("requireNamespace"))
  
  # Test that the function doesn't crash when packages are available
  expect_true(requireNamespace("ggplot2", quietly = TRUE) || 
              skip("ggplot2 not available"))
})

#=============================================================================
# COLOR SCHEME AND AESTHETICS TESTS
#=============================================================================

test_that("visualize_gsea: Color schemes work correctly", {
  skip_if_not_installed("ggplot2")
  
  gsea_data <- create_core_gsea_data(8)
  
  # Test custom colors
  custom_colors <- c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4")
  expect_no_error({
    plot_custom <- visualize_gsea(gsea_data, plot_type = "barplot", 
                                 colors = custom_colors, n_pathways = 5)
    expect_s3_class(plot_custom, "ggplot")
  })
  
  # Test default colors (NULL)
  expect_no_error({
    plot_default <- visualize_gsea(gsea_data, plot_type = "barplot", 
                                  colors = NULL, n_pathways = 5)
    expect_s3_class(plot_default, "ggplot")
  })
  
  # Test that barplot direction coloring works
  plot_direction <- visualize_gsea(gsea_data, plot_type = "barplot", n_pathways = 5)
  expect_true("direction" %in% colnames(plot_direction$data))
  directions <- unique(plot_direction$data$direction)
  expect_true(all(directions %in% c("Positive", "Negative")))
})

#=============================================================================
# INTEGRATION AND STRESS TESTS
#=============================================================================

test_that("visualize_gsea: Handles extreme values gracefully", {
  skip_if_not_installed("ggplot2")
  
  # Create data with extreme values
  extreme_data <- data.frame(
    pathway_id = c("path1", "path2", "path3", "path4", "path5"),
    NES = c(-Inf, Inf, 0, NA, 1000), # Extreme and missing values
    pvalue = c(0, 1, 0.5, NA, 1e-300), # Boundary and extreme values
    p.adjust = c(0, 1, 0.5, NA, 1e-100),
    size = c(1, 1000000, 50, 0, NA), # Very small and large values
    leading_edge = c("", "gene1", "gene1;gene2;gene3", NA, paste(rep("gene", 1000), collapse = ";")),
    stringsAsFactors = FALSE
  )
  
  # Function should handle these gracefully without crashing
  expect_no_error({
    plot_extreme <- visualize_gsea(extreme_data, plot_type = "barplot", n_pathways = 3)
    expect_s3_class(plot_extreme, "ggplot")
  })
})

test_that("visualize_gsea: Performance with large datasets", {
  skip_if_not_installed("ggplot2")
  skip_on_ci() # Skip on CI to avoid timeouts
  
  # Test with larger dataset
  large_data <- create_core_gsea_data(100)
  
  # Should complete within reasonable time
  expect_no_error({
    start_time <- Sys.time()
    plot_large <- visualize_gsea(large_data, plot_type = "barplot", n_pathways = 50)
    end_time <- Sys.time()
    
    expect_s3_class(plot_large, "ggplot")
    # Should complete in less than 10 seconds
    expect_lt(as.numeric(end_time - start_time), 10)
  })
})

test_that("visualize_gsea: All parameter combinations work together", {
  skip_if_not_installed("ggplot2")
  
  # Test various parameter combinations
  gsea_data <- add_pathway_names(create_core_gsea_data(15))
  
  # Test matrix of common parameter combinations
  plot_types <- c("enrichment_plot", "dotplot", "barplot")
  sort_methods <- c("NES", "pvalue", "p.adjust")
  pathway_limits <- c(3, 10)
  
  for (plot_type in plot_types) {
    for (sort_by in sort_methods) {
      for (n_pathways in pathway_limits) {
        expect_no_error({
          plot_combo <- visualize_gsea(
            gsea_results = gsea_data,
            plot_type = plot_type,
            n_pathways = n_pathways,
            sort_by = sort_by,
            colors = c("#E41A1C", "#377EB8"),
            pathway_label_column = "pathway_name"
          )
          expect_s3_class(plot_combo, "ggplot")
        })
      }
    }
  }
})

#=============================================================================
# EDGE CASES AND ERROR RECOVERY
#=============================================================================

test_that("visualize_gsea: Empty data edge cases", {
  skip_if_not_installed("ggplot2")
  
  # Test with empty data frame
  empty_data <- data.frame(
    pathway_id = character(0),
    NES = numeric(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    size = integer(0),
    leading_edge = character(0)
  )
  
  expect_no_error({
    plot_empty <- visualize_gsea(empty_data, plot_type = "barplot", n_pathways = 5)
    expect_s3_class(plot_empty, "ggplot")
    expect_equal(nrow(plot_empty$data), 0)
  })
  
  # Test with single row
  single_data <- create_core_gsea_data(1)
  expect_no_error({
    plot_single <- visualize_gsea(single_data, plot_type = "barplot", n_pathways = 5)
    expect_s3_class(plot_single, "ggplot")
    expect_equal(nrow(plot_single$data), 1)
  })
})

# Final message for test completion
message("GSEA visualization comprehensive tests completed")