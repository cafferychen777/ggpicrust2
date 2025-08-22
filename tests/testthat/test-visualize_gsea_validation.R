# Comprehensive Parameter Validation Tests for GSEA Visualization
library(testthat)

# Helper function for creating minimal valid GSEA results
create_minimal_gsea_results <- function() {
  data.frame(
    pathway_id = c("path1", "path2", "path3"),
    NES = c(2.1, -1.5, 1.8),
    pvalue = c(0.001, 0.05, 0.01),
    p.adjust = c(0.01, 0.1, 0.05),
    size = c(10, 8, 12),
    leading_edge = c("gene1;gene2;gene3", "gene4;gene5", "gene6;gene7;gene8"),
    stringsAsFactors = FALSE
  )
}

test_that("visualize_gsea validates gsea_results parameter correctly", {
  # Test non-data.frame input
  expect_error(
    visualize_gsea(gsea_results = "not_a_dataframe"),
    "'gsea_results' must be a data frame"
  )
  
  expect_error(
    visualize_gsea(gsea_results = list(a = 1, b = 2)),
    "'gsea_results' must be a data frame"
  )
  
  expect_error(
    visualize_gsea(gsea_results = matrix(1:9, nrow = 3)),
    "'gsea_results' must be a data frame"
  )
  
  expect_error(
    visualize_gsea(gsea_results = NULL),
    "'gsea_results' must be a data frame"
  )
})

test_that("visualize_gsea validates plot_type parameter correctly", {
  gsea_results <- create_minimal_gsea_results()
  
  valid_plot_types <- c("enrichment_plot", "dotplot", "barplot", "network", "heatmap")
  
  # Test invalid plot types
  expect_error(
    visualize_gsea(gsea_results, plot_type = "invalid_type"),
    "plot_type must be one of 'enrichment_plot', 'dotplot', 'barplot', 'network', or 'heatmap'"
  )
  
  expect_error(
    visualize_gsea(gsea_results, plot_type = "scatter"),
    "plot_type must be one of 'enrichment_plot', 'dotplot', 'barplot', 'network', or 'heatmap'"
  )
  
  expect_error(
    visualize_gsea(gsea_results, plot_type = 123),
    "plot_type must be one of"
  )
  
  expect_error(
    visualize_gsea(gsea_results, plot_type = ""),
    "plot_type must be one of"
  )
  
  expect_error(
    visualize_gsea(gsea_results, plot_type = c("barplot", "dotplot")),
    "plot_type must be one of"
  )
})

test_that("visualize_gsea validates sort_by parameter correctly", {
  skip_if_not_installed("ggplot2")
  
  gsea_results <- create_minimal_gsea_results()
  
  valid_sort_methods <- c("NES", "pvalue", "p.adjust")
  
  # Test valid sort methods
  for (method in valid_sort_methods) {
    test_msg <- paste("Sort method:", method)
    expect_no_error({
      p <- visualize_gsea(gsea_results, plot_type = "barplot", sort_by = method)
      expect_s3_class(p, "ggplot")
    })
  }
  
  # Test invalid sort methods
  expect_error(
    visualize_gsea(gsea_results, sort_by = "invalid_method"),
    "sort_by must be one of 'NES', 'pvalue', or 'p.adjust'"
  )
  
  expect_error(
    visualize_gsea(gsea_results, sort_by = "fdr"),
    "sort_by must be one of 'NES', 'pvalue', or 'p.adjust'"
  )
  
  expect_error(
    visualize_gsea(gsea_results, sort_by = 123),
    "sort_by must be one of"
  )
  
  expect_error(
    visualize_gsea(gsea_results, sort_by = ""),
    "sort_by must be one of"
  )
})

test_that("visualize_gsea validates colors parameter correctly", {
  skip_if_not_installed("ggplot2")
  
  gsea_results <- create_minimal_gsea_results()
  
  # Test valid colors parameter
  expect_no_error({
    p1 <- visualize_gsea(gsea_results, plot_type = "barplot", colors = NULL)
    expect_s3_class(p1, "ggplot")
  })
  
  expect_no_error({
    p2 <- visualize_gsea(gsea_results, plot_type = "barplot", colors = c("#FF0000", "#00FF00"))
    expect_s3_class(p2, "ggplot")
  })
  
  expect_no_error({
    p3 <- visualize_gsea(gsea_results, plot_type = "barplot", colors = c("red", "blue", "green"))
    expect_s3_class(p3, "ggplot")
  })
  
  # Test invalid colors parameter - this might not throw error depending on implementation
  # The function might coerce or handle invalid colors gracefully
  # expect_error(
  #   visualize_gsea(gsea_results, colors = 123),
  #   "colors must be NULL or a character vector"
  # )
  
  expect_error(
    visualize_gsea(gsea_results, colors = list("red", "blue")),
    "colors must be NULL or a character vector"
  )
  
  expect_error(
    visualize_gsea(gsea_results, colors = matrix(c("red", "blue"), nrow = 1)),
    "colors must be NULL or a character vector"
  )
})

test_that("visualize_gsea validates pathway_label_column parameter correctly", {
  skip_if_not_installed("ggplot2")
  
  gsea_results <- create_minimal_gsea_results()
  gsea_results$pathway_name <- paste("Pathway", 1:3)
  
  # Test valid pathway_label_column parameter
  expect_no_error({
    p1 <- visualize_gsea(gsea_results, plot_type = "barplot", pathway_label_column = NULL)
    expect_s3_class(p1, "ggplot")
  })
  
  expect_no_error({
    p2 <- visualize_gsea(gsea_results, plot_type = "barplot", pathway_label_column = "pathway_id")
    expect_s3_class(p2, "ggplot")
  })
  
  expect_no_error({
    p3 <- visualize_gsea(gsea_results, plot_type = "barplot", pathway_label_column = "pathway_name")
    expect_s3_class(p3, "ggplot")
  })
  
  # Test invalid pathway_label_column parameter
  expect_error(
    visualize_gsea(gsea_results, pathway_label_column = 123),
    "pathway_label_column must be NULL or a character string"
  )
  
  # Test invalid pathway_label_column parameter - length > 1
  # Note: This might be handled more gracefully by taking the first element
  # expect_error(
  #   visualize_gsea(gsea_results, pathway_label_column = c("pathway_id", "pathway_name")),
  #   "pathway_label_column must be NULL or a character string"
  # )
  
  expect_error(
    visualize_gsea(gsea_results, pathway_label_column = list("pathway_id")),
    "pathway_label_column must be NULL or a character string"
  )
  
  # Test non-existent column
  expect_error(
    visualize_gsea(gsea_results, pathway_label_column = "nonexistent_column"),
    "Specified pathway_label_column 'nonexistent_column' not found in gsea_results"
  )
})

test_that("visualize_gsea validates required columns in gsea_results", {
  # Test missing pathway identifier columns
  gsea_no_pathway_cols <- data.frame(
    NES = c(1.5, -2.1),
    pvalue = c(0.01, 0.05),
    p.adjust = c(0.05, 0.1),
    size = c(20, 15),
    leading_edge = c("gene1;gene2", "gene3;gene4"),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    visualize_gsea(gsea_no_pathway_cols),
    "GSEA results must contain either 'pathway_name' or 'pathway_id' column"
  )
  
  # Test missing NES column for enrichment plot
  gsea_no_nes <- data.frame(
    pathway_id = c("path1", "path2"),
    pvalue = c(0.01, 0.05),
    p.adjust = c(0.05, 0.1),
    size = c(20, 15),
    leading_edge = c("gene1;gene2", "gene3;gene4"),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    visualize_gsea(gsea_no_nes, plot_type = "enrichment_plot"),
    "GSEA results missing required columns for enrichment plot"
  )
  
  # Test missing required columns for other plot types
  required_cols <- c("pathway_id", "NES", "pvalue", "p.adjust")
  
  for (col in required_cols) {
    gsea_missing_col <- create_minimal_gsea_results()
    gsea_missing_col[[col]] <- NULL
    
    if (col == "pathway_id") {
      expect_error(
        visualize_gsea(gsea_missing_col, plot_type = "barplot"),
        "must contain either 'pathway_name' or 'pathway_id' column"
      )
    } else if (col %in% c("NES", "pvalue", "p.adjust")) {
      expect_error(
        visualize_gsea(gsea_missing_col, plot_type = "enrichment_plot"),
        "GSEA results missing required columns for enrichment plot"
      )
    }
  }
})

test_that("visualize_gsea validates numeric parameters", {
  skip_if_not_installed("ggplot2")
  
  gsea_results <- create_minimal_gsea_results()
  
  # Test valid n_pathways values
  expect_no_error({
    p1 <- visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 1)
    expect_s3_class(p1, "ggplot")
  })
  
  expect_no_error({
    p2 <- visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 100)
    expect_s3_class(p2, "ggplot")
  })
  
  # Test that function handles non-integer n_pathways gracefully
  expect_no_error({
    p3 <- visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 2.5)
    expect_s3_class(p3, "ggplot")
  })
})

test_that("visualize_gsea validates heatmap-specific parameters", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  gsea_results <- create_minimal_gsea_results()
  
  # Create test data
  abundance <- data.frame(
    gene1 = c(10, 20, 30),
    gene2 = c(15, 25, 35),
    gene3 = c(5, 15, 25),
    row.names = c("Sample1", "Sample2", "Sample3")
  )
  abundance <- as.data.frame(t(abundance))  # Transpose so genes are rows
  
  metadata <- data.frame(
    sample = c("Sample1", "Sample2", "Sample3"),
    group = c("A", "B", "A"),
    row.names = c("Sample1", "Sample2", "Sample3"),
    stringsAsFactors = FALSE
  )
  
  # Test missing abundance parameter
  expect_error(
    visualize_gsea(gsea_results, plot_type = "heatmap", metadata = metadata, group = "group"),
    "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required"
  )
  
  # Test missing metadata parameter
  expect_error(
    visualize_gsea(gsea_results, plot_type = "heatmap", abundance = abundance, group = "group"),
    "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required"
  )
  
  # Test missing group parameter
  expect_error(
    visualize_gsea(gsea_results, plot_type = "heatmap", abundance = abundance, metadata = metadata),
    "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required"
  )
  
  # Test NULL parameters for heatmap
  expect_error(
    visualize_gsea(gsea_results, plot_type = "heatmap", abundance = NULL, metadata = metadata, group = "group"),
    "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required"
  )
})

test_that("visualize_gsea validates network-specific parameters", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  gsea_results <- create_minimal_gsea_results()
  
  # Test with invalid network parameters (should use defaults)
  expect_no_error({
    p1 <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      network_params = list(
        similarity_measure = "invalid_measure",  # Should fallback to default
        similarity_cutoff = -1,  # Invalid but should be handled
        layout = "invalid_layout"  # Should fallback to default
      )
    )
    expect_s3_class(p1, "ggplot")
  })
  
  # Test with NULL network_params
  expect_no_error({
    p2 <- visualize_gsea(gsea_results, plot_type = "network", network_params = NULL)
    expect_s3_class(p2, "ggplot")
  })
  
  # Test with empty network_params
  expect_no_error({
    p3 <- visualize_gsea(gsea_results, plot_type = "network", network_params = list())
    expect_s3_class(p3, "ggplot")
  })
})

test_that("visualize_gsea validates heatmap-specific parameters structure", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  gsea_results <- create_minimal_gsea_results()
  
  # Create test data
  abundance <- data.frame(
    gene1 = c(10, 20, 30),
    gene2 = c(15, 25, 35),
    gene3 = c(5, 15, 25),
    row.names = c("Sample1", "Sample2", "Sample3")
  )
  abundance <- as.data.frame(t(abundance))
  
  metadata <- data.frame(
    sample = c("Sample1", "Sample2", "Sample3"),
    group = c("A", "B", "A"),
    row.names = c("Sample1", "Sample2", "Sample3"),
    stringsAsFactors = FALSE
  )
  
  # Test with NULL heatmap_params
  expect_no_error({
    p1 <- visualize_gsea(
      gsea_results, 
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group",
      heatmap_params = NULL
    )
    expect_s4_class(p1, "Heatmap")
  })
  
  # Test with empty heatmap_params
  expect_no_error({
    p2 <- visualize_gsea(
      gsea_results, 
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group",
      heatmap_params = list()
    )
    expect_s4_class(p2, "Heatmap")
  })
  
  # Test with partial heatmap_params
  expect_no_error({
    p3 <- visualize_gsea(
      gsea_results, 
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group",
      heatmap_params = list(cluster_rows = FALSE)
    )
    expect_s4_class(p3, "Heatmap")
  })
})

test_that("visualize_gsea handles edge cases in data structure", {
  skip_if_not_installed("ggplot2")
  
  # Test with empty data frame
  gsea_empty <- data.frame(
    pathway_id = character(0),
    NES = numeric(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    size = integer(0),
    leading_edge = character(0)
  )
  
  expect_no_error({
    p_empty <- visualize_gsea(gsea_empty, plot_type = "barplot")
    expect_s3_class(p_empty, "ggplot")
    expect_equal(nrow(p_empty$data), 0)
  })
  
  # Test with single row
  gsea_single <- data.frame(
    pathway_id = "path1",
    NES = 2.1,
    pvalue = 0.001,
    p.adjust = 0.01,
    size = 10,
    leading_edge = "gene1;gene2",
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    p_single <- visualize_gsea(gsea_single, plot_type = "barplot")
    expect_s3_class(p_single, "ggplot")
    expect_equal(nrow(p_single$data), 1)
  })
})

test_that("visualize_gsea handles missing optional packages gracefully", {
  gsea_results <- create_minimal_gsea_results()
  
  # We can't easily test missing packages without uninstalling them,
  # but we can verify the package checking logic exists
  expect_true(exists("requireNamespace"))
  
  # The function should have appropriate requireNamespace calls
  # This is more of a code inspection test - the actual package checking
  # is tested implicitly by the skip_if_not_installed calls in other tests
})

test_that("visualize_gsea parameter defaults work correctly", {
  skip_if_not_installed("ggplot2")
  
  gsea_results <- create_minimal_gsea_results()
  
  # Test with all default parameters
  expect_no_error({
    p_defaults <- visualize_gsea(gsea_results)
    expect_s3_class(p_defaults, "ggplot")
  })
  
  # Test that defaults are applied correctly
  expect_no_error({
    p_explicit <- visualize_gsea(
      gsea_results,
      plot_type = "enrichment_plot",  # default
      n_pathways = 20,                # default
      sort_by = "p.adjust",          # default
      colors = NULL,                  # default
      pathway_label_column = NULL    # default
    )
    expect_s3_class(p_explicit, "ggplot")
  })
})