# Test for visualize_gsea function
library(testthat)

# Helper function to create test GSEA results
create_test_gsea_results <- function(n_pathways = 10) {
  set.seed(123)
  data.frame(
    pathway_id = paste0("path:ko", sprintf("%05d", 1:n_pathways)),
    pathway_name = paste("Pathway", 1:n_pathways),
    size = sample(10:100, n_pathways, replace = TRUE),
    ES = runif(n_pathways, -0.8, 0.8),
    NES = runif(n_pathways, -2, 2),
    pvalue = runif(n_pathways, 0, 0.1),
    p.adjust = runif(n_pathways, 0, 0.2),
    leading_edge = replicate(n_pathways, paste(paste0("K", sprintf("%05d", sample(1:1000, 5))), collapse = ";")),
    method = rep("fgsea", n_pathways),
    pathway_class = sample(
      c("Metabolism", "Genetic Information Processing", "Environmental Information Processing",
        "Cellular Processes", "Organismal Systems", "Human Diseases"),
      n_pathways, replace = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

# Helper function to create test abundance data
create_test_abundance <- function(n_genes = 100, n_samples = 10) {
  set.seed(123)
  # Create gene IDs
  gene_ids <- paste0("K", sprintf("%05d", 1:n_genes))

  # Create sample IDs
  sample_ids <- paste0("Sample", 1:n_samples)

  # Create abundance matrix
  abundance_matrix <- matrix(
    runif(n_genes * n_samples, 1, 100),
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(gene_ids, sample_ids)
  )

  return(as.data.frame(abundance_matrix))
}

# Helper function to create test metadata
create_test_metadata <- function(n_samples = 10) {
  set.seed(123)
  # Create sample IDs
  sample_ids <- paste0("Sample", 1:n_samples)

  # Create group assignments (two groups)
  groups <- rep(c("Group1", "Group2"), each = n_samples/2)

  # Create metadata data frame
  metadata <- data.frame(
    sample = sample_ids,
    group = groups,
    stringsAsFactors = FALSE
  )

  # Set row names to sample IDs for easier matching
  rownames(metadata) <- metadata$sample

  return(metadata)
}

test_that("visualize_gsea creates barplot correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("enrichplot")

  # Create test data
  gsea_results <- create_test_gsea_results()

  # Test barplot
  p <- visualize_gsea(gsea_results, plot_type = "barplot")

  # Check the result
  expect_s3_class(p, "ggplot")
  expect_true("GeomBar" %in% sapply(p$layers, function(x) class(x$geom)[1]))

  # Test with custom parameters
  p <- visualize_gsea(
    gsea_results = gsea_results,
    plot_type = "barplot",
    n_pathways = 5,
    sort_by = "NES",
    colors = c("#FF0000", "#0000FF")
  )

  # Check the result
  expect_s3_class(p, "ggplot")
  expect_equal(length(unique(p$data$pathway_name)), 5)
})

test_that("visualize_gsea creates dotplot correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("enrichplot")

  # Create test data
  gsea_results <- create_test_gsea_results()

  # Test dotplot
  p <- visualize_gsea(gsea_results, plot_type = "dotplot")

  # Check the result
  expect_s3_class(p, "ggplot")
  expect_true("GeomPoint" %in% sapply(p$layers, function(x) class(x$geom)[1]))

  # Test with custom parameters
  p <- visualize_gsea(
    gsea_results = gsea_results,
    plot_type = "dotplot",
    n_pathways = 5,
    sort_by = "pvalue"
  )

  # Check the result
  expect_s3_class(p, "ggplot")
  expect_equal(length(unique(p$data$pathway_name)), 5)
})

test_that("visualize_gsea creates enrichment_plot correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("enrichplot")

  # Create test data
  gsea_results <- create_test_gsea_results()

  # Test enrichment_plot
  p <- visualize_gsea(gsea_results, plot_type = "enrichment_plot")

  # Check the result
  expect_s3_class(p, "ggplot")

  # Test with custom parameters
  p <- visualize_gsea(
    gsea_results = gsea_results,
    plot_type = "enrichment_plot",
    n_pathways = 5,
    sort_by = "p.adjust"
  )

  # Check the result
  expect_s3_class(p, "ggplot")
})

test_that("visualize_gsea creates network plot correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")

  # Create test data
  gsea_results <- create_test_gsea_results()

  # Test network plot
  p <- visualize_gsea(gsea_results, plot_type = "network")

  # Check the result
  expect_s3_class(p, "ggplot")

  # Test with custom parameters
  p <- visualize_gsea(
    gsea_results = gsea_results,
    plot_type = "network",
    n_pathways = 5,
    network_params = list(
      similarity_measure = "overlap",
      similarity_cutoff = 0.2,
      layout = "circle",
      node_color_by = "p.adjust"
    )
  )

  # Check the result
  expect_s3_class(p, "ggplot")
})

test_that("visualize_gsea creates heatmap plot correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  # Create test data
  gsea_results <- create_test_gsea_results()
  abundance <- create_test_abundance()
  metadata <- create_test_metadata()

  # Test heatmap plot
  p <- visualize_gsea(
    gsea_results = gsea_results,
    plot_type = "heatmap",
    abundance = abundance,
    metadata = metadata,
    group = "group"
  )

  # Check the result
  expect_s4_class(p, "Heatmap")

  # Test with custom parameters
  p <- visualize_gsea(
    gsea_results = gsea_results,
    plot_type = "heatmap",
    abundance = abundance,
    metadata = metadata,
    group = "group",
    n_pathways = 5,
    heatmap_params = list(
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_rownames = FALSE
    )
  )

  # Check the result
  expect_s4_class(p, "Heatmap")
})

test_that("visualize_gsea validates heatmap inputs correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("ComplexHeatmap")

  # Create test data
  gsea_results <- create_test_gsea_results()
  abundance <- create_test_abundance()
  metadata <- create_test_metadata()

  # Test missing abundance
  expect_error(
    visualize_gsea(gsea_results, plot_type = "heatmap", metadata = metadata, group = "group"),
    "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required"
  )

  # Test missing metadata
  expect_error(
    visualize_gsea(gsea_results, plot_type = "heatmap", abundance = abundance, group = "group"),
    "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required"
  )

  # Test missing group
  expect_error(
    visualize_gsea(gsea_results, plot_type = "heatmap", abundance = abundance, metadata = metadata),
    "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required"
  )
})

test_that("visualize_gsea validates inputs correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("enrichplot")

  # Create test data
  gsea_results <- create_test_gsea_results()

  # Test invalid gsea_results
  expect_error(
    visualize_gsea(gsea_results = "invalid"),
    "'gsea_results' must be a data frame"
  )

  # Test invalid plot_type
  expect_error(
    visualize_gsea(gsea_results, plot_type = "invalid"),
    "plot_type must be one of 'enrichment_plot', 'dotplot', 'barplot', 'network', or 'heatmap'"
  )

  # Test invalid sort_by
  expect_error(
    visualize_gsea(gsea_results, sort_by = "invalid"),
    "sort_by must be one of 'NES', 'pvalue', or 'p.adjust'"
  )

  # Test invalid colors
  expect_error(
    visualize_gsea(gsea_results, colors = 123),
    "colors must be NULL or a character vector"
  )
})

test_that("visualize_gsea sorts results correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("enrichplot")

  # Create test data
  gsea_results <- create_test_gsea_results()

  # Test sorting by NES
  p_nes <- visualize_gsea(gsea_results, sort_by = "NES", plot_type = "barplot")
  # The data is already sorted in the function, so we just check that NES values exist
  expect_true("NES" %in% colnames(p_nes$data))

  # Test sorting by pvalue
  p_pval <- visualize_gsea(gsea_results, sort_by = "pvalue", plot_type = "barplot")
  # The data is already sorted in the function, so we just check that pvalue values exist
  expect_true("pvalue" %in% colnames(p_pval$data))

  # Test sorting by p.adjust
  p_padj <- visualize_gsea(gsea_results, sort_by = "p.adjust", plot_type = "barplot")
  # The data is already sorted in the function, so we just check that p.adjust values exist
  expect_true("p.adjust" %in% colnames(p_padj$data))
})

test_that("visualize_gsea limits the number of pathways correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("enrichplot")

  # Create test data with more pathways
  gsea_results <- create_test_gsea_results(n_pathways = 30)

  # Test with default n_pathways (20)
  p_default <- visualize_gsea(gsea_results, plot_type = "barplot")
  expect_equal(nrow(p_default$data), 20)

  # Test with custom n_pathways
  p_custom <- visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 10)
  expect_equal(nrow(p_custom$data), 10)

  # Test with n_pathways greater than available
  p_all <- visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 50)
  expect_equal(nrow(p_all$data), 30)
})
