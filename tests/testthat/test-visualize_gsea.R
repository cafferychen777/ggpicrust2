# Test for visualize_gsea function

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
    pathway_class = sample(c("Metabolism", "Genetic Information Processing", "Cellular Processes"), n_pathways, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# Helper for test abundance and metadata
create_test_abundance <- function(n_genes = 100, n_samples = 10) {
  set.seed(123)
  abundance_matrix <- matrix(runif(n_genes * n_samples, 1, 100), nrow = n_genes, ncol = n_samples,
    dimnames = list(paste0("K", sprintf("%05d", 1:n_genes)), paste0("Sample", 1:n_samples)))
  as.data.frame(abundance_matrix)
}

create_test_metadata <- function(n_samples = 10) {
  metadata <- data.frame(
    sample = paste0("Sample", 1:n_samples),
    group = rep(c("Group1", "Group2"), each = n_samples/2),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample
  metadata
}

test_that("visualize_gsea creates different plot types correctly", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("enrichplot")

  gsea_results <- create_test_gsea_results()

  # Test barplot
  p_bar <- visualize_gsea(gsea_results, plot_type = "barplot")
  expect_s3_class(p_bar, "ggplot")

  # Test dotplot
  p_dot <- visualize_gsea(gsea_results, plot_type = "dotplot")
  expect_s3_class(p_dot, "ggplot")

  # Test enrichment_plot
  p_enrich <- visualize_gsea(gsea_results, plot_type = "enrichment_plot")
  expect_s3_class(p_enrich, "ggplot")
})

test_that("visualize_gsea creates network plot correctly", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")

  gsea_results <- create_test_gsea_results()
  p <- visualize_gsea(gsea_results, plot_type = "network")
  expect_s3_class(p, "ggplot")
})

test_that("visualize_gsea creates heatmap plot correctly", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  gsea_results <- create_test_gsea_results()
  abundance <- create_test_abundance()
  metadata <- create_test_metadata()

  p <- visualize_gsea(
    gsea_results = gsea_results,
    plot_type = "heatmap",
    abundance = abundance,
    metadata = metadata,
    group = "group"
  )
  expect_s4_class(p, "Heatmap")
})

test_that("visualize_gsea validates inputs correctly", {
  skip_if_not_installed("enrichplot")

  gsea_results <- create_test_gsea_results()

  expect_error(visualize_gsea(gsea_results = "invalid"), "'gsea_results' must be a data frame")
  expect_error(visualize_gsea(gsea_results, plot_type = "invalid"), "plot_type must be one of")
  expect_error(visualize_gsea(gsea_results, sort_by = "invalid"), "sort_by must be one of")
})

test_that("visualize_gsea limits pathways correctly", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("enrichplot")

  gsea_results <- create_test_gsea_results(n_pathways = 30)

  p_default <- visualize_gsea(gsea_results, plot_type = "barplot")
  expect_equal(nrow(p_default$data), 20)

  p_custom <- visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 10)
  expect_equal(nrow(p_custom$data), 10)
})
