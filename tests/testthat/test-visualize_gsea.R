# Test for visualize_gsea function
# create_gsea_test_results() is defined in helper-gsea.R

test_that("visualize_gsea creates different plot types correctly", {
  skip_if_not_installed("enrichplot")

  gsea_results <- create_gsea_test_results()

  p_bar <- visualize_gsea(gsea_results, plot_type = "barplot")
  expect_s3_class(p_bar, "ggplot")

  p_dot <- visualize_gsea(gsea_results, plot_type = "dotplot")
  expect_s3_class(p_dot, "ggplot")

  p_enrich <- visualize_gsea(gsea_results, plot_type = "enrichment_plot")
  expect_s3_class(p_enrich, "ggplot")
})

test_that("visualize_gsea creates network plot correctly", {
  skip_if_not_installed("igraph")

  gsea_results <- create_gsea_test_results()
  p <- visualize_gsea(gsea_results, plot_type = "network")
  expect_s3_class(p, "ggplot")
})

test_that("visualize_gsea creates heatmap plot correctly", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  gsea_results <- create_gsea_test_results()

  set.seed(123)
  abundance <- as.data.frame(matrix(
    runif(100 * 10, 1, 100), nrow = 100, ncol = 10,
    dimnames = list(paste0("K", sprintf("%05d", 1:100)), paste0("Sample", 1:10))
  ))

  metadata <- data.frame(
    sample = paste0("Sample", 1:10),
    group = rep(c("Group1", "Group2"), each = 5),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample

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

  gsea_results <- create_gsea_test_results()

  expect_error(visualize_gsea(gsea_results = "invalid"), "'gsea_results' must be a data frame")
  expect_error(visualize_gsea(gsea_results, plot_type = "invalid"), "plot_type must be one of")
  expect_error(visualize_gsea(gsea_results, sort_by = "invalid"), "sort_by must be one of")
})

test_that("visualize_gsea limits pathways correctly", {
  skip_if_not_installed("enrichplot")

  gsea_results <- create_gsea_test_results(n_pathways = 30)

  p_default <- visualize_gsea(gsea_results, plot_type = "barplot")
  expect_equal(nrow(p_default$data), 20)

  p_custom <- visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 10)
  expect_equal(nrow(p_custom$data), 10)
})
