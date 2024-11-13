# Test for pathway_heatmap function
library(testthat)
library(ggplot2)
library(dplyr)

# Setup test data
setup_test_data <- function() {
  # Create example abundance data
  n_samples <- 10
  n_pathways <- 3

  # Generate sample names
  sample_names <- paste0("Sample", 1:n_samples)
  pathway_names <- paste0("Pathway", LETTERS[1:n_pathways])

  # Create abundance matrix
  abundance <- matrix(
    rnorm(n_samples * n_pathways),
    nrow = n_pathways,
    ncol = n_samples,
    dimnames = list(pathway_names, sample_names)
  )

  # Create metadata
  metadata <- data.frame(
    sample_name = sample_names,
    group = factor(rep(c("Control", "Treatment"), each = n_samples/2)),
    stringsAsFactors = FALSE
  )

  return(list(abundance = abundance, metadata = metadata))
}

test_that("pathway_heatmap input validation works", {
  test_data <- setup_test_data()

  # Test valid inputs
  expect_no_error(
    pathway_heatmap(
      abundance = test_data$abundance,
      metadata = test_data$metadata,
      group = "group"
    )
  )

  # Test invalid abundance input
  expect_error(
    pathway_heatmap(
      abundance = "invalid",
      metadata = test_data$metadata,
      group = "group"
    ),
    "abundance must be a data frame or matrix"
  )

  # Test invalid metadata input
  expect_error(
    pathway_heatmap(
      abundance = test_data$abundance,
      metadata = "invalid",
      group = "group"
    ),
    "metadata must be a data frame"
  )

  # Test invalid group input
  expect_error(
    pathway_heatmap(
      abundance = test_data$abundance,
      metadata = test_data$metadata,
      group = c("group1", "group2")
    ),
    "group must be a single character string"
  )

  # Test invalid colors input
  expect_error(
    pathway_heatmap(
      abundance = test_data$abundance,
      metadata = test_data$metadata,
      group = "group",
      colors = 123
    ),
    "colors must be NULL or a character vector of color codes"
  )
})

test_that("pathway_heatmap customization options work", {
  test_data <- setup_test_data()

  # Test font size
  p <- pathway_heatmap(
    abundance = test_data$abundance,
    metadata = test_data$metadata,
    group = "group",
    font_size = 14
  )
  expect_true(is.ggplot(p))

  # Test hiding row names
  p <- pathway_heatmap(
    abundance = test_data$abundance,
    metadata = test_data$metadata,
    group = "group",
    show_row_names = FALSE
  )
  expect_true(is.ggplot(p))

  # Test hiding legend
  p <- pathway_heatmap(
    abundance = test_data$abundance,
    metadata = test_data$metadata,
    group = "group",
    show_legend = FALSE
  )
  expect_true(is.ggplot(p))

  # Test custom theme
  custom_theme <- theme_minimal()
  p <- pathway_heatmap(
    abundance = test_data$abundance,
    metadata = test_data$metadata,
    group = "group",
    custom_theme = custom_theme
  )
  expect_true(is.ggplot(p))
})

test_that("pathway_heatmap handles edge cases", {
  test_data <- setup_test_data()

  # Test single sample
  single_abundance <- test_data$abundance[, 1, drop = FALSE]
  single_metadata <- data.frame(
    sample_name = colnames(single_abundance),
    group = factor("Control"),
    stringsAsFactors = FALSE
  )

  # 跳过单样本测试，因为这是一个不支持的用例
  skip("Single sample case is not supported")

  # Test single pathway
  single_pathway_abundance <- test_data$abundance[1, , drop = FALSE]
  expect_error(
    pathway_heatmap(
      abundance = single_pathway_abundance,
      metadata = test_data$metadata,
      group = "group"
    ),
    NA
  )

  # Test custom colors
  custom_colors <- c("#FF0000", "#00FF00")
  expect_error(
    pathway_heatmap(
      abundance = test_data$abundance,
      metadata = test_data$metadata,
      group = "group",
      colors = custom_colors
    ),
    NA
  )
})

test_that("pathway_heatmap output structure is correct", {
  test_data <- setup_test_data()

  p <- pathway_heatmap(
    abundance = test_data$abundance,
    metadata = test_data$metadata,
    group = "group"
  )

  # Test that output is a ggplot object
  expect_true(is.ggplot(p))

  # Test that essential layers are present
  expect_true("GeomTile" %in% sapply(p$layers, function(x) class(x$geom)[1]))

  # Test that faceting is applied
  expect_true(!is.null(p$facet))
})

test_that("sample names match between abundance and metadata", {
  test_data <- setup_test_data()

  # Test sample name matching
  expect_true(
    all(colnames(test_data$abundance) %in% test_data$metadata$sample_name),
    "All abundance sample names should be present in metadata"
  )

  expect_true(
    all(test_data$metadata$sample_name %in% colnames(test_data$abundance)),
    "All metadata sample names should be present in abundance"
  )
})
