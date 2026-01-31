# Helper: prepare GSEA data for ridgeplot integration tests
create_ridgeplot_test_data <- function() {
  data(ko_abundance, package = "ggpicrust2")
  data(metadata, package = "ggpicrust2")

  abundance_data <- as.data.frame(ko_abundance)
  rownames(abundance_data) <- abundance_data[, "#NAME"]
  abundance_data <- abundance_data[, -1]

  set.seed(42)
  gsea_results <- suppressMessages(pathway_gsea(
    abundance = abundance_data,
    metadata = metadata,
    group = "Environment",
    pathway_type = "KEGG",
    method = "camera",
    min_size = 5
  ))

  list(gsea_results = gsea_results, abundance = abundance_data, metadata = metadata)
}

test_that("pathway_ridgeplot creates a ggplot object", {
  skip_if_not_installed("ggridges")
  skip_on_cran()

  td <- create_ridgeplot_test_data()

  p <- pathway_ridgeplot(
    gsea_results = td$gsea_results,
    abundance = td$abundance,
    metadata = td$metadata,
    group = "Environment",
    n_pathways = 5
  )

  expect_s3_class(p, "ggplot")
})

test_that("pathway_ridgeplot handles custom parameters", {
  skip_if_not_installed("ggridges")
  skip_on_cran()

  td <- create_ridgeplot_test_data()

  p <- pathway_ridgeplot(
    gsea_results = td$gsea_results,
    abundance = td$abundance,
    metadata = td$metadata,
    group = "Environment",
    n_pathways = 10,
    sort_by = "pvalue",
    show_direction = TRUE,
    colors = c("Down" = "blue", "Up" = "red"),
    title = "Custom Title",
    scale_height = 0.8,
    alpha = 0.5
  )

  expect_s3_class(p, "ggplot")
})

test_that("pathway_ridgeplot handles show_direction = FALSE", {
  skip_if_not_installed("ggridges")
  skip_on_cran()

  td <- create_ridgeplot_test_data()

  p <- pathway_ridgeplot(
    gsea_results = td$gsea_results,
    abundance = td$abundance,
    metadata = td$metadata,
    group = "Environment",
    n_pathways = 5,
    show_direction = FALSE
  )

  expect_s3_class(p, "ggplot")
})

test_that("pathway_ridgeplot errors on invalid input", {
  skip_if_not_installed("ggridges")

  expect_error(
    pathway_ridgeplot(
      gsea_results = "not a data frame",
      abundance = data.frame(),
      metadata = data.frame(),
      group = "group"
    ),
    "must be a data frame"
  )

  expect_error(
    pathway_ridgeplot(
      gsea_results = data.frame(pathway_id = "test"),
      abundance = "not a data frame",
      metadata = data.frame(),
      group = "group"
    ),
    "must be a data frame or matrix"
  )
})

test_that("pathway_ridgeplot errors on missing group column", {
  skip_if_not_installed("ggridges")

  gsea_results <- data.frame(
    pathway_id = "test",
    pvalue = 0.01,
    direction = "Up"
  )

  abundance <- matrix(1:10, nrow = 5)
  rownames(abundance) <- paste0("K", 1:5)
  colnames(abundance) <- c("S1", "S2")

  metadata <- data.frame(
    sample = c("S1", "S2"),
    other_col = c("A", "B"),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample

  expect_error(
    pathway_ridgeplot(
      gsea_results = gsea_results,
      abundance = abundance,
      metadata = metadata,
      group = "nonexistent"
    ),
    "not found in metadata"
  )
})
