# Helper for compact heatmap test data
create_heatmap_test_data <- function() {
  abundance <- matrix(
    c(
      10, 20, 30, 40,
      5, 15, 25, 35,
      8, 18, 28, 38
    ),
    nrow = 3,
    byrow = TRUE
  )
  rownames(abundance) <- c("Pathway1", "Pathway2", "Pathway3")
  colnames(abundance) <- c("S1", "S2", "S3", "S4")

  metadata <- data.frame(
    sample = c("S1", "S2", "S3", "S4"),
    group = c("A", "A", "B", "B"),
    batch = c("X", "Y", "X", "Y"),
    stringsAsFactors = FALSE
  )

  list(abundance = abundance, metadata = metadata)
}

test_that("pathway_heatmap basic functionality works", {
  td <- create_heatmap_test_data()
  p <- pathway_heatmap(
    abundance = td$abundance,
    metadata = td$metadata,
    group = "group"
  )
  expect_s3_class(p, "ggplot")
})

test_that("pathway_heatmap supports secondary_groups", {
  td <- create_heatmap_test_data()
  p <- pathway_heatmap(
    abundance = td$abundance,
    metadata = td$metadata,
    group = "group",
    secondary_groups = "batch"
  )
  expect_s3_class(p, "ggplot")
})

test_that("pathway_heatmap warns for deprecated facet_by", {
  td <- create_heatmap_test_data()
  expect_warning(
    pathway_heatmap(
      abundance = td$abundance,
      metadata = td$metadata,
      group = "group",
      facet_by = "batch"
    ),
    "deprecated"
  )
})

test_that("pathway_heatmap fails fast when sample IDs do not match", {
  td <- create_heatmap_test_data()
  bad_metadata <- td$metadata
  bad_metadata$sample <- paste0("X", seq_len(nrow(bad_metadata)))

  expect_error(
    pathway_heatmap(
      abundance = td$abundance,
      metadata = bad_metadata,
      group = "group"
    ),
    "Cannot find matching sample identifiers"
  )
})
