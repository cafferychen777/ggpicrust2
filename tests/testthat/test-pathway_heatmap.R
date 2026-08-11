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

test_that("pathway_heatmap preserves a leading feature ID column", {
  abundance <- data.frame(
    feature = c("PathA", "PathB", "PathC"),
    S1 = c(10, 5, 8),
    S2 = c(20, 15, 18),
    S3 = c(30, 25, 28),
    S4 = c(40, 35, 38),
    check.names = FALSE
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  p <- pathway_heatmap(
    abundance = abundance,
    metadata = metadata,
    group = "group"
  )

  expect_s3_class(p, "ggplot")
  expect_setequal(as.character(unique(p$data$rowname)), abundance$feature)
  expect_false(any(grepl("^Pathway[0-9]+$", as.character(unique(p$data$rowname)))))
})

test_that("pathway_heatmap accepts finite transformed zero-sum sample columns", {
  abundance <- matrix(
    c(
      -1, 1, 2, 4,
       0, 2, 3, 5,
       1, 3, 4, 6
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Pathway", 1:3), paste0("S", 1:4))
  )
  metadata <- data.frame(
    sample = colnames(abundance),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  p <- pathway_heatmap(abundance, metadata, "group")
  expect_s3_class(p, "ggplot")
  expect_setequal(as.character(unique(p$data$Sample)), colnames(abundance))
})

test_that("pathway_heatmap rejects missing and non-finite abundance values", {
  td <- create_heatmap_test_data()

  abundance_na <- td$abundance
  abundance_na[1, 1] <- NA_real_
  expect_error(
    pathway_heatmap(abundance_na, td$metadata, "group"),
    "missing values"
  )

  abundance_inf <- td$abundance
  abundance_inf[1, 1] <- Inf
  expect_error(
    pathway_heatmap(abundance_inf, td$metadata, "group"),
    "non-finite values"
  )
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

test_that("pathway_heatmap revalidates grouping variables after sample alignment", {
  abundance <- matrix(
    seq_len(3 * 4),
    nrow = 3,
    dimnames = list(paste0("Pathway", 1:3), paste0("S", 1:4))
  )

  metadata_secondary_collapses <- data.frame(
    sample = paste0("S", 1:5),
    group = c("A", "A", "B", "B", "B"),
    batch = c("X", "X", "X", "X", "Y"),
    stringsAsFactors = FALSE
  )
  expect_error(
    pathway_heatmap(
      abundance = abundance,
      metadata = metadata_secondary_collapses,
      group = "group",
      secondary_groups = "batch"
    ),
    "must contain at least 2 group"
  )

  metadata_group_na <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", NA),
    stringsAsFactors = FALSE
  )
  expect_error(
    pathway_heatmap(
      abundance = abundance,
      metadata = metadata_group_na,
      group = "group"
    ),
    "non-missing, non-empty group labels"
  )

  metadata_group_blank <- metadata_group_na
  metadata_group_blank$group <- c("A", "A", "B", " ")
  expect_error(
    pathway_heatmap(
      abundance = abundance,
      metadata = metadata_group_blank,
      group = "group"
    ),
    "non-missing, non-empty group labels"
  )
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

test_that("pathway_heatmap handles constant rows with correlation clustering", {
  abundance <- matrix(
    c(
      1, 1, 1, 1,
      2, 3, 4, 5,
      5, 4, 3, 2
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("constant", "up", "down"), paste0("S", 1:4))
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  expect_message(
    p <- pathway_heatmap(
      abundance = abundance,
      metadata = metadata,
      group = "group",
      cluster_rows = TRUE,
      clustering_distance = "correlation"
    ),
    "Undefined pearson correlation"
  )
  expect_s3_class(p, "ggplot")

  expect_message(
    p_spearman <- pathway_heatmap(
      abundance = abundance,
      metadata = metadata,
      group = "group",
      cluster_rows = TRUE,
      clustering_distance = "spearman"
    ),
    "Undefined spearman correlation"
  )
  expect_s3_class(p_spearman, "ggplot")
})

test_that("pathway_heatmap validates Ward linkage distance semantics", {
  td <- create_heatmap_test_data()

  expect_error(
    pathway_heatmap(
      abundance = td$abundance,
      metadata = td$metadata,
      group = "group",
      cluster_rows = TRUE,
      clustering_method = "ward.D2",
      clustering_distance = "correlation"
    ),
    "Ward clustering methods.*euclidean"
  )

  expect_message(
    p_correlation <- pathway_heatmap(
      abundance = td$abundance,
      metadata = td$metadata,
      group = "group",
      cluster_rows = TRUE,
      clustering_method = "average",
      clustering_distance = "correlation"
    ),
    "Pathways ordered by hierarchical clustering"
  )
  expect_s3_class(p_correlation, "ggplot")

  expect_message(
    p_ward <- pathway_heatmap(
      abundance = td$abundance,
      metadata = td$metadata,
      group = "group",
      cluster_rows = TRUE,
      clustering_method = "ward.D2",
      clustering_distance = "euclidean"
    ),
    "Pathways ordered by hierarchical clustering"
  )
  expect_s3_class(p_ward, "ggplot")
})

test_that("pathway_heatmap handles zero-variance sample profiles with correlation column clustering", {
  abundance <- matrix(
    c(
      1, 1, 1, 1,
      1, 1, 1, 1,
      1, 1, 1, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Pathway", 1:3), paste0("S", 1:4))
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  expect_message(
    p <- pathway_heatmap(
      abundance = abundance,
      metadata = metadata,
      group = "group",
      cluster_cols = TRUE,
      clustering_distance = "correlation"
    ),
    "Undefined pearson correlation"
  )
  expect_s3_class(p, "ggplot")
})

test_that("pathway_heatmap does not overwrite a group named sample_name", {
  td <- create_heatmap_test_data()
  metadata <- data.frame(
    sample = td$metadata$sample,
    sample_name = td$metadata$group,
    stringsAsFactors = FALSE
  )

  p <- pathway_heatmap(td$abundance, metadata, group = "sample_name")

  expect_s3_class(p, "ggplot")
  expect_setequal(p$data$sample_name, c("A", "B"))
  expect_error(ggplot2::ggplot_build(p), NA)
})

test_that("pathway_heatmap validates logical flags and grouping hierarchy", {
  td <- create_heatmap_test_data()

  expect_error(
    pathway_heatmap(td$abundance, td$metadata, "group", show_legend = NA),
    "show_legend.*TRUE or FALSE"
  )
  expect_error(
    pathway_heatmap(
      td$abundance,
      td$metadata,
      "group",
      secondary_groups = "group"
    ),
    "Grouping columns must be unique"
  )
})

test_that("pathway_heatmap requires explicit sample identity", {
  td <- create_heatmap_test_data()
  unnamed_abundance <- td$abundance
  colnames(unnamed_abundance) <- NULL
  misleading_metadata <- td$metadata
  misleading_metadata$sample <- paste0("Sample", seq_len(nrow(misleading_metadata)))

  expect_error(
    pathway_heatmap(unnamed_abundance, misleading_metadata, "group"),
    "sample identifiers in column names"
  )
})

test_that("pathway_heatmap validates clustering cardinality and display parameters", {
  td <- create_heatmap_test_data()
  one_pathway <- td$abundance[1, , drop = FALSE]

  expect_error(
    pathway_heatmap(
      one_pathway,
      td$metadata,
      "group",
      cluster_rows = TRUE
    ),
    "At least two pathways"
  )
  expect_error(
    pathway_heatmap(td$abundance, td$metadata, "group", font_size = NA_real_),
    "font_size.*finite"
  )
  expect_error(
    pathway_heatmap(td$abundance, td$metadata, "group", low_color = "not-a-color"),
    "invalid R color"
  )
  expect_error(
    pathway_heatmap(
      td$abundance,
      td$metadata,
      "group",
      colorbar_position = "somewhere"
    ),
    "colorbar_position.*must be one of"
  )
})
