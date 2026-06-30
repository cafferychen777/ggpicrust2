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

test_that("pathway_ridgeplot aligns metadata to abundance columns before fold-change calculation", {
  skip_if_not_installed("ggridges")

  abundance <- matrix(
    c(1, 2, 8, 10,
      4, 4, 4, 4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("K00001", "K00002"), c("S1", "S2", "S3", "S4"))
  )
  metadata <- data.frame(
    sample = c("S3", "S4", "S1", "S2"),
    group = c("Treatment", "Treatment", "Control", "Control"),
    stringsAsFactors = FALSE
  )
  pathway_reference <- data.frame(
    pathway_id = "ko00010",
    pathway_name = "Example pathway",
    ko_members = "K00001;K00002",
    stringsAsFactors = FALSE
  )
  gsea_results <- data.frame(
    pathway_id = "ko00010",
    pathway_name = "Example pathway",
    NES = 1.5,
    pvalue = 0.01,
    p.adjust = 0.02,
    stringsAsFactors = FALSE
  )

  plot <- pathway_ridgeplot(
    gsea_results = gsea_results,
    abundance = abundance,
    metadata = metadata,
    group = "group",
    pathway_reference = pathway_reference
  )

  k00001_fc <- plot$data$log2fc[plot$data$gene_id == "K00001"]
  expect_gt(k00001_fc, 0)
})

test_that("pathway_ridgeplot uses factor levels rather than sample order for default two-group comparison", {
  skip_if_not_installed("ggridges")

  abundance <- matrix(
    c(8, 10, 1, 2,
      4, 4, 4, 4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("K00001", "K00002"), c("S3", "S4", "S1", "S2"))
  )
  metadata <- data.frame(
    sample = c("S1", "S2", "S3", "S4"),
    group = factor(c("Control", "Control", "Treatment", "Treatment"),
                   levels = c("Control", "Treatment")),
    stringsAsFactors = FALSE
  )
  pathway_reference <- data.frame(
    pathway_id = "ko00010",
    pathway_name = "Example pathway",
    ko_members = "K00001;K00002",
    stringsAsFactors = FALSE
  )
  gsea_results <- data.frame(
    pathway_id = "ko00010",
    pathway_name = "Example pathway",
    NES = 1.5,
    pvalue = 0.01,
    p.adjust = 0.02,
    stringsAsFactors = FALSE
  )

  plot <- pathway_ridgeplot(
    gsea_results = gsea_results,
    abundance = abundance,
    metadata = metadata,
    group = "group",
    pathway_reference = pathway_reference
  )

  k00001_fc <- plot$data$log2fc[plot$data$gene_id == "K00001"]
  expect_gt(k00001_fc, 0)
})

test_that("pathway_ridgeplot requires explicit comparison for multi-group metadata", {
  skip_if_not_installed("ggridges")

  abundance <- matrix(
    c(1, 1, 8, 8, 20, 20),
    nrow = 1,
    dimnames = list("K00001", paste0("S", 1:6))
  )
  metadata <- data.frame(
    sample = paste0("S", 1:6),
    group = factor(c("A", "A", "B", "B", "C", "C"),
                   levels = c("A", "B", "C")),
    stringsAsFactors = FALSE
  )
  pathway_reference <- data.frame(
    pathway_id = "ko00010",
    ko_members = "K00001",
    stringsAsFactors = FALSE
  )
  gsea_results <- data.frame(
    pathway_id = "ko00010",
    NES = 1,
    pvalue = 0.01,
    p.adjust = 0.02,
    stringsAsFactors = FALSE
  )

  expect_error(
    pathway_ridgeplot(
      gsea_results = gsea_results,
      abundance = abundance,
      metadata = metadata,
      group = "group",
      pathway_reference = pathway_reference
    ),
    "requires 'comparison = c\\(group1, group2\\)'"
  )

  plot <- pathway_ridgeplot(
    gsea_results = gsea_results,
    abundance = abundance,
    metadata = metadata,
    group = "group",
    comparison = c("A", "C"),
    pathway_reference = pathway_reference
  )

  expect_gt(plot$data$log2fc[plot$data$gene_id == "K00001"], 0)
})

test_that("pathway_ridgeplot validates GSEA statistics and display parameters", {
  skip_if_not_installed("ggridges")

  abundance <- matrix(
    c(1, 2, 8, 10),
    nrow = 1,
    dimnames = list("K00001", c("S1", "S2", "S3", "S4"))
  )
  metadata <- data.frame(
    sample = c("S1", "S2", "S3", "S4"),
    group = c("Control", "Control", "Treatment", "Treatment"),
    stringsAsFactors = FALSE
  )
  pathway_reference <- data.frame(
    pathway_id = "ko00010",
    ko_members = "K00001",
    stringsAsFactors = FALSE
  )
  gsea_results <- data.frame(
    pathway_id = "ko00010",
    NES = 1,
    pvalue = 0.01,
    p.adjust = 2,
    stringsAsFactors = FALSE
  )

  expect_error(
    pathway_ridgeplot(
      gsea_results = gsea_results,
      abundance = abundance,
      metadata = metadata,
      group = "group",
      pathway_reference = pathway_reference
    ),
    "between 0 and 1"
  )

  gsea_results$p.adjust <- 0.01
  expect_error(
    pathway_ridgeplot(
      gsea_results = gsea_results,
      abundance = abundance,
      metadata = metadata,
      group = "group",
      pathway_reference = pathway_reference,
      sort_by = "invalid"
    ),
    "'sort_by' must be one of"
  )

  expect_error(
    pathway_ridgeplot(
      gsea_results = gsea_results,
      abundance = abundance,
      metadata = metadata,
      group = "group",
      pathway_reference = pathway_reference,
      n_pathways = 0
    ),
    "n_pathways"
  )
})

test_that("pathway_ridgeplot validates pathway identifiers and label fallback", {
  skip_if_not_installed("ggridges")

  abundance <- matrix(
    c(1, 2, 8, 10),
    nrow = 1,
    dimnames = list("K00001", c("S1", "S2", "S3", "S4"))
  )
  metadata <- data.frame(
    sample = c("S1", "S2", "S3", "S4"),
    group = c("Control", "Control", "Treatment", "Treatment"),
    stringsAsFactors = FALSE
  )
  pathway_reference <- data.frame(
    pathway_id = "ko00010",
    ko_members = "K00001",
    stringsAsFactors = FALSE
  )
  gsea_results <- data.frame(
    pathway_id = "ko00010",
    pathway_name = NA_character_,
    NES = 1,
    pvalue = 0.01,
    p.adjust = 0.02,
    stringsAsFactors = FALSE
  )

  plot <- pathway_ridgeplot(
    gsea_results = gsea_results,
    abundance = abundance,
    metadata = metadata,
    group = "group",
    pathway_reference = pathway_reference
  )
  expect_identical(as.character(unique(plot$data$pathway)), "ko00010")

  gsea_bad_id <- gsea_results
  gsea_bad_id$pathway_id[1] <- NA_character_
  expect_error(
    pathway_ridgeplot(
      gsea_results = gsea_bad_id,
      abundance = abundance,
      metadata = metadata,
      group = "group",
      pathway_reference = pathway_reference
    ),
    "pathway_id.*non-empty"
  )

  gsea_duplicate <- rbind(
    transform(gsea_results, pathway_name = "First"),
    transform(gsea_results, pathway_name = "Second")
  )
  expect_error(
    pathway_ridgeplot(
      gsea_results = gsea_duplicate,
      abundance = abundance,
      metadata = metadata,
      group = "group",
      pathway_reference = pathway_reference
    ),
    "duplicate pathway_id"
  )

  expect_error(
    pathway_ridgeplot(
      gsea_results = gsea_results,
      abundance = abundance,
      metadata = metadata,
      group = "group",
      pathway_reference = pathway_reference,
      pathway_type = NA_character_
    ),
    "pathway_type must be one of"
  )
})

test_that("pathway_ridgeplot sorts NES by absolute effect size", {
  skip_if_not_installed("ggridges")

  abundance <- matrix(
    c(1, 2, 8, 10,
      1, 1, 5, 5,
      5, 5, 1, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("K00001", "K00002", "K00003"), c("S1", "S2", "S3", "S4"))
  )
  metadata <- data.frame(
    sample = c("S1", "S2", "S3", "S4"),
    group = c("Control", "Control", "Treatment", "Treatment"),
    stringsAsFactors = FALSE
  )
  pathway_reference <- data.frame(
    pathway_id = c("ko00010", "ko00020", "ko00030"),
    pathway_name = c("Small negative", "Largest positive", "Middle negative"),
    ko_members = c("K00001", "K00002", "K00003"),
    stringsAsFactors = FALSE
  )
  gsea_results <- data.frame(
    pathway_id = c("ko00010", "ko00020", "ko00030"),
    pathway_name = c("Small negative", "Largest positive", "Middle negative"),
    NES = c(-1, 3, -2),
    pvalue = c(0.03, 0.02, 0.01),
    p.adjust = c(0.03, 0.02, 0.01),
    stringsAsFactors = FALSE
  )

  plot <- pathway_ridgeplot(
    gsea_results = gsea_results,
    abundance = abundance,
    metadata = metadata,
    group = "group",
    pathway_reference = pathway_reference,
    n_pathways = 1,
    sort_by = "NES"
  )

  expect_identical(as.character(unique(plot$data$pathway)), "Largest positive")
})
