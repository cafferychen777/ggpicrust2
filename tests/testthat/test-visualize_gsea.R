# Test for visualize_gsea function
# create_gsea_test_results() is defined in helper-gsea.R

test_that("visualize_gsea creates different plot types correctly", {
  # enrichment_plot/dotplot/barplot are pure ggplot2 and do not require
  # enrichplot; asserting that here locks the dependency boundary.
  gsea_results <- create_gsea_test_results()

  p_bar <- visualize_gsea(gsea_results, plot_type = "barplot")
  expect_s3_class(p_bar, "ggplot")

  p_dot <- visualize_gsea(gsea_results, plot_type = "dotplot")
  expect_s3_class(p_dot, "ggplot")

  p_enrich <- visualize_gsea(gsea_results, plot_type = "enrichment_plot")
  expect_s3_class(p_enrich, "ggplot")
})

test_that("visualize_gsea labels signed p-value scores from camera/fry distinctly from NES", {
  gsea_results <- create_gsea_test_results()
  gsea_results$signed_log10_pvalue <- gsea_results$NES
  gsea_results$score_type <- "signed_log10_pvalue"
  gsea_results$score_label <- "Signed -log10(p-value)"

  p_bar <- visualize_gsea(gsea_results, plot_type = "barplot")
  expect_equal(p_bar$labels$y, "Signed -log10(p-value)")

  p_dot <- visualize_gsea(gsea_results, plot_type = "dotplot")
  expect_equal(p_dot$labels$x, "Signed -log10(p-value)")

  mixed_scores <- gsea_results
  mixed_scores$score_label[1] <- "Normalized Enrichment Score (NES)"
  expect_error(
    visualize_gsea(mixed_scores, plot_type = "barplot"),
    "multiple score labels"
  )
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
  gsea_results$leading_edge <- vapply(seq_len(nrow(gsea_results)), function(i) {
    paste(paste0("K", sprintf("%05d", ((i - 1) * 5 + 1):(i * 5))), collapse = ";")
  }, character(1))

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

test_that("visualize_gsea heatmap handles PICRUSt2 #NAME abundance input", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  gsea_results <- data.frame(
    pathway_id = "path:ko00010",
    pathway_name = "Test pathway",
    NES = 1.5,
    p.adjust = 0.01,
    leading_edge = "K00001;K00002",
    stringsAsFactors = FALSE
  )
  abundance <- data.frame(
    `#NAME` = c("K00001", "K00002", "K00003"),
    S1 = c(1, 2, 10),
    S2 = c(2, 3, 10),
    S3 = c(8, 9, 10),
    S4 = c(9, 10, 10),
    check.names = FALSE
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  p <- visualize_gsea(
    gsea_results = gsea_results,
    plot_type = "heatmap",
    abundance = abundance,
    metadata = metadata,
    group = "group"
  )

  expect_s4_class(p, "Heatmap")
  expect_true(any(p@matrix != 0))
})

test_that("visualize_gsea heatmap handles generic leading feature ID abundance input", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  gsea_results <- data.frame(
    pathway_id = "path:ko00010",
    pathway_name = "Test pathway",
    NES = 1.5,
    p.adjust = 0.01,
    leading_edge = "K00001;K00002",
    stringsAsFactors = FALSE
  )
  abundance <- data.frame(
    feature = c("K00001", "K00002", "K00003"),
    S1 = c(1, 2, 10),
    S2 = c(2, 3, 10),
    S3 = c(8, 9, 10),
    S4 = c(9, 10, 10),
    check.names = FALSE
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  p <- visualize_gsea(
    gsea_results = gsea_results,
    plot_type = "heatmap",
    abundance = abundance,
    metadata = metadata,
    group = "group"
  )

  expect_s4_class(p, "Heatmap")
  expect_true(any(p@matrix != 0))
})

test_that("visualize_gsea heatmap rejects unmatched leading-edge genes", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  gsea_results <- data.frame(
    pathway_id = "path:ko00010",
    pathway_name = "Test pathway",
    NES = 1.5,
    p.adjust = 0.01,
    leading_edge = "K99999",
    stringsAsFactors = FALSE
  )
  abundance <- matrix(
    seq_len(3 * 4),
    nrow = 3,
    dimnames = list(paste0("K", sprintf("%05d", 1:3)), paste0("S", 1:4))
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    visualize_gsea(
      gsea_results = gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group"
    ),
    "None of the non-empty leading-edge genes"
  )
})

test_that("visualize_gsea heatmap rejects missing group annotations after alignment", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  gsea_results <- data.frame(
    pathway_id = "path:ko00010",
    pathway_name = "Test pathway",
    NES = 1.5,
    p.adjust = 0.01,
    leading_edge = "K00001",
    stringsAsFactors = FALSE
  )
  abundance <- matrix(
    seq_len(3 * 4),
    nrow = 3,
    dimnames = list(paste0("K", sprintf("%05d", 1:3)), paste0("S", 1:4))
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", NA),
    stringsAsFactors = FALSE
  )

  expect_error(
    visualize_gsea(
      gsea_results = gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group"
    ),
    "contains NA values after sample alignment"
  )
})

test_that("visualize_gsea validates inputs correctly", {
  gsea_results <- create_gsea_test_results()

  expect_error(visualize_gsea(gsea_results = "invalid"), "'gsea_results' must be a data frame")
  expect_error(visualize_gsea(gsea_results, plot_type = "invalid"), "'plot_type' must be one of")
  expect_error(visualize_gsea(gsea_results, sort_by = "invalid"), "'sort_by' must be one of")
  expect_error(
    visualize_gsea(gsea_results, pathway_label_column = c("pathway_id", "pathway_name")),
    "single non-empty character string"
  )
  expect_error(
    visualize_gsea(gsea_results, pathway_label_column = ""),
    "single non-empty character string"
  )
})

test_that("visualize_gsea validates selected pathway identifiers and labels", {
  gsea_results <- create_gsea_test_results(n_pathways = 3)

  bad_id <- gsea_results
  bad_id$pathway_id[1] <- NA_character_
  expect_error(
    visualize_gsea(bad_id, plot_type = "barplot"),
    "pathway_id.*non-empty"
  )

  duplicate_id <- gsea_results
  duplicate_id$pathway_id[2] <- duplicate_id$pathway_id[1]
  expect_error(
    visualize_gsea(duplicate_id, plot_type = "barplot", n_pathways = 3),
    "one row per pathway_id"
  )

  duplicate_outside_selection <- gsea_results
  duplicate_outside_selection$pathway_id[3] <- duplicate_outside_selection$pathway_id[1]
  duplicate_outside_selection$p.adjust <- c(0.01, 0.02, 0.99)
  p_selected <- visualize_gsea(
    duplicate_outside_selection,
    plot_type = "barplot",
    n_pathways = 2
  )
  expect_s3_class(p_selected, "ggplot")

  missing_labels <- gsea_results
  missing_labels$pathway_name <- c(NA_character_, "", "Named pathway")
  p_labeled <- visualize_gsea(missing_labels, plot_type = "barplot", n_pathways = 3)
  labels <- as.character(p_labeled$data$pathway_label)

  expect_false(any(is.na(labels)))
  expect_true(all(c(missing_labels$pathway_id[1],
                    missing_labels$pathway_id[2],
                    "Named pathway") %in% labels))
})

test_that("visualize_gsea validates GSEA statistic columns before plotting", {
  gsea_results <- create_gsea_test_results()

  bad_p <- gsea_results
  bad_p$p.adjust[1] <- 1.2
  expect_error(
    visualize_gsea(bad_p, plot_type = "enrichment_plot"),
    "between 0 and 1"
  )

  bad_nes <- gsea_results
  bad_nes$NES[1] <- Inf
  expect_error(
    visualize_gsea(bad_nes, plot_type = "barplot"),
    "finite numeric"
  )

  bad_size <- gsea_results
  bad_size$size[1] <- 1.5
  expect_error(
    visualize_gsea(bad_size, plot_type = "dotplot"),
    "positive integer"
  )

  huge_size <- gsea_results
  huge_size$size[1] <- 1e20
  expect_error(
    expect_warning(visualize_gsea(huge_size, plot_type = "dotplot"), NA),
    "positive integer"
  )

  expect_error(
    visualize_gsea(gsea_results, n_pathways = 0),
    "n_pathways"
  )
  expect_error(
    expect_warning(visualize_gsea(gsea_results, n_pathways = 1e20), NA),
    "n_pathways.*single finite integer"
  )
})

test_that("visualize_gsea limits pathways correctly", {
  gsea_results <- create_gsea_test_results(n_pathways = 30)

  p_default <- visualize_gsea(gsea_results, plot_type = "barplot")
  expect_equal(nrow(p_default$data), 20)

  p_custom <- visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 10)
  expect_equal(nrow(p_custom$data), 10)
})

test_that("visualize_gsea network treats missing leading edges as empty sets", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")

  gsea_results <- data.frame(
    pathway_id = c("p1", "p2"),
    pathway_name = c("Pathway 1", "Pathway 2"),
    size = c(10, 20),
    NES = c(1.5, -1.2),
    pvalue = c(0.01, 0.02),
    p.adjust = c(0.03, 0.04),
    leading_edge = c(NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )

  p <- visualize_gsea(
    gsea_results,
    plot_type = "network",
    network_params = list(similarity_cutoff = 0)
  )

  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "No pathways to display for network")
})

test_that("visualize_gsea network rejects duplicate pathway IDs", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")

  gsea_results <- data.frame(
    pathway_id = c("p1", "p1", "p2"),
    pathway_name = c("Pathway 1a", "Pathway 1b", "Pathway 2"),
    size = c(10, 12, 20),
    NES = c(1.5, 1.2, -1.2),
    pvalue = c(0.01, 0.02, 0.03),
    p.adjust = c(0.03, 0.04, 0.05),
    leading_edge = c("K00001;K00002", "K00002;K00003", "K00002;K00004"),
    stringsAsFactors = FALSE
  )

  expect_error(
    visualize_gsea(
      gsea_results,
      plot_type = "network",
      network_params = list(similarity_cutoff = 0)
    ),
    "Duplicate pathway_id values"
  )
})
