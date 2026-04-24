create_followup_errorbar_test_data <- function(n_features = 5, p_adjust = NULL) {
  set.seed(123)
  n_samples <- 10
  abundance <- matrix(runif(n_features * n_samples), nrow = n_features, ncol = n_samples)
  rownames(abundance) <- paste0("pathway", seq_len(n_features))
  colnames(abundance) <- paste0("sample", seq_len(n_samples))

  if (is.null(p_adjust)) p_adjust <- rep(0.01, n_features)

  daa_results_df <- data.frame(
    feature = paste0("pathway", seq_len(n_features)),
    pathway_name = paste0("Pathway ", seq_len(n_features)),
    p_adjust = p_adjust,
    method = rep("ALDEx2_Welch's t test", n_features),
    group1 = rep("GroupA", n_features),
    group2 = rep("GroupB", n_features),
    stringsAsFactors = FALSE
  )

  group <- factor(rep(c("GroupA", "GroupB"), each = n_samples / 2))
  names(group) <- colnames(abundance)

  list(abundance = abundance, daa_results_df = daa_results_df, Group = group)
}

test_that("import_MicrobiomeAnalyst_daa_results defaults and validates shape", {
  daa <- data.frame(
    feature = c("p1", "p2"),
    p_values = c(0.01, 0.02),
    p_adjust = c(0.03, 0.04),
    Statistics = c(1.2, -0.5),
    stringsAsFactors = FALSE
  )

  imported <- import_MicrobiomeAnalyst_daa_results(data = daa)
  expect_equal(imported$group1, c("control", "control"))
  expect_equal(imported$group2, c("treatment", "treatment"))
  expect_false(anyNA(names(imported)))

  too_short <- daa[, 1:3]
  expect_error(
    import_MicrobiomeAnalyst_daa_results(data = too_short),
    "at least four columns"
  )

  extra_col <- cbind(daa, extra = 1)
  imported_extra <- import_MicrobiomeAnalyst_daa_results(data = extra_col)
  expect_false(anyNA(names(imported_extra)))
  expect_true("extra" %in% names(imported_extra))
})

test_that("import_MicrobiomeAnalyst_daa_results warns that data wins over file_path", {
  daa <- data.frame(
    feature = "p1",
    p_values = 0.01,
    p_adjust = 0.03,
    Statistics = 1,
    stringsAsFactors = FALSE
  )

  expect_warning(
    imported <- import_MicrobiomeAnalyst_daa_results(
      file_path = tempfile(fileext = ".csv"),
      data = daa
    ),
    "Using data and ignoring file_path"
  )
  expect_equal(imported$feature, "p1")
})

test_that("pathway_errorbar order='group' resolves tied groups deterministically", {
  abundance <- matrix(
    c(10, 20, 30, 40,
      20, 10, 40, 30),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("f1", "f2"), paste0("S", 1:4))
  )
  daa <- data.frame(
    feature = c("f1", "f2"),
    method = "m",
    group1 = "A",
    group2 = "B",
    p_adjust = c(0.01, 0.02),
    description = c("d1", "d2"),
    stringsAsFactors = FALSE
  )
  group <- c("A", "A", "B", "B")
  names(group) <- colnames(abundance)

  expect_warning(
    plot <- pathway_errorbar(
      abundance,
      daa,
      group,
      order = "group",
      x_lab = "description",
      p_value_bar = "FALSE"
    ),
    NA
  )
  expect_s3_class(plot, "patchwork")
})

test_that("pathway_errorbar normalizes logical flags and rejects invalid p_adjust", {
  td <- create_followup_errorbar_test_data(n_features = 2, p_adjust = c(0, 0.01))
  expect_error(
    plot <- pathway_errorbar(
      td$abundance,
      td$daa_results_df,
      td$Group,
      x_lab = "pathway_name",
      p_value_bar = "FALSE",
      smart_colors = "FALSE",
      accessibility_mode = "FALSE"
    ),
    NA
  )
  expect_s3_class(plot, "patchwork")

  td$daa_results_df$p_adjust[1] <- -0.1
  expect_error(
    pathway_errorbar(td$abundance, td$daa_results_df, td$Group, x_lab = "pathway_name"),
    "non-negative"
  )
})

test_that("pathway_errorbar_table gives ko_to_kegg an explicit contract", {
  td <- create_followup_errorbar_test_data(n_features = 2, p_adjust = c(0.01, 0.02))

  expect_error(
    pathway_errorbar_table(
      abundance = td$abundance,
      daa_results_df = td$daa_results_df,
      Group = td$Group,
      ko_to_kegg = TRUE
    ),
    "pathway_class"
  )

  td$daa_results_df$pathway_class <- c("Metabolism", "Metabolism")
  expect_error(
    res <- pathway_errorbar_table(
      abundance = td$abundance,
      daa_results_df = td$daa_results_df,
      Group = td$Group,
      ko_to_kegg = "TRUE"
    ),
    NA
  )
  expect_true("pathway_class" %in% names(res))
})

test_that("gsea_pathway_annotation preserves input order", {
  kegg_ref <- ggpicrust2:::load_reference_data("KEGG")
  ids <- rev(head(kegg_ref$pathway, 3))
  gsea_results <- data.frame(
    pathway_id = ids,
    NES = c(1, 2, 3),
    pvalue = c(0.03, 0.02, 0.01),
    p.adjust = c(0.03, 0.02, 0.01),
    stringsAsFactors = FALSE
  )

  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")
  expect_equal(annotated$pathway_id, ids)
})

test_that("visualize_gsea validates network contract before delayed plot build", {
  gsea_results <- create_gsea_test_results()

  missing_leading_edge <- gsea_results[, setdiff(names(gsea_results), "leading_edge")]
  expect_error(
    visualize_gsea(missing_leading_edge, plot_type = "network"),
    "leading_edge"
  )

  expect_error(
    visualize_gsea(
      gsea_results,
      plot_type = "network",
      network_params = list(similarity_measure = "bad")
    ),
    "network_params\\$similarity_measure"
  )

  expect_error(
    visualize_gsea(
      gsea_results,
      plot_type = "network",
      network_params = list(node_color_by = "missing")
    ),
    "network_params\\$node_color_by"
  )
})

test_that("visualize_gsea honors edge_width_by", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")

  gsea_results <- create_gsea_test_results(n_pathways = 3)
  gsea_results$leading_edge <- c("K1;K2", "K2;K3", "K3;K4")

  similarity_plot <- visualize_gsea(
    gsea_results,
    plot_type = "network",
    network_params = list(edge_width_by = "similarity", similarity_cutoff = 0)
  )
  constant_plot <- visualize_gsea(
    gsea_results,
    plot_type = "network",
    network_params = list(edge_width_by = "constant", similarity_cutoff = 0)
  )

  expect_true("edge_width" %in% names(similarity_plot$layers[[1]]$mapping))
  expect_false("edge_width" %in% names(constant_plot$layers[[1]]$mapping))
})

test_that("visualize_gsea heatmap treats annotation_colors=NULL as default", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  gsea_results <- create_gsea_test_results(n_pathways = 2)
  gsea_results$pathway_id <- c("p1", "p2")
  gsea_results$leading_edge <- c("K00001;K00002", "K00002;K00003")
  abundance <- matrix(
    1:12,
    nrow = 3,
    dimnames = list(c("K00001", "K00002", "K00003"), paste0("S", 1:4))
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    heatmap <- visualize_gsea(
      gsea_results,
      plot_type = "heatmap",
      abundance = abundance,
      metadata = metadata,
      group = "group",
      heatmap_params = list(annotation_colors = NULL)
    ),
    NA
  )
  expect_s4_class(heatmap, "Heatmap")
})

test_that("compare_metagenome_results checks heatmap namespaces explicitly", {
  body_src <- paste(deparse(body(ggpicrust2::compare_metagenome_results)), collapse = "\n")
  expect_true(grepl("require_package\\(\"circlize\"", body_src))
  expect_true(grepl("require_package\\(\"ComplexHeatmap\"", body_src))
})
