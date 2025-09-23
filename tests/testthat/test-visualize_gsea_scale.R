# Tests for `scale` parameter in visualize_gsea

library(testthat)

# Load the current development version of visualize_gsea into this session
suppressWarnings(suppressMessages(source(file.path("..","..","R","visualize_gsea.R"))))

# Minimal helpers local to this file
create_test_gsea_results_scale <- function(n_pathways = 10) {
  set.seed(42)
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
    pathway_class = sample(c("A","B","C"), n_pathways, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

create_test_abundance_scale <- function(n_genes = 100, n_samples = 10) {
  set.seed(42)
  gene_ids <- paste0("K", sprintf("%05d", 1:n_genes))
  sample_ids <- paste0("Sample", 1:n_samples)
  mat <- matrix(runif(n_genes * n_samples, 1, 100), nrow = n_genes, ncol = n_samples,
                dimnames = list(gene_ids, sample_ids))
  as.data.frame(mat)
}

create_test_metadata_scale <- function(n_samples = 10) {
  set.seed(42)
  sample_ids <- paste0("Sample", 1:n_samples)
  groups <- rep(c("Group1","Group2"), each = n_samples/2)
  md <- data.frame(sample = sample_ids, group = groups, stringsAsFactors = FALSE)
  rownames(md) <- md$sample
  md
}

# Utility to get scale classes from ggplot object
.get_scale_classes <- function(p) {
  if (is.null(p$scales) || length(p$scales$scales) == 0) return(character())
  unique(unlist(lapply(p$scales$scales, class)))
}


test_that("enrichment_plot honors scale palette/vector/function", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("enrichplot")

  gsea_results <- create_test_gsea_results_scale()

  # vector
  pal_vec <- grDevices::colorRampPalette(c("black","white","yellow"))(5)
  p_vec <- visualize_gsea(gsea_results, plot_type = "enrichment_plot", scale = pal_vec)
  expect_s3_class(p_vec, "ggplot")

  # function
  pal_fun <- grDevices::terrain.colors
  p_fun <- visualize_gsea(gsea_results, plot_type = "enrichment_plot", scale = pal_fun)
  expect_s3_class(p_fun, "ggplot")
})


test_that("dotplot honors scale palette/function", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("enrichplot")

  gsea_results <- create_test_gsea_results_scale()

  pal_fun <- grDevices::heat.colors
  p <- visualize_gsea(gsea_results, plot_type = "dotplot", scale = pal_fun)
  expect_s3_class(p, "ggplot")
})


test_that("barplot honors discrete two-color palette for direction", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("enrichplot")

  gsea_results <- create_test_gsea_results_scale()

  pal_vec <- c("#111111", "#FFEE00")
  p <- visualize_gsea(gsea_results, plot_type = "barplot", scale = pal_vec)
  expect_s3_class(p, "ggplot")
  # manual discrete scale
  classes <- .get_scale_classes(p)
  expect_true(any(grepl("Scale.*Discrete", classes)))
})


test_that("network accepts custom scale (diverging or sequential)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")

  gsea_results <- create_test_gsea_results_scale()

  pal_vec <- c("#440154", "#21908C", "#FDE725")
  p <- visualize_gsea(gsea_results, plot_type = "network", scale = pal_vec)
  expect_s3_class(p, "ggplot")
})


test_that("heatmap accepts scale or explicit col_fun via heatmap_params", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  gsea_results <- create_test_gsea_results_scale()
  abundance <- create_test_abundance_scale()
  metadata <- create_test_metadata_scale()

  # via scale
  pal_vec <- grDevices::colorRampPalette(c("navy","white","firebrick"))(5)
  p1 <- visualize_gsea(gsea_results, plot_type = "heatmap",
                       abundance = abundance, metadata = metadata, group = "group",
                       scale = pal_vec)
  expect_s4_class(p1, "Heatmap")

  # via explicit col_fun (takes precedence). Construct col_fun inline to avoid relying on exported helper.
  skip_if_not_installed("circlize")
  col_fun <- circlize::colorRamp2(c(-3, 0, 3), c("#2166AC", "#FFFFFF", "#B2182B"))
  p2 <- visualize_gsea(gsea_results, plot_type = "heatmap",
                       abundance = abundance, metadata = metadata, group = "group",
                       heatmap_params = list(col_fun = col_fun))
  expect_s4_class(p2, "Heatmap")
})

