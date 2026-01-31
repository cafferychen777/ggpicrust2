# Tests for pathway_gsea function

# Helper: create standard GSEA test data
create_gsea_test_data <- function(n_features = 50, n_samples = 10) {
  set.seed(123)
  abundance <- matrix(rnorm(n_features * n_samples), nrow = n_features, ncol = n_samples)
  colnames(abundance) <- paste0("Sample", 1:n_samples)
  rownames(abundance) <- paste0("K", sprintf("%05d", 1:n_features))

  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("Control", "Treatment"), each = n_samples / 2)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample_name

  list(abundance = abundance, metadata = metadata)
}

test_that("pathway_gsea works with fgsea method", {
  skip_if_not_installed("fgsea")

  test_data <- create_gsea_test_data()

  local_mocked_bindings(
    prepare_gene_sets = function(...) {
      list("path:ko00001" = c("K00001", "K00002"), "path:ko00002" = c("K00002", "K00003"))
    },
    run_fgsea = function(...) {
      data.frame(
        pathway_id = c("path:ko00001", "path:ko00002"),
        pathway_name = c("path:ko00001", "path:ko00002"),
        size = c(2, 2), ES = c(0.5, -0.3), NES = c(1.2, -0.8),
        pvalue = c(0.01, 0.05), p.adjust = c(0.02, 0.1),
        leading_edge = c("K00001;K00002", "K00003"),
        stringsAsFactors = FALSE
      )
    }
  )

  result <- pathway_gsea(
    abundance = test_data$abundance,
    metadata = test_data$metadata,
    group = "group",
    method = "fgsea"
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("pathway_id", "NES", "pvalue", "method") %in% colnames(result)))
})

test_that("pathway_gsea works with camera method", {
  skip_if_not_installed("limma")

  test_data <- create_gsea_test_data()
  abundance <- abs(test_data$abundance) * 100

  local_mocked_bindings(
    prepare_gene_sets = function(...) {
      list(
        "ko00010" = paste0("K", sprintf("%05d", 1:10)),
        "ko00020" = paste0("K", sprintf("%05d", 11:20))
      )
    }
  )

  result <- pathway_gsea(
    abundance = abundance,
    metadata = test_data$metadata,
    group = "group",
    method = "camera",
    pathway_type = "KEGG"
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("pathway_id", "direction", "pvalue") %in% colnames(result)))
  expect_equal(unique(result$method), "camera")
})

test_that("pathway_gsea validates inputs correctly", {
  test_data <- create_gsea_test_data()

  expect_error(pathway_gsea(abundance = "invalid", metadata = test_data$metadata, group = "group"),
               "'abundance' must be a data frame or matrix")
  expect_error(pathway_gsea(abundance = test_data$abundance, metadata = "invalid", group = "group"),
               "'metadata' must be a data frame")
  expect_error(pathway_gsea(abundance = test_data$abundance, metadata = test_data$metadata, group = "invalid_group"),
               "not found in metadata")
  expect_error(pathway_gsea(abundance = test_data$abundance, metadata = test_data$metadata, group = "group", method = "invalid"),
               "method must be one of")
})

test_that("prepare_gene_sets works for KEGG and MetaCyc pathway types", {
  gene_sets_kegg <- prepare_gene_sets("KEGG")
  expect_type(gene_sets_kegg, "list")
  expect_true(length(gene_sets_kegg) > 0)

  gene_sets_metacyc <- prepare_gene_sets("MetaCyc")
  expect_type(gene_sets_metacyc, "list")
})
