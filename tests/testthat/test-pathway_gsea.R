# Helper function to create test data for GSEA tests
create_test_data <- function(n_features = 3, n_samples = 10) {
  set.seed(123)
  abundance <- matrix(rnorm(n_features * n_samples), nrow = n_features, ncol = n_samples)
  colnames(abundance) <- paste0("Sample", 1:n_samples)
  rownames(abundance) <- paste0("K", sprintf("%05d", 1:n_features))

  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = factor(rep(c("Control", "Treatment"), each = n_samples/2))
  )
  rownames(metadata) <- metadata$sample_name

  return(list(abundance = abundance, metadata = metadata))
}

test_that("pathway_gsea handles basic input correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("fgsea")

  # Create test data
  test_data <- create_test_data()
  abundance <- test_data$abundance
  metadata <- test_data$metadata

  # Mock the prepare_gene_sets function to return a simple gene set
  mockery::stub(pathway_gsea, "prepare_gene_sets", function(...) {
    list(
      "path:ko00001" = c("K00001", "K00002"),
      "path:ko00002" = c("K00002", "K00003")
    )
  })

  # Mock the run_fgsea function to return a simple result
  mockery::stub(pathway_gsea, "run_fgsea", function(...) {
    data.frame(
      pathway_id = c("path:ko00001", "path:ko00002"),
      pathway_name = c("path:ko00001", "path:ko00002"),
      size = c(2, 2),
      ES = c(0.5, -0.3),
      NES = c(1.2, -0.8),
      pvalue = c(0.01, 0.05),
      p.adjust = c(0.02, 0.1),
      leading_edge = c("K00001;K00002", "K00003"),
      stringsAsFactors = FALSE
    )
  })

  # Test the function
  result <- pathway_gsea(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    method = "fgsea"
  )

  # Check the result
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(result$pathway_id, c("path:ko00001", "path:ko00002"))
  expect_equal(result$method, c("fgsea", "fgsea"))
})

test_that("pathway_gsea validates inputs correctly", {
  # Create test data
  test_data <- create_test_data()
  abundance <- test_data$abundance
  metadata <- test_data$metadata

  # Test invalid abundance
  expect_error(
    pathway_gsea(abundance = "invalid", metadata = metadata, group = "group"),
    "'abundance' must be a data frame or matrix"
  )

  # Test invalid metadata
  expect_error(
    pathway_gsea(abundance = abundance, metadata = "invalid", group = "group"),
    "'metadata' must be a data frame"
  )

  # Test invalid group
  expect_error(
    pathway_gsea(abundance = abundance, metadata = metadata, group = "invalid_group"),
    "Group variable invalid_group not found in metadata"
  )

  # Test invalid pathway_type
  expect_error(
    pathway_gsea(abundance = abundance, metadata = metadata, group = "group", pathway_type = "invalid"),
    "pathway_type must be one of 'KEGG', 'MetaCyc', or 'GO'"
  )

  # Test invalid method
  expect_error(
    pathway_gsea(abundance = abundance, metadata = metadata, group = "group", method = "invalid"),
    "method must be one of 'fgsea', 'GSEA', or 'clusterProfiler'"
  )

  # Test invalid rank_method
  expect_error(
    pathway_gsea(abundance = abundance, metadata = metadata, group = "group", rank_method = "invalid"),
    "rank_method must be one of 'signal2noise', 't_test', 'log2_ratio', or 'diff_abundance'"
  )
})

test_that("calculate_rank_metric works correctly", {
  # Create test data
  test_data <- create_test_data(n_features = 5, n_samples = 10)
  abundance <- test_data$abundance
  metadata <- test_data$metadata

  # Test signal2noise method
  metric_s2n <- calculate_rank_metric(abundance, metadata, "group", "signal2noise")
  expect_type(metric_s2n, "double")
  expect_length(metric_s2n, nrow(abundance))
  expect_named(metric_s2n, rownames(abundance))

  # Test t_test method
  metric_ttest <- calculate_rank_metric(abundance, metadata, "group", "t_test")
  expect_type(metric_ttest, "double")
  expect_length(metric_ttest, nrow(abundance))

  # Test log2_ratio method
  metric_log2 <- calculate_rank_metric(abundance, metadata, "group", "log2_ratio")
  expect_type(metric_log2, "double")
  expect_length(metric_log2, nrow(abundance))

  # Test diff_abundance method
  metric_diff <- calculate_rank_metric(abundance, metadata, "group", "diff_abundance")
  expect_type(metric_diff, "double")
  expect_length(metric_diff, nrow(abundance))

  # Test error with more than two groups
  metadata_multi <- metadata
  metadata_multi$group <- factor(rep(c("A", "B", "C"), length.out = nrow(metadata)))
  expect_error(
    calculate_rank_metric(abundance, metadata_multi, "group", "signal2noise"),
    "GSEA currently only supports two-group comparisons"
  )
})

test_that("prepare_gene_sets works correctly", {
  # Skip detailed testing of prepare_gene_sets since it depends on package data
  # Just test that it returns a list and issues warnings for unsupported pathway types

  # Test KEGG pathway - just check it returns a list with some elements
  gene_sets <- prepare_gene_sets("KEGG")
  expect_type(gene_sets, "list")
  expect_true(length(gene_sets) > 0)

  # Test MetaCyc pathway (not implemented yet)
  expect_warning(
    prepare_gene_sets("MetaCyc"),
    "MetaCyc pathway gene sets not yet implemented"
  )

  # Test GO pathway (not implemented yet)
  expect_warning(
    prepare_gene_sets("GO"),
    "GO pathway gene sets not yet implemented"
  )
})

test_that("run_fgsea works correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("fgsea")

  # Create test data
  ranked_list <- setNames(rnorm(100), paste0("K", sprintf("%05d", 1:100)))
  gene_sets <- list(
    "path:ko00001" = c("K00001", "K00002", "K00003"),
    "path:ko00002" = c("K00004", "K00005", "K00006")
  )

  # Mock the fgsea function
  mockery::stub(run_fgsea, "fgsea::fgsea", function(...) {
    data.frame(
      pathway = c("path:ko00001", "path:ko00002"),
      pval = c(0.01, 0.05),
      padj = c(0.02, 0.1),
      ES = c(0.5, -0.3),
      NES = c(1.2, -0.8),
      size = c(3, 3),
      leadingEdge = I(list(c("K00001", "K00002"), c("K00004"))),
      stringsAsFactors = FALSE
    )
  })

  # Test the function
  result <- run_fgsea(ranked_list, gene_sets)

  # Check the result
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(result$pathway_id, c("path:ko00001", "path:ko00002"))
  expect_equal(result$pvalue, c(0.01, 0.05))
  expect_equal(result$p.adjust, c(0.02, 0.1))
})




