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
    "method must be one of:"
  )

  # Test invalid rank_method (only for preranked methods)
  expect_error(
    pathway_gsea(abundance = abundance, metadata = metadata, group = "group",
                 method = "fgsea", rank_method = "invalid"),
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
  # Just test that it returns a list for all supported pathway types

  # Test KEGG pathway - just check it returns a list with some elements
  gene_sets <- prepare_gene_sets("KEGG")
  expect_type(gene_sets, "list")
  expect_true(length(gene_sets) > 0)

  # Test MetaCyc pathway - now supported
  gene_sets_metacyc <- prepare_gene_sets("MetaCyc")
  expect_type(gene_sets_metacyc, "list")
  # MetaCyc should return gene sets (may be empty if no data)

  # Test GO pathway - now supported
  gene_sets_go <- prepare_gene_sets("GO")
  expect_type(gene_sets_go, "list")
  # GO should return gene sets (may be empty if no data)
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

# ============================================================================
# Tests for limma-based methods (camera/fry) - Issue #193
# ============================================================================

test_that("build_design_matrix works correctly", {
  # Create test metadata
  metadata <- data.frame(
    sample = paste0("S", 1:10),
    group = factor(rep(c("Control", "Treatment"), each = 5)),
    age = c(25, 30, 35, 28, 32, 27, 33, 29, 31, 26),
    sex = factor(rep(c("M", "F"), 5)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample

  # Test simple design (no covariates)
  design_simple <- build_design_matrix(metadata, "group")
  expect_true(is.matrix(design_simple))
  expect_equal(nrow(design_simple), 10)
  expect_equal(ncol(design_simple), 2)  # Intercept + group
  expect_true("Intercept" %in% colnames(design_simple))

  # Test design with covariates
  design_with_cov <- build_design_matrix(metadata, "group", covariates = c("age", "sex"))
  expect_true(is.matrix(design_with_cov))
  expect_equal(nrow(design_with_cov), 10)
  expect_true(ncol(design_with_cov) >= 3)  # Intercept + group + age + sex
})

test_that("run_limma_gsea validates inputs correctly", {
  skip_if_not_installed("limma")

  # Create minimal test data
  set.seed(123)
  abundance_mat <- matrix(abs(rnorm(50 * 10)), nrow = 50, ncol = 10)
  colnames(abundance_mat) <- paste0("S", 1:10)
  rownames(abundance_mat) <- paste0("K", sprintf("%05d", 1:50))

  metadata <- data.frame(
    sample = paste0("S", 1:10),
    group = factor(rep(c("Control", "Treatment"), each = 5)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample

  gene_sets <- list(
    "pathway1" = paste0("K", sprintf("%05d", 1:10)),
    "pathway2" = paste0("K", sprintf("%05d", 11:20))
  )

  # Test that it runs without error
  result <- run_limma_gsea(
    abundance_mat = abundance_mat,
    metadata = metadata,
    group = "group",
    gene_sets = gene_sets,
    method = "camera",
    min_size = 5,
    max_size = 100
  )

  expect_s3_class(result, "data.frame")
  expect_true("pathway_id" %in% colnames(result))
  expect_true("direction" %in% colnames(result))
  expect_true("pvalue" %in% colnames(result))
  expect_true("method" %in% colnames(result))
})

test_that("camera method works with covariates", {
  skip_if_not_installed("limma")

  # Create test data with covariates
  set.seed(456)
  n_features <- 100
  n_samples <- 20

  abundance_mat <- matrix(abs(rnorm(n_features * n_samples) * 100),
                          nrow = n_features, ncol = n_samples)
  colnames(abundance_mat) <- paste0("S", 1:n_samples)
  rownames(abundance_mat) <- paste0("K", sprintf("%05d", 1:n_features))

  metadata <- data.frame(
    sample = paste0("S", 1:n_samples),
    group = factor(rep(c("Control", "Treatment"), each = n_samples/2)),
    age = rnorm(n_samples, mean = 50, sd = 10),
    sex = factor(sample(c("M", "F"), n_samples, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample

  gene_sets <- list(
    "pathway1" = paste0("K", sprintf("%05d", 1:15)),
    "pathway2" = paste0("K", sprintf("%05d", 16:30)),
    "pathway3" = paste0("K", sprintf("%05d", 31:50))
  )

  # Test camera without covariates
  result_no_cov <- run_limma_gsea(
    abundance_mat = abundance_mat,
    metadata = metadata,
    group = "group",
    covariates = NULL,
    gene_sets = gene_sets,
    method = "camera",
    min_size = 5,
    max_size = 100
  )

  expect_s3_class(result_no_cov, "data.frame")
  expect_true(nrow(result_no_cov) > 0)
  expect_equal(unique(result_no_cov$method), "camera")

  # Test camera with covariates
  result_with_cov <- run_limma_gsea(
    abundance_mat = abundance_mat,
    metadata = metadata,
    group = "group",
    covariates = c("age", "sex"),
    gene_sets = gene_sets,
    method = "camera",
    min_size = 5,
    max_size = 100
  )

  expect_s3_class(result_with_cov, "data.frame")
  expect_true(nrow(result_with_cov) > 0)
  expect_equal(unique(result_with_cov$method), "camera")

  # Results may differ when adjusting for covariates
  # (but this is expected behavior, not a test failure)
})

test_that("fry method works correctly", {
  skip_if_not_installed("limma")

  # Create test data
  set.seed(789)
  n_features <- 100
  n_samples <- 16

  abundance_mat <- matrix(abs(rnorm(n_features * n_samples) * 100),
                          nrow = n_features, ncol = n_samples)
  colnames(abundance_mat) <- paste0("S", 1:n_samples)
  rownames(abundance_mat) <- paste0("K", sprintf("%05d", 1:n_features))

  metadata <- data.frame(
    sample = paste0("S", 1:n_samples),
    group = factor(rep(c("Control", "Treatment"), each = n_samples/2)),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample

  gene_sets <- list(
    "pathway1" = paste0("K", sprintf("%05d", 1:20)),
    "pathway2" = paste0("K", sprintf("%05d", 21:40))
  )

  # Test fry method
  result <- run_limma_gsea(
    abundance_mat = abundance_mat,
    metadata = metadata,
    group = "group",
    gene_sets = gene_sets,
    method = "fry",
    min_size = 5,
    max_size = 100
  )

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_equal(unique(result$method), "fry")
  expect_true(all(result$direction %in% c("Up", "Down")))
})

test_that("pathway_gsea with camera method works end-to-end", {
  skip_if_not_installed("limma")

  # Create test data
  set.seed(101)
  test_data <- create_test_data(n_features = 50, n_samples = 12)
  abundance <- abs(test_data$abundance) * 100  # Make positive for counts
  metadata <- test_data$metadata

  # Mock prepare_gene_sets to return testable gene sets
  mockery::stub(pathway_gsea, "prepare_gene_sets", function(...) {
    list(
      "ko00010" = paste0("K", sprintf("%05d", 1:10)),
      "ko00020" = paste0("K", sprintf("%05d", 11:20)),
      "ko00030" = paste0("K", sprintf("%05d", 21:30))
    )
  })

  # Test with camera method
  result <- pathway_gsea(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    method = "camera",
    pathway_type = "KEGG"
  )

  expect_s3_class(result, "data.frame")
  expect_true("pathway_id" %in% colnames(result))
  expect_true("direction" %in% colnames(result))
  expect_true("pvalue" %in% colnames(result))
  expect_true("p.adjust" %in% colnames(result))
  expect_true("method" %in% colnames(result))
  expect_equal(unique(result$method), "camera")
})

test_that("pathway_gsea warns about covariates for preranked methods", {
  skip_if_not_installed("fgsea")

  test_data <- create_test_data()
  abundance <- test_data$abundance
  metadata <- test_data$metadata
  metadata$age <- rnorm(nrow(metadata))

  # Mock to avoid actual computation
  mockery::stub(pathway_gsea, "prepare_gene_sets", function(...) {
    list("ko00010" = c("K00001", "K00002"))
  })
  mockery::stub(pathway_gsea, "run_fgsea", function(...) {
    data.frame(
      pathway_id = "ko00010",
      pathway_name = "ko00010",
      size = 2, ES = 0.5, NES = 1.2,
      pvalue = 0.01, p.adjust = 0.02,
      leading_edge = "K00001",
      stringsAsFactors = FALSE
    )
  })

  # Should warn about covariates being ignored for fgsea
  expect_warning(
    pathway_gsea(
      abundance = abundance,
      metadata = metadata,
      group = "group",
      method = "fgsea",
      covariates = c("age")
    ),
    "Covariates are only supported for 'camera' and 'fry' methods"
  )
})

test_that("pathway_gsea validates covariate existence", {
  test_data <- create_test_data()
  abundance <- test_data$abundance
  metadata <- test_data$metadata

  # Should error when covariate doesn't exist
  expect_error(
    pathway_gsea(
      abundance = abundance,
      metadata = metadata,
      group = "group",
      method = "camera",
      covariates = c("nonexistent_covariate")
    ),
    "Covariates not found in metadata"
  )
})

test_that("inter.gene.cor parameter validation works", {
  test_data <- create_test_data()
  abundance <- test_data$abundance
  metadata <- test_data$metadata

  # Should error for invalid inter.gene.cor
  expect_error(
    pathway_gsea(
      abundance = abundance,
      metadata = metadata,
      group = "group",
      method = "camera",
      inter.gene.cor = 1.5  # Invalid: > 1
    ),
    "inter.gene.cor must be a numeric value between 0 and 1"
  )

  expect_error(
    pathway_gsea(
      abundance = abundance,
      metadata = metadata,
      group = "group",
      method = "camera",
      inter.gene.cor = -0.1  # Invalid: < 0
    ),
    "inter.gene.cor must be a numeric value between 0 and 1"
  )
})
