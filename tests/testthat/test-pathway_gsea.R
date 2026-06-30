# Tests for pathway_gsea function

# Helper: create standard GSEA test data
create_gsea_test_data <- function(n_features = 50, n_samples = 10) {
  set.seed(123)
  abundance <- matrix(abs(rnorm(n_features * n_samples)) * 100,
                      nrow = n_features, ncol = n_samples)
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

test_that("run_fgsea honors the requested p-value adjustment method", {
  skip_if_not_installed("fgsea")

  run <- getFromNamespace("run_fgsea", "ggpicrust2")
  ranked_list <- c(
    K00001 = 5,
    K00002 = 4,
    K00003 = -4,
    K00004 = -5,
    K00005 = 1,
    K00006 = -1
  )
  gene_sets <- list(
    set1 = c("K00001", "K00002"),
    set2 = c("K00003", "K00004"),
    set3 = c("K00005", "K00006")
  )

  set.seed(1)
  result <- run(
    ranked_list,
    gene_sets,
    min_size = 2,
    max_size = 10,
    p_adjust_method = "bonferroni"
  )

  expect_gt(nrow(result), 0)
  expect_equal(
    result$p.adjust,
    stats::p.adjust(result$pvalue, method = "bonferroni")
  )
})

test_that("run_fgsea rejects invalid preranked statistic vectors before backend execution", {
  run <- getFromNamespace("run_fgsea", "ggpicrust2")
  gene_sets <- list(set1 = c("K00001", "K00002"))

  expect_error(
    run(c(1, 2), gene_sets, min_size = 1, max_size = 10),
    "one feature name per ranking statistic"
  )

  expect_error(
    run(c(K00001 = 1, K00001 = -1, K00002 = 0),
        gene_sets, min_size = 1, max_size = 10),
    "duplicated feature names"
  )

  expect_error(
    run(c(K00001 = 1, K00002 = Inf),
        gene_sets, min_size = 1, max_size = 10),
    "finite, non-missing"
  )

  expect_error(
    run(c(K00001 = 0, K00002 = 0),
        gene_sets, min_size = 1, max_size = 10),
    "at least two distinct"
  )
})

test_that("pathway_gsea rejects negative abundance values instead of coercing them to zero", {
  skip_if_not_installed("limma")

  test_data <- create_gsea_test_data()
  test_data$abundance[1, 1] <- -1

  expect_error(
    pathway_gsea(
      abundance = test_data$abundance,
      metadata = test_data$metadata,
      group = "group",
      method = "camera"
    ),
    "non-negative values"
  )

  calc <- getFromNamespace("calculate_rank_metric", "ggpicrust2")
  expect_error(
    calc(test_data$abundance, test_data$metadata, "group"),
    "non-negative values"
  )
})

test_that("pathway_gsea rejects MetaCyc pathway-level input for GSEA", {
  skip_if_not_installed("limma")

  test_data <- create_gsea_test_data(n_features = 6, n_samples = 6)
  rownames(test_data$abundance) <- c(
    "PWY-7219",
    "1CMET2-PWY",
    "EC:1.1.1.1",
    "EC:2.2.2.2",
    "EC:3.3.3.3",
    "EC:4.4.4.4"
  )

  expect_error(
    pathway_gsea(
      abundance = abs(test_data$abundance),
      metadata = test_data$metadata,
      group = "group",
      pathway_type = "MetaCyc",
      method = "camera",
      min_size = 1,
      max_size = 10
    ),
    "requires gene/enzyme-level EC abundance input"
  )
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
  expect_equal(unique(result$score_type), "signed_log10_pvalue")
  expect_equal(unique(result$score_label), "Signed -log10(p-value)")
  expect_equal(result$NES, result$signed_log10_pvalue)
})

test_that("run_limma_gsea stops when voom fails", {
  skip_if_not_installed("limma")

  test_data <- create_gsea_test_data(n_features = 6, n_samples = 6)
  run_limma <- getFromNamespace("run_limma_gsea", "ggpicrust2")

  local_mocked_bindings(
    voom = function(...) {
      stop("synthetic voom failure")
    },
    .package = "limma"
  )

  expect_error(
    run_limma(
      abundance_mat = test_data$abundance,
      metadata = test_data$metadata,
      group = "group",
      gene_sets = list(set1 = rownames(test_data$abundance)[1:3]),
      method = "camera",
      min_size = 2,
      max_size = 10
    ),
    paste0(
      "limma::voom\\(\\) failed, so camera gene-set testing was not run.*",
      "not a statistically equivalent fallback.*synthetic voom failure"
    )
  )
})

test_that("pathway_gsea validates inputs correctly", {
  test_data <- create_gsea_test_data()
  nonnegative_abundance <- abs(test_data$abundance)

  expect_error(pathway_gsea(abundance = "invalid", metadata = test_data$metadata, group = "group"),
               "'abundance' must be a data frame or matrix")
  expect_error(pathway_gsea(abundance = test_data$abundance, metadata = "invalid", group = "group"),
               "'metadata' must be a data frame")
  expect_error(pathway_gsea(abundance = test_data$abundance, metadata = test_data$metadata, group = "invalid_group"),
               "not found in metadata")
  expect_error(pathway_gsea(abundance = test_data$abundance, metadata = test_data$metadata, group = "group", method = "invalid"),
               "'method' must be one of")
  expect_error(
    pathway_gsea(abundance = nonnegative_abundance, metadata = test_data$metadata,
                 group = "group", min_size = 0),
    "'min_size' must be positive"
  )
  expect_error(
    pathway_gsea(abundance = nonnegative_abundance, metadata = test_data$metadata,
                 group = "group", min_size = 10, max_size = 5),
    "'min_size' must be less than or equal to 'max_size'"
  )
  expect_error(
    pathway_gsea(abundance = nonnegative_abundance, metadata = test_data$metadata,
                 group = "group", nperm = 1.5),
    "'nperm' must be a single finite integer"
  )
  expect_error(
    pathway_gsea(abundance = nonnegative_abundance, metadata = test_data$metadata,
                 group = "group", seed = NA_real_),
    "'seed' must be a single non-negative integer"
  )
  expect_error(
    pathway_gsea(abundance = nonnegative_abundance, metadata = test_data$metadata,
                 group = "group", p_adjust_method = "bogus"),
    "'p_adjust_method' must be one of"
  )
  expect_error(
    pathway_gsea(abundance = nonnegative_abundance, metadata = test_data$metadata,
                 group = "group", method = "camera", inter.gene.cor = c(0.1, 0.2)),
    "inter.gene.cor must be a numeric value between -1 and 1"
  )
})

test_that("pathway_gsea validates scalar choices and gene-set identifiers", {
  test_data <- create_gsea_test_data()
  nonnegative_abundance <- abs(test_data$abundance)

  expect_error(
    pathway_gsea(
      abundance = nonnegative_abundance,
      metadata = test_data$metadata,
      group = "group",
      method = c("camera", "fry")
    ),
    "'method' must be one of"
  )
  expect_error(
    pathway_gsea(
      abundance = nonnegative_abundance,
      metadata = test_data$metadata,
      group = "group",
      pathway_type = c("KEGG", "GO")
    ),
    "'pathway_type' must be one of"
  )
  expect_error(
    pathway_gsea(
      abundance = nonnegative_abundance,
      metadata = test_data$metadata,
      group = "group",
      method = "fgsea",
      rank_method = c("t_test", "log2_ratio")
    ),
    "'rank_method' must be one of"
  )
  expect_error(
    prepare_gene_sets(NULL),
    "'pathway_type' must be one of"
  )
  expect_error(
    prepare_gene_sets("GO", go_category = c("MF", "CC")),
    "'go_category' must be a single non-empty character value"
  )

  validate_gene_sets <- getFromNamespace("validate_gene_sets", "ggpicrust2")
  duplicated_set_names <- list(
    set1 = c("K00001", "K00002"),
    set1 = c("K00003", "K00004")
  )
  expect_error(
    validate_gene_sets(duplicated_set_names, "test gene_sets"),
    "duplicated gene-set names"
  )
  expect_error(
    validate_gene_sets(list(set1 = c("K00001", NA_character_)), "test gene_sets"),
    "must contain non-empty values"
  )

  sanitized <- validate_gene_sets(
    list(set1 = c("K00001", "K00001", "K00002")),
    "test gene_sets"
  )
  expect_equal(sanitized$set1, c("K00001", "K00002"))

  run_limma <- getFromNamespace("run_limma_gsea", "ggpicrust2")
  expect_error(
    run_limma(
      nonnegative_abundance,
      test_data$metadata,
      "group",
      gene_sets = duplicated_set_names,
      min_size = 2,
      max_size = 10
    ),
    "duplicated gene-set names"
  )
})

test_that("pathway_gsea rejects nonnumeric matrices and missing feature rownames", {
  skip_if_not_installed("limma")

  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  no_feature_ids <- matrix(
    seq_len(8 * 4),
    nrow = 8,
    dimnames = list(NULL, paste0("S", 1:4))
  )
  expect_error(
    pathway_gsea(
      abundance = no_feature_ids,
      metadata = metadata,
      group = "group",
      method = "camera"
    ),
    "feature identifiers|row names"
  )

  nonnumeric <- matrix(
    as.character(seq_len(8 * 4)),
    nrow = 8,
    dimnames = list(paste0("K", sprintf("%05d", 1:8)), paste0("S", 1:4))
  )
  expect_error(
    pathway_gsea(
      abundance = nonnumeric,
      metadata = metadata,
      group = "group",
      method = "camera"
    ),
    "must be numeric"
  )
})

test_that("pathway_gsea revalidates groups and design variables after alignment", {
  skip_if_not_installed("limma")

  abundance <- matrix(
    seq_len(10 * 4),
    nrow = 10,
    dimnames = list(paste0("K", sprintf("%05d", 1:10)), paste0("S", 1:4))
  )

  metadata_single_after_align <- data.frame(
    sample = paste0("S", 1:6),
    group = c("A", "A", "A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  expect_error(
    pathway_gsea(
      abundance = abundance,
      metadata = metadata_single_after_align,
      group = "group",
      method = "camera"
    ),
    "after sample alignment"
  )

  metadata_group_na <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", NA),
    stringsAsFactors = FALSE
  )
  expect_error(
    pathway_gsea(
      abundance = abundance,
      metadata = metadata_group_na,
      group = "group",
      method = "camera"
    ),
    "contains missing or empty values after sample alignment"
  )

  metadata_cov_na <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    age = c(20, NA, 30, 31),
    stringsAsFactors = FALSE
  )
  expect_error(
    pathway_gsea(
      abundance = abundance,
      metadata = metadata_cov_na,
      group = "group",
      method = "camera",
      covariates = "age"
    ),
    "Design variables contain missing values"
  )

  build_design <- getFromNamespace("build_design_matrix", "ggpicrust2")
  expect_error(
    build_design(metadata_cov_na, "group", covariates = "age"),
    "Design variables contain missing values"
  )
})

test_that("calculate_rank_metric rejects missing group assignments", {
  calc <- getFromNamespace("calculate_rank_metric", "ggpicrust2")

  abundance <- matrix(
    c(1, 2, 10, 11,
      3, 4, 12, 13),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("K00001", "K00002"), paste0("S", 1:4))
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", NA),
    stringsAsFactors = FALSE
  )

  expect_error(
    calc(abundance, metadata, "group", method = "signal2noise"),
    "contains missing or empty values after sample alignment"
  )

  abundance_with_extra_missing <- matrix(
    c(1, 2, 10, 11, 5,
      3, 4, 12, 13, 6),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("K00001", "K00002"), paste0("S", 1:5))
  )
  metadata_with_extra_missing <- data.frame(
    sample = paste0("S", 1:5),
    group = c("A", "A", "B", "B", NA),
    stringsAsFactors = FALSE
  )

  expect_error(
    calc(
      abundance_with_extra_missing,
      metadata_with_extra_missing,
      "group",
      method = "diff_abundance",
      comparison = c("B", "A")
    ),
    "contains missing or empty values after sample alignment"
  )
})

test_that("calculate_rank_metric requires explicit feature IDs and non-tied ranks", {
  calc <- getFromNamespace("calculate_rank_metric", "ggpicrust2")
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  no_feature_ids <- matrix(
    c(1, 2, 10, 11,
      3, 4, 12, 13),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, paste0("S", 1:4))
  )
  expect_error(
    calc(no_feature_ids, metadata, "group", method = "diff_abundance"),
    "feature identifiers|row names"
  )

  tied_abundance <- matrix(
    5,
    nrow = 3,
    ncol = 4,
    dimnames = list(paste0("K0000", 1:3), paste0("S", 1:4))
  )
  expect_error(
    calc(tied_abundance, metadata, "group", method = "diff_abundance"),
    "at least two distinct"
  )
})

test_that("limma GSEA contrast resolution is explicit and exact", {
  build_design <- getFromNamespace("build_design_matrix", "ggpicrust2")
  resolve_contrast <- getFromNamespace("resolve_limma_contrast", "ggpicrust2")

  metadata_two_group <- data.frame(
    sample = paste0("S", 1:6),
    group = factor(rep(c("A", "B"), each = 3)),
    batch = factor(rep(c("batch1", "B"), 3)),
    stringsAsFactors = FALSE
  )
  design_two_group <- build_design(metadata_two_group, "group", covariates = "batch")

  expect_equal(
    resolve_contrast(design_two_group, metadata_two_group, "group", NULL),
    match("groupB", colnames(design_two_group))
  )
  expect_equal(
    resolve_contrast(design_two_group, metadata_two_group, "group", "B"),
    match("groupB", colnames(design_two_group))
  )
  covariate_col <- setdiff(colnames(design_two_group), c("Intercept", "groupB"))
  expect_equal(
    resolve_contrast(design_two_group, metadata_two_group, "group", covariate_col),
    match(covariate_col, colnames(design_two_group))
  )

  metadata_no_b_level <- metadata_two_group
  metadata_no_b_level$group <- factor(rep(c("Control", "Treatment"), each = 3))
  design_no_b_level <- build_design(metadata_no_b_level, "group", covariates = "batch")
  expect_error(
    resolve_contrast(design_no_b_level, metadata_no_b_level, "group", "B"),
    "exact design column or non-reference group level"
  )

  metadata_three_group <- data.frame(
    sample = paste0("S", 1:9),
    group = factor(rep(c("A", "B", "C"), each = 3)),
    stringsAsFactors = FALSE
  )
  design_three_group <- build_design(metadata_three_group, "group")
  expect_error(
    resolve_contrast(design_three_group, metadata_three_group, "group", NULL),
    "Specify 'contrast' explicitly"
  )
  expect_equal(
    resolve_contrast(design_three_group, metadata_three_group, "group", "C"),
    match("groupC", colnames(design_three_group))
  )

  expect_error(
    resolve_contrast(design_three_group, metadata_three_group, "group", "A"),
    "reference level"
  )
  expect_error(
    resolve_contrast(design_three_group, metadata_three_group, "group", c(0, 1)),
    "length equal to the number of design columns"
  )
  expect_equal(
    resolve_contrast(design_three_group, metadata_three_group, "group", c(0, -1, 1)),
    c(0, -1, 1)
  )

  named_out_of_order <- c(groupC = 1, groupB = -1, Intercept = 0)
  expect_equal(
    resolve_contrast(design_three_group, metadata_three_group, "group", named_out_of_order),
    c(Intercept = 0, groupB = -1, groupC = 1)
  )
  expect_error(
    resolve_contrast(
      design_three_group,
      metadata_three_group,
      "group",
      c(Intercept = 0, groupB = -1, wrong = 1)
    ),
    "names must match the design columns exactly"
  )
})

test_that("limma GSEA design supports non-syntactic group and covariate names", {
  build_design <- getFromNamespace("build_design_matrix", "ggpicrust2")
  resolve_contrast <- getFromNamespace("resolve_limma_contrast", "ggpicrust2")

  metadata <- data.frame(
    sample = paste0("S", 1:6),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  metadata[["treatment group"]] <- factor(
    rep(c("baseline group", "treated group"), each = 3)
  )
  metadata[["age years"]] <- c(20, 23, 26, 21, 24, 27)

  design <- build_design(
    metadata,
    "treatment group",
    covariates = "age years"
  )

  expect_equal(
    colnames(design),
    c("Intercept", "`treatment group`treated group", "`age years`")
  )
  expect_equal(
    resolve_contrast(design, metadata, "treatment group", NULL),
    match("`treatment group`treated group", colnames(design))
  )
  expect_equal(
    resolve_contrast(design, metadata, "treatment group", "treated group"),
    match("`treatment group`treated group", colnames(design))
  )
})

test_that("pathway_gsea camera handles non-syntactic design variables", {
  skip_if_not_installed("limma")

  test_data <- create_gsea_test_data(n_features = 30, n_samples = 6)
  abundance <- abs(test_data$abundance) * 100
  metadata <- data.frame(
    sample_name = colnames(abundance),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  metadata[["treatment group"]] <- factor(
    rep(c("baseline group", "treated group"), each = 3)
  )
  metadata[["age years"]] <- c(20, 23, 26, 21, 24, 27)

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
    metadata = metadata,
    group = "treatment group",
    covariates = "age years",
    method = "camera",
    pathway_type = "KEGG",
    contrast = "treated group"
  )

  expect_s3_class(result, "data.frame")
  expect_equal(unique(result$method), "camera")
  expect_true(all(c("pathway_id", "direction", "pvalue") %in% colnames(result)))
})

test_that("prepare_gene_sets works for KEGG and MetaCyc pathway types", {
  gene_sets_kegg <- prepare_gene_sets("KEGG")
  expect_type(gene_sets_kegg, "list")
  expect_true(length(gene_sets_kegg) > 0)
  expect_false("ko01001" %in% names(gene_sets_kegg))
  expect_false("ko99980" %in% names(gene_sets_kegg))

  gene_sets_metacyc <- prepare_gene_sets("MetaCyc")
  expect_type(gene_sets_metacyc, "list")
})

test_that("calculate_rank_metric aligns by sample column and handles constant t-test rows", {
  calc <- getFromNamespace("calculate_rank_metric", "ggpicrust2")
  abundance <- matrix(
    c(
      1, 1, 1, 1,
      1, 2, 8, 9
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("K00001", "K00002"), paste0("S", 1:4))
  )
  metadata <- data.frame(
    sample_name = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  metric <- calc(abundance, metadata, "group", method = "t_test")
  expect_named(metric, rownames(abundance))
  expect_equal(unname(metric["K00001"]), 0)
  expect_true(is.finite(metric["K00002"]))
})

test_that("calculate_rank_metric uses explicit preranked comparison direction", {
  calc <- getFromNamespace("calculate_rank_metric", "ggpicrust2")

  abundance <- matrix(
    c(
      10, 11, 2, 3, 5, 5,
      1, 2, 12, 13, 5, 6
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("K00001", "K00002"), paste0("S", 1:6))
  )
  metadata <- data.frame(
    sample_name = paste0("S", 1:6),
    group = factor(
      c("A", "A", "B", "B", "C", "C"),
      levels = c("A", "B", "C")
    ),
    stringsAsFactors = FALSE
  )

  expect_error(
    calc(abundance, metadata, "group", method = "diff_abundance"),
    "requires exactly two group levels"
  )

  metric <- calc(
    abundance,
    metadata,
    "group",
    method = "diff_abundance",
    comparison = c("B", "A")
  )
  expect_equal(attr(metric, "comparison"), c("B", "A"))
  expect_gt(metric["K00002"], 0)
  expect_lt(metric["K00001"], 0)

  reversed <- calc(
    abundance,
    metadata,
    "group",
    method = "diff_abundance",
    comparison = c("A", "B")
  )
  expect_equal(unname(reversed["K00002"]), -unname(metric["K00002"]))
})

test_that("pathway_gsea passes explicit preranked comparison to fgsea", {
  skip_if_not_installed("fgsea")

  abundance <- matrix(
    c(
      1, 2, 20, 21,
      20, 21, 1, 2,
      5, 5, 5, 5,
      6, 6, 6, 6,
      7, 7, 7, 7,
      8, 8, 8, 8,
      9, 9, 9, 9,
      10, 10, 10, 10
    ),
    nrow = 8,
    byrow = TRUE,
    dimnames = list(paste0("K", sprintf("%05d", 1:8)), paste0("S", 1:4))
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = factor(c("Control", "Control", "Treatment", "Treatment"),
                   levels = c("Control", "Treatment")),
    stringsAsFactors = FALSE
  )
  captured_ranked_list <- NULL

  local_mocked_bindings(
    prepare_gene_sets = function(...) {
      list("ko00010" = rownames(abundance))
    },
    run_fgsea = function(ranked_list, ...) {
      captured_ranked_list <<- ranked_list
      data.frame(
        pathway_id = "ko00010",
        pathway_name = "ko00010",
        size = 8,
        ES = 0.5,
        NES = 1.2,
        pvalue = 0.01,
        p.adjust = 0.02,
        leading_edge = "K00001",
        stringsAsFactors = FALSE
      )
    }
  )

  result <- pathway_gsea(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    method = "fgsea",
    rank_method = "diff_abundance",
    comparison = c("Treatment", "Control")
  )

  expect_gt(captured_ranked_list["K00001"], 0)
  expect_lt(captured_ranked_list["K00002"], 0)
  expect_equal(result$group1, "Treatment")
  expect_equal(result$group2, "Control")
})

test_that("pathway_gsea rejects missing groups with explicit preranked comparison", {
  skip_if_not_installed("fgsea")

  abundance <- matrix(
    rep(seq_len(8), 5),
    nrow = 8,
    dimnames = list(paste0("K", sprintf("%05d", 1:8)), paste0("S", 1:5))
  )
  metadata <- data.frame(
    sample = paste0("S", 1:5),
    group = c("A", "A", "B", "B", NA),
    stringsAsFactors = FALSE
  )

  expect_error(
    pathway_gsea(
      abundance = abundance,
      metadata = metadata,
      group = "group",
      method = "fgsea",
      rank_method = "diff_abundance",
      comparison = c("B", "A")
    ),
    "contains missing or empty values after sample alignment"
  )
})

test_that("pathway_gsea rejects comparison for limma-based GSEA methods", {
  test_data <- create_gsea_test_data()

  expect_error(
    pathway_gsea(
      abundance = test_data$abundance,
      metadata = test_data$metadata,
      group = "group",
      method = "camera",
      comparison = c("Treatment", "Control")
    ),
    "only supported for preranked GSEA methods"
  )
})

test_that("pathway_gsea rejects design-only arguments for preranked methods", {
  test_data <- create_gsea_test_data()
  metadata <- test_data$metadata
  metadata$age <- seq_len(nrow(metadata)) + 20

  expect_error(
    pathway_gsea(
      abundance = test_data$abundance,
      metadata = metadata,
      group = "group",
      method = "fgsea",
      covariates = "age"
    ),
    "Covariates are only supported for 'camera' and 'fry'"
  )

  expect_error(
    pathway_gsea(
      abundance = test_data$abundance,
      metadata = metadata,
      group = "group",
      method = "fgsea",
      contrast = "groupTreatment"
    ),
    "'contrast' is only supported for limma-based GSEA methods"
  )
})

test_that("build_design_matrix validates covariate names strictly", {
  build_design <- getFromNamespace("build_design_matrix", "ggpicrust2")
  metadata <- data.frame(
    sample = paste0("S", 1:6),
    group = rep(c("A", "B"), each = 3),
    age = seq_len(6) + 20,
    batch = rep(c("X", "Y"), 3),
    stringsAsFactors = FALSE
  )

  expect_error(
    build_design(metadata, "group", covariates = c("age", NA_character_)),
    "'covariates' must be a character vector"
  )
  expect_error(
    build_design(metadata, "group", covariates = c("age", "age")),
    "'covariates' must contain unique column names"
  )
  expect_error(
    build_design(metadata, "group", covariates = "group"),
    "must not include the group column"
  )
  expect_error(
    build_design(metadata, c("group", "batch"), covariates = "age"),
    "'group' must be a single non-empty character string"
  )
})

test_that("build_design_matrix rejects non-finite or rank-deficient designs", {
  build_design <- getFromNamespace("build_design_matrix", "ggpicrust2")
  metadata <- data.frame(
    sample = paste0("S", 1:6),
    group = rep(c("A", "B"), each = 3),
    age = seq_len(6) + 20,
    stringsAsFactors = FALSE
  )

  metadata_constant <- metadata
  metadata_constant$constant_covariate <- 1
  expect_error(
    build_design(metadata_constant, "group", covariates = "constant_covariate"),
    "rank deficient"
  )

  metadata_confounded <- metadata
  metadata_confounded$batch <- metadata_confounded$group
  expect_error(
    build_design(metadata_confounded, "group", covariates = "batch"),
    "rank deficient"
  )

  metadata_infinite <- metadata
  metadata_infinite$age[1] <- Inf
  expect_error(
    build_design(metadata_infinite, "group", covariates = "age"),
    "non-finite values"
  )
})

test_that("pathway_gsea validates PICRUSt2 #NAME input after stripping the ID column", {
  skip_if_not_installed("fgsea")

  abundance <- data.frame(
    `#NAME` = paste0("K", sprintf("%05d", 1:8)),
    S1 = c(1:8),
    S2 = c(2:9),
    S3 = c(10:17),
    S4 = c(11:18),
    check.names = FALSE
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    prepare_gene_sets = function(...) {
      list("ko00010" = abundance$`#NAME`)
    },
    run_fgsea = function(...) {
      data.frame(
        pathway_id = "ko00010",
        pathway_name = "ko00010",
        size = 8,
        ES = 0.5,
        NES = 1.2,
        pvalue = 0.01,
        p.adjust = 0.02,
        leading_edge = "K00001",
        stringsAsFactors = FALSE
      )
    }
  )

  result <- pathway_gsea(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    method = "fgsea"
  )
  expect_s3_class(result, "data.frame")
  expect_equal(result$pathway_id, "ko00010")
})

test_that("clusterProfiler GSEA keeps the standard schema for empty results", {
  skip_if_not_installed("clusterProfiler")

  abundance <- matrix(
    c(
      2, 2, 1, 1,
      1, 1, 2, 2,
      1, 1, 1, 1,
      2, 2, 2, 2
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("K", sprintf("%05d", 1:4)), paste0("S", 1:4))
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    prepare_gene_sets = function(...) {
      list(set1 = c("K00001", "K00002"))
    }
  )

  expect_warning(
    result <- pathway_gsea(
      abundance = abundance,
      metadata = metadata,
      group = "group",
      method = "clusterProfiler",
      min_size = 3,
      max_size = 10
    ),
    "No gene sets overlapped the ranked feature list"
  )

  expect_equal(nrow(result), 0)
  expect_true(all(c(
    "pathway_id", "pathway_name", "size", "ES", "NES",
    "pvalue", "p.adjust", "leading_edge", "method"
  ) %in% colnames(result)))
})

test_that("pathway_gsea accepts a generic leading feature ID column", {
  skip_if_not_installed("fgsea")

  abundance <- data.frame(
    feature = paste0("K", sprintf("%05d", 1:8)),
    S1 = c(1:8),
    S2 = c(2:9),
    S3 = c(10:17),
    S4 = c(11:18),
    check.names = FALSE
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    prepare_gene_sets = function(...) {
      list("ko00010" = abundance$feature)
    },
    run_fgsea = function(...) {
      data.frame(
        pathway_id = "ko00010",
        pathway_name = "ko00010",
        size = 8,
        ES = 0.5,
        NES = 1.2,
        pvalue = 0.01,
        p.adjust = 0.02,
        leading_edge = "K00001",
        stringsAsFactors = FALSE
      )
    }
  )

  result <- pathway_gsea(
    abundance = abundance,
    metadata = metadata,
    group = "group",
    method = "fgsea"
  )
  expect_s3_class(result, "data.frame")
  expect_equal(result$pathway_id, "ko00010")
})

test_that("prepare_gene_sets warns when `organism` is set to a non-default value but still returns the same KO gene sets", {
  # Regression: `organism` was documented on both `prepare_gene_sets()`
  # and `pathway_gsea()` but the KEGG and GO branches both read
  # KO-based reference tables (`ko_to_kegg_reference`,
  # `ko_to_go_reference`) that are organism-independent by
  # construction. A user passing `organism = "hsa"` silently got the
  # same gene sets as the default -- a promise/implementation gap --
  # which we now surface as a deprecation warning. The gene-set
  # content must remain identical across organism values (the
  # argument truly has no effect), both for KEGG and GO.
  expect_warning(
    gs_hsa <- prepare_gene_sets("KEGG", organism = "hsa"),
    regexp = "organism.*deprecated.*no effect"
  )
  gs_default <- suppressWarnings(prepare_gene_sets("KEGG", organism = "ko"))
  expect_identical(gs_hsa, gs_default)

  # Default value does not warn.
  expect_warning(
    prepare_gene_sets("KEGG", organism = "ko"),
    regexp = NA
  )

  # GO branch: same contract.
  expect_warning(
    go_hsa <- prepare_gene_sets("GO", organism = "hsa"),
    regexp = "organism.*deprecated.*no effect"
  )
  go_default <- suppressWarnings(prepare_gene_sets("GO", organism = "ko"))
  expect_identical(go_hsa, go_default)
})
