# Tests for pathway_errorbar_table

test_that("pathway_errorbar_table auto-detects the sample column from metadata", {
  # Regression: pathway_errorbar_table() used to hardcode
  # `sample_col = "sample_name"`, so users whose metadata used the far more
  # common `sample` column got an immediate
  # "Column 'sample_name' not found in metadata". Align behavior with the
  # rest of the package (align_samples / find_sample_column) by
  # auto-detecting the sample column when `sample_col = NULL`.

  abundance <- matrix(
    c(10, 20, 30, 40,
      15, 25, 35, 45,
      12, 22, 32, 42),
    nrow = 3, byrow = TRUE,
    dimnames = list(paste0("p", 1:3), paste0("S", 1:4))
  )

  metadata <- data.frame(
    sample = paste0("S", 1:4),
    Env    = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  daa_results_df <- data.frame(
    feature   = paste0("p", 1:3),
    method    = "ALDEx2",
    group1    = "A",
    group2    = "B",
    p_values  = c(0.01, 0.2, 0.04),
    p_adjust  = c(0.03, 0.2, 0.04),
    stringsAsFactors = FALSE
  )

  # Sample column auto-detection should succeed without passing sample_col.
  res <- pathway_errorbar_table(
    abundance      = abundance,
    daa_results_df = daa_results_df,
    Group          = metadata$Env,
    metadata       = metadata,
    p_values_threshold = 0.05
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("feature", "group1", "group2",
                    "mean_rel_abundance_group1",
                    "mean_rel_abundance_group2",
                    "p_adjust") %in% colnames(res)))
  # Only the features passing p_values_threshold should appear.
  expect_setequal(res$feature, c("p1", "p3"))
})

test_that("pathway_errorbar_table honors an explicit sample_col", {
  abundance <- matrix(
    c(10, 20, 30, 40,
      15, 25, 35, 45,
      12, 22, 32, 42),
    nrow = 3, byrow = TRUE,
    dimnames = list(paste0("p", 1:3), paste0("S", 1:4))
  )

  # Non-standard column name -- user must pass it explicitly.
  metadata <- data.frame(
    my_sample_id = paste0("S", 1:4),
    Env          = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  daa_results_df <- data.frame(
    feature   = paste0("p", 1:3),
    method    = "ALDEx2",
    group1    = "A",
    group2    = "B",
    p_values  = c(0.01, 0.2, 0.04),
    p_adjust  = c(0.03, 0.2, 0.04),
    stringsAsFactors = FALSE
  )

  res <- pathway_errorbar_table(
    abundance      = abundance,
    daa_results_df = daa_results_df,
    Group          = metadata$Env,
    metadata       = metadata,
    sample_col     = "my_sample_id",
    p_values_threshold = 0.05
  )

  expect_s3_class(res, "data.frame")
  expect_true("feature" %in% colnames(res))
})

test_that("pathway_errorbar_table aligns named Group vector to abundance columns", {
  abundance <- matrix(
    c(1, 1, 9, 9,
      9, 9, 1, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("f1", "f2"), paste0("S", 1:4))
  )
  group <- c(S1 = "A", S2 = "A", S3 = "B", S4 = "B")
  group_reordered <- group[rev(names(group))]

  daa_results_df <- data.frame(
    feature = c("f1", "f2"),
    method = "method",
    group1 = "A",
    group2 = "B",
    p_adjust = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  res <- pathway_errorbar_table(
    abundance = abundance,
    daa_results_df = daa_results_df,
    Group = group_reordered,
    p_values_threshold = 0.05
  )

  expect_gt(res$log2_fold_change[res$feature == "f1"], 0)
  expect_lt(res$log2_fold_change[res$feature == "f2"], 0)
})

test_that("pathway_errorbar_table rejects multi-contrast input instead of silently truncating", {
  # Regression: the previous implementation derived the two group names via
  #   group1_name <- unique(daa_results_filtered_sub_df$group1)[1]
  #   group2_name <- unique(daa_results_filtered_sub_df$group2)[1]
  # which LOOKS like "silently pick the first of many" behavior. In reality
  # the upstream validate_daa_results() was supposed to guarantee a single
  # pair, but the idiom would mask any future validator bypass by quietly
  # collapsing a multi-contrast DAA result down to its first (group1, group2)
  # block -- a dangerous "default success" where the returned table reflects
  # only one contrast out of several without any warning.
  #
  # The multi-contrast shape is exactly what pathway_daa() now emits for
  # >=3 groups (one row-block per non-reference level), so this is a real
  # user input, not a synthetic edge case.
  abundance <- matrix(
    c(10, 20, 30, 40, 50, 60,
      15, 25, 35, 45, 55, 65,
      12, 22, 32, 42, 52, 62),
    nrow = 3, byrow = TRUE,
    dimnames = list(paste0("p", 1:3), paste0("S", 1:6))
  )

  Group <- c("A", "A", "B", "B", "C", "C")

  # Two contrasts against reference "A": (A vs B) and (A vs C).
  daa_results_df <- data.frame(
    feature  = rep(paste0("p", 1:3), times = 2),
    method   = "DESeq2",
    group1   = "A",
    group2   = rep(c("B", "C"), each = 3),
    p_values = c(0.01, 0.2, 0.04, 0.03, 0.5, 0.02),
    p_adjust = c(0.02, 0.2, 0.05, 0.04, 0.5, 0.03),
    stringsAsFactors = FALSE
  )

  # Must hard-fail, not silently return a table for only the first contrast.
  expect_error(
    pathway_errorbar_table(
      abundance      = abundance,
      daa_results_df = daa_results_df,
      Group          = Group,
      p_values_threshold = 0.05
    ),
    regexp = "multiple group pairs"
  )
})

test_that("pathway_errorbar_table rejects duplicated features in one DAA contrast", {
  abundance <- matrix(
    c(10, 20, 30, 40,
      15, 25, 35, 45),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("p1", "p2"), paste0("S", 1:4))
  )
  Group <- c("A", "A", "B", "B")
  daa_results_df <- data.frame(
    feature = c("p1", "p1", "p2"),
    method = "DESeq2",
    group1 = "A",
    group2 = "B",
    p_adjust = c(0.01, 0.02, 0.03),
    stringsAsFactors = FALSE
  )

  expect_error(
    pathway_errorbar_table(
      abundance = abundance,
      daa_results_df = daa_results_df,
      Group = Group,
      p_values_threshold = 0.05
    ),
    "duplicated feature"
  )
})

test_that("pathway_errorbar_table rejects missing or incompatible Group labels", {
  abundance <- matrix(
    c(1, 1, 9, 9,
      9, 9, 1, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("f1", "f2"), paste0("S", 1:4))
  )
  daa_results_df <- data.frame(
    feature = c("f1", "f2"),
    method = "method",
    group1 = "A",
    group2 = "B",
    p_adjust = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  expect_error(
    pathway_errorbar_table(
      abundance = abundance,
      daa_results_df = daa_results_df,
      Group = c("A", NA, "B", "B"),
      p_values_threshold = 0.05
    ),
    "non-missing, non-empty group labels"
  )

  expect_error(
    pathway_errorbar_table(
      abundance = abundance,
      daa_results_df = daa_results_df,
      Group = c("X", "X", "Y", "Y"),
      p_values_threshold = 0.05
    ),
    "required DAA group"
  )
})

test_that("pathway_errorbar_table rejects significant DAA features missing from abundance", {
  abundance <- matrix(
    c(1, 1, 9, 9,
      9, 9, 1, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("f1", "f2"), paste0("S", 1:4))
  )
  daa_results_df <- data.frame(
    feature = c("f1", "missing_feature"),
    method = "method",
    group1 = "A",
    group2 = "B",
    p_adjust = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  expect_error(
    pathway_errorbar_table(
      abundance = abundance,
      daa_results_df = daa_results_df,
      Group = c("A", "A", "B", "B"),
      p_values_threshold = 0.05
    ),
    "missing from abundance row names: missing_feature"
  )
})
