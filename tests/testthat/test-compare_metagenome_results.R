# Tests for compare_metagenome_results alignment semantics.
#
# Regression: the per-feature Spearman correlation used to index rows by
# position (metagenomes[[i]][k, ] vs metagenomes[[j]][k, ]). Two
# metagenomes with identical row names in different orders therefore
# compared *different* features and produced nonsense (often negative)
# correlations. We now align all metagenomes on the shared feature set by
# name before both the DAA cbind step and the correlation loop.

test_that("compare_metagenome_results aligns metagenomes by feature name, not by row position", {
  set.seed(1)
  n_feat <- 20
  n_samp <- 6

  m1 <- abs(matrix(rnorm(n_feat * n_samp), nrow = n_feat, ncol = n_samp))
  rownames(m1) <- paste0("K", sprintf("%04d", seq_len(n_feat)))
  colnames(m1) <- paste0("S", seq_len(n_samp))

  # m2 is m1 with rows reversed -- identical *content per feature name*,
  # but different row positions. A by-name aligned correlation must be
  # exactly 1.0 on the diagonal comparison.
  m2 <- m1[rev(rownames(m1)), , drop = FALSE]

  skip_if_not_installed("ComplexHeatmap")

  res <- suppressWarnings(suppressMessages(
    compare_metagenome_results(list(m1, m2),
                               names = c("m1", "m2"),
                               daa_method = "paired Wilcoxon",
                               correlation_permutations = 19)
  ))

  # Diagonal is trivially 1.
  expect_equal(res$correlation$cor_matrix["m1", "m1"], 1)
  expect_equal(res$correlation$cor_matrix["m2", "m2"], 1)
  # Off-diagonal should be 1 (identical features), not the ~ -0.5 the
  # by-position bug produced.
  expect_equal(res$correlation$cor_matrix["m1", "m2"], 1)
  expect_equal(res$correlation$cor_matrix["m2", "m1"], 1)
})

test_that("compare_metagenome_results errors cleanly when metagenomes have no row names", {
  set.seed(2)
  m1 <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  m2 <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  colnames(m1) <- paste0("S", 1:4)
  colnames(m2) <- paste0("S", 1:4)
  # deliberately no rownames on m2
  rownames(m1) <- paste0("K", 1:10)

  expect_error(
    compare_metagenome_results(list(m1, m2),
                               names = c("a", "b"),
                               daa_method = "paired Wilcoxon"),
    regexp = "row names"
  )
})

test_that("compare_metagenome_results errors when metagenomes share no features", {
  set.seed(3)
  m1 <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  m2 <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  rownames(m1) <- paste0("K", 1:10)
  rownames(m2) <- paste0("J", 1:10)  # disjoint feature namespace
  colnames(m1) <- paste0("S", 1:4)
  colnames(m2) <- paste0("S", 1:4)

  expect_error(
    compare_metagenome_results(list(m1, m2),
                               names = c("a", "b"),
                               daa_method = "paired Wilcoxon"),
    regexp = "No shared feature"
  )
})

test_that("compare_metagenome_results aligns metagenomes by sample name, not by column position", {
  # Regression: the per-feature Spearman correlation used to compute
  #   cor(m[i][k, ], m[j][k, ])
  # indexing columns by POSITION. Two metagenomes with identical
  # feature/sample content but columns reordered therefore correlated
  # "sample 1 of metagenome A" against "sample 10 of metagenome B" --
  # biologically meaningless -- and produced median correlations as low
  # as -0.6 for matrices that are really identical under by-name
  # alignment. Sample columns are now intersected by name before the
  # correlation loop so this is guaranteed to be correlation = 1.
  skip_if_not_installed("ComplexHeatmap")

  set.seed(4)
  n_feat <- 20
  n_samp <- 6
  m1 <- abs(matrix(rnorm(n_feat * n_samp), nrow = n_feat, ncol = n_samp))
  rownames(m1) <- paste0("K", sprintf("%04d", seq_len(n_feat)))
  colnames(m1) <- paste0("S", seq_len(n_samp))

  # m2 is m1 with columns in reverse order. Same biological content,
  # only the sample-axis positions differ.
  m2 <- m1[, rev(colnames(m1)), drop = FALSE]

  res <- suppressWarnings(suppressMessages(
    compare_metagenome_results(list(m1, m2),
                               names = c("m1", "m2"),
                               daa_method = "paired Wilcoxon",
                               correlation_permutations = 19)
  ))

  # By-name alignment must recover correlation = 1, not the ~-0.6 the
  # position-indexed bug produced on identical content.
  expect_equal(res$correlation$cor_matrix["m1", "m2"], 1)
  expect_equal(res$correlation$cor_matrix["m2", "m1"], 1)
})

test_that("compare_metagenome_results errors cleanly when metagenomes have no sample names", {
  # Parallel to the missing-rownames test: without colnames we cannot
  # align samples by name, and position alignment is the bug we are
  # preventing.
  set.seed(5)
  m1 <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  m2 <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  rownames(m1) <- paste0("K", 1:10)
  rownames(m2) <- paste0("K", 1:10)
  colnames(m1) <- paste0("S", 1:4)
  # deliberately no colnames on m2

  expect_error(
    compare_metagenome_results(list(m1, m2),
                               names = c("a", "b"),
                               daa_method = "paired Wilcoxon"),
    regexp = "column names"
  )
})

test_that("compare_metagenome_results validates metagenome labels and matrix keys before downstream analysis", {
  m1 <- matrix(
    seq_len(20),
    nrow = 5,
    dimnames = list(paste0("K", 1:5), paste0("S", 1:4))
  )
  m2 <- m1

  expect_error(
    compare_metagenome_results(list(m1), names = "m1"),
    "at least two"
  )

  expect_error(
    compare_metagenome_results(list(m1, m2), names = c("m", "m")),
    "duplicated labels"
  )

  duplicated_samples <- m2
  colnames(duplicated_samples)[2] <- colnames(duplicated_samples)[1]
  expect_error(
    compare_metagenome_results(list(m1, duplicated_samples),
                               names = c("m1", "m2")),
    "duplicated sample identifiers"
  )

  duplicated_features <- m2
  rownames(duplicated_features)[2] <- rownames(duplicated_features)[1]
  expect_error(
    compare_metagenome_results(list(m1, duplicated_features),
                               names = c("m1", "m2")),
    "duplicated feature identifiers"
  )

  negative_values <- m2
  negative_values[1, 1] <- -1
  expect_error(
    compare_metagenome_results(list(m1, negative_values),
                               names = c("m1", "m2")),
    "non-negative values"
  )
})

test_that("compare_metagenome_results rejects metagenomes whose sample sets do not overlap", {
  # Mismatched column counts used to fall through to stats::cor() and
  # abort mid-loop with "incompatible dimensions". Disjoint sample sets
  # now fail fast at the boundary with an actionable message.
  set.seed(6)
  m1 <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
  m2 <- abs(matrix(rnorm(60), nrow = 10, ncol = 6))
  rownames(m1) <- paste0("K", 1:10)
  rownames(m2) <- paste0("K", 1:10)
  colnames(m1) <- paste0("S",  1:4)
  colnames(m2) <- paste0("S2_", 1:6)  # disjoint sample namespace

  expect_error(
    compare_metagenome_results(list(m1, m2),
                               names = c("a", "b"),
                               daa_method = "paired Wilcoxon"),
    regexp = "No shared sample"
  )
})

test_that("compare_metagenome_results warns when sample intersection drops samples", {
  # A partial sample overlap is almost always a user mistake (parallel
  # quantifications should share the same biological samples). Proceed
  # on the shared set but make the drop visible.
  skip_if_not_installed("ComplexHeatmap")

  set.seed(7)
  # m1 has samples S1..S4, m2 has samples S1..S3 plus a unique S_extra.
  # Shared set is S1..S3; m1 drops S4, m2 drops S_extra.
  n_feat <- 15
  feat_names <- paste0("K", sprintf("%03d", seq_len(n_feat)))

  m1 <- abs(matrix(rnorm(n_feat * 4), nrow = n_feat, ncol = 4))
  rownames(m1) <- feat_names
  colnames(m1) <- c(paste0("S", 1:3), "S4")

  m2 <- abs(matrix(rnorm(n_feat * 4), nrow = n_feat, ncol = 4))
  rownames(m2) <- feat_names
  colnames(m2) <- c(paste0("S", 1:3), "S_extra")

  warnings <- character()
  withCallingHandlers(
    suppressMessages(
      compare_metagenome_results(list(m1, m2),
                                 names = c("m1", "m2"),
                                 daa_method = "paired Wilcoxon",
                                 correlation_permutations = 19)
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("Metagenome 'm1'", warnings)))
  expect_true(any(grepl("Metagenome 'm2'", warnings)))
})

test_that("compare_metagenome_results errors when no finite feature correlations are available", {
  skip_if_not_installed("ComplexHeatmap")

  m1 <- matrix(
    c(1, 1, 1, 1,
      2, 2, 2, 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("const1", "const2"), paste0("S", 1:4))
  )
  m2 <- m1

  expect_error(
    compare_metagenome_results(
      metagenomes = list(m1, m2),
      names = c("m1", "m2"),
      daa_method = "paired Wilcoxon",
      correlation_permutations = 19
    ),
    "No finite per-feature Spearman correlations"
  )
})

test_that("compare_metagenome_results rejects independent-group DAA methods", {
  m1 <- matrix(
    seq_len(30),
    nrow = 5,
    dimnames = list(paste0("K", 1:5), paste0("S", 1:6))
  )
  m2 <- m1 + 1

  expect_error(
    compare_metagenome_results(
      list(m1, m2),
      names = c("m1", "m2"),
      daa_method = "LinDA"
    ),
    "'daa_method' must be one of"
  )
})

test_that("paired Wilcoxon DAA retains sample pairing", {
  paired_test <- getFromNamespace(
    "run_paired_wilcoxon_comparison",
    "ggpicrust2"
  )
  feature1_group1 <- seq(0.10, 0.60, length.out = 6)
  feature1_group2 <- feature1_group1 + 0.05
  group1 <- rbind(
    feature1 = feature1_group1,
    feature2 = 1 - feature1_group1
  )
  group2 <- rbind(
    feature1 = feature1_group2,
    feature2 = 1 - feature1_group2
  )
  colnames(group1) <- colnames(group2) <- paste0("S", 1:6)

  result <- paired_test(
    group1,
    group2,
    group1 = "method1",
    group2 = "method2",
    p_adjust_method = "BH"
  )
  unpaired_p <- stats::wilcox.test(
    group2["feature1", ],
    group1["feature1", ],
    paired = FALSE,
    exact = FALSE
  )$p.value

  expect_true(all(result$method == "Paired Wilcoxon signed-rank test"))
  expect_true(all(result$group1 == "method1"))
  expect_true(all(result$group2 == "method2"))
  expect_lt(result$p_values[result$feature == "feature1"], unpaired_p)
  expect_gt(result$median_difference[result$feature == "feature1"], 0)
  expect_gt(result$log2_fold_change[result$feature == "feature1"], 0)
})

test_that("paired ALDEx2 DAA returns paired tests and effect direction", {
  skip_if_not_installed("ALDEx2")

  paired_aldex2 <- getFromNamespace(
    "run_paired_aldex2_comparison",
    "ggpicrust2"
  )
  set.seed(11)
  group1 <- matrix(
    sample(5:40, 24, replace = TRUE),
    nrow = 4,
    dimnames = list(paste0("K", 1:4), paste0("S", 1:6))
  )
  group2 <- group1
  group2["K1", ] <- group2["K1", ] + 10

  result <- suppressWarnings(suppressMessages(paired_aldex2(
    group1,
    group2,
    group1 = "method1",
    group2 = "method2",
    p_adjust_method = "BH"
  )))

  expect_setequal(
    unique(result$method),
    c("ALDEx2 paired t test",
      "ALDEx2 paired Wilcoxon signed-rank test")
  )
  expect_true(all(result$group1 == "method1"))
  expect_true(all(result$group2 == "method2"))
  expect_true(all(result$adj_method == "BH (ALDEx2 expected)"))
  expect_true(all(
    result$log2_fold_change[result$feature == "K1"] > 0
  ))
})

test_that("paired ALDEx2 aligns effect rows by feature identifier", {
  skip_if_not_installed("ALDEx2")

  paired_aldex2 <- getFromNamespace(
    "run_paired_aldex2_comparison",
    "ggpicrust2"
  )
  set.seed(12)
  feature_ids <- paste0("K", 1:4)
  group1 <- matrix(
    sample(5:40, 24, replace = TRUE),
    nrow = 4,
    dimnames = list(feature_ids, paste0("S", 1:6))
  )
  group2 <- group1 + 2

  local_mocked_bindings(
    aldex.effect = function(...) {
      output_ids <- rev(feature_ids)
      feature_number <- match(output_ids, feature_ids)
      data.frame(
        effect = feature_number,
        diff.btw = feature_number * 10,
        rab.all = feature_number * 100,
        overlap = feature_number / 10,
        row.names = output_ids
      )
    },
    .package = "ALDEx2"
  )

  result <- suppressWarnings(suppressMessages(paired_aldex2(
    group1,
    group2,
    group1 = "method1",
    group2 = "method2",
    p_adjust_method = "BH"
  )))

  first_test <- result[result$method == "ALDEx2 paired t test", ]
  expect_equal(first_test$feature, feature_ids)
  expect_equal(first_test$effect_size, seq_along(feature_ids))
  expect_equal(first_test$diff_btw, seq_along(feature_ids) * 10)
})

test_that("correlation permutation p-values preserve pairing and RNG state", {
  skip_if_not_installed("ComplexHeatmap")

  set.seed(101)
  m1 <- matrix(
    runif(20 * 6),
    nrow = 20,
    dimnames = list(paste0("K", 1:20), paste0("S", 1:6))
  )
  m2 <- m1

  set.seed(202)
  expected_next_random_value <- runif(1)
  set.seed(202)
  result <- suppressMessages(compare_metagenome_results(
    list(m1, m2),
    names = c("m1", "m2"),
    daa_method = "paired Wilcoxon",
    correlation_permutations = 39,
    correlation_seed = 7
  ))
  observed_next_random_value <- runif(1)

  expect_equal(observed_next_random_value, expected_next_random_value)
  expect_true(all(is.na(diag(result$correlation$p_matrix))))
  expect_true(all(is.na(diag(result$correlation$p_adjust_matrix))))
  expect_gte(result$correlation$p_matrix["m1", "m2"], 1 / 40)
  expect_lte(result$correlation$p_matrix["m1", "m2"], 1)
  expect_equal(
    result$correlation$p_adjust_matrix["m1", "m2"],
    result$correlation$p_matrix["m1", "m2"]
  )
  expect_equal(result$correlation$n_features_matrix["m1", "m2"], 20)
})

test_that("compare_metagenome_results validates correlation inference controls", {
  m1 <- matrix(
    seq_len(30),
    nrow = 5,
    dimnames = list(paste0("K", 1:5), paste0("S", 1:6))
  )
  m2 <- m1 + 1

  expect_error(
    compare_metagenome_results(
      list(m1, m2),
      names = c("m1", "m2"),
      daa_method = "paired Wilcoxon",
      correlation_permutations = -1
    ),
    "correlation_permutations.*non-negative"
  )
  expect_error(
    compare_metagenome_results(
      list(m1, m2),
      names = c("m1", "m2"),
      daa_method = "paired Wilcoxon",
      correlation_p_adjust_method = "invalid"
    ),
    "correlation_p_adjust_method.*must be one of"
  )

  two_sample <- m1[, 1:2, drop = FALSE]
  expect_error(
    compare_metagenome_results(
      list(two_sample, two_sample),
      names = c("m1", "m2"),
      daa_method = "paired Wilcoxon"
    ),
    "At least 3 shared samples"
  )
})
