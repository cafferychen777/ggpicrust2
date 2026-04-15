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

  skip_if_not_installed("MicrobiomeStat")
  skip_if_not_installed("ComplexHeatmap")

  res <- suppressWarnings(suppressMessages(
    compare_metagenome_results(list(m1, m2),
                               names = c("m1", "m2"),
                               daa_method = "LinDA")
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
                               daa_method = "LinDA"),
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
                               daa_method = "LinDA"),
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
  skip_if_not_installed("MicrobiomeStat")
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
                               daa_method = "LinDA")
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
                               daa_method = "LinDA"),
    regexp = "column names"
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
                               daa_method = "LinDA"),
    regexp = "No shared sample"
  )
})

test_that("compare_metagenome_results warns when sample intersection drops samples", {
  # A partial sample overlap is almost always a user mistake (parallel
  # quantifications should share the same biological samples). Proceed
  # on the shared set but make the drop visible.
  skip_if_not_installed("MicrobiomeStat")
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

  expect_warning(
    suppressMessages(
      compare_metagenome_results(list(m1, m2),
                                 names = c("m1", "m2"),
                                 daa_method = "LinDA")
    ),
    regexp = "sample"
  )
})
