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
