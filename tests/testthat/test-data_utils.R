# Tests for abundance-quality helpers in data_utils.R:
#   - compute_relative_abundance()
#   - validate_abundance(check_zero_columns = TRUE)
#
# Regression: x / sum(x) used to produce NaN on zero-sum sample columns,
# which then got silently dropped by downstream mean(..., na.rm = TRUE)
# aggregations. That computed group stats from the wrong sample size
# without any warning. Both the helper and the input validator now
# refuse zero-sum columns with an actionable error that names the
# offending sample(s).

test_that("compute_relative_abundance matches the old apply idiom on clean input", {
  cra <- getFromNamespace("compute_relative_abundance", "ggpicrust2")

  set.seed(1)
  M <- matrix(sample(1:20, 12), nrow = 3, ncol = 4,
              dimnames = list(paste0("f", 1:3), paste0("S", 1:4)))

  old_way <- apply(t(M), 1, function(x) x / sum(x))
  new_way <- cra(M)
  expect_equal(new_way, old_way)
  # Columns sum to 1 after normalization.
  expect_equal(unname(colSums(new_way)), rep(1, ncol(M)))
})

test_that("compute_relative_abundance rejects zero-sum columns with a named error", {
  cra <- getFromNamespace("compute_relative_abundance", "ggpicrust2")

  M <- matrix(c(10, 20, 0, 40,
                15, 25, 0, 45,
                12, 22, 0, 42),
              nrow = 3, byrow = TRUE,
              dimnames = list(paste0("f", 1:3), c("S1", "S2", "BAD", "S4")))

  expect_error(cra(M, context = "unit-test"),
               regexp = "BAD",
               class = "simpleError")
  expect_error(cra(M), regexp = "unit-test|abundance", fixed = FALSE)
})

test_that("compute_relative_abundance rejects columns that sum to NA", {
  cra <- getFromNamespace("compute_relative_abundance", "ggpicrust2")

  M <- matrix(c(10, NA, 30,
                15, NA, 35),
              nrow = 2, byrow = TRUE,
              dimnames = list(paste0("f", 1:2), c("S1", "S2", "S3")))

  expect_error(cra(M), regexp = "S2")
})

test_that("compute_relative_abundance falls back to index labels when colnames are missing", {
  cra <- getFromNamespace("compute_relative_abundance", "ggpicrust2")

  M <- matrix(c(1, 0, 3,
                2, 0, 4),
              nrow = 2, byrow = TRUE)  # no colnames

  # Position 2 is zero-sum; error must still identify it.
  expect_error(cra(M), regexp = "\\b2\\b")
})

test_that("validate_abundance rejects zero-sum sample columns by default", {
  va <- getFromNamespace("validate_abundance", "ggpicrust2")

  M <- matrix(c(10, 20, 0, 40,
                15, 25, 0, 45),
              nrow = 2, byrow = TRUE,
              dimnames = list(paste0("f", 1:2), c("S1", "S2", "BAD", "S4")))

  expect_error(va(M), regexp = "BAD")

  # Data-frame form (abundance-only columns) also checked.
  df <- as.data.frame(M)
  expect_error(va(df), regexp = "BAD")
})

test_that("validate_abundance can be asked to skip the zero-column check", {
  va <- getFromNamespace("validate_abundance", "ggpicrust2")

  M <- matrix(c(10, 0,
                15, 0),
              nrow = 2, byrow = TRUE,
              dimnames = list(paste0("f", 1:2), c("S1", "BAD")))

  expect_true(va(M, check_zero_columns = FALSE))
})

test_that("validate_abundance tolerates a non-numeric ID column in data frames", {
  # abundance_to_numeric_matrix() must strip a character ID column before
  # summing, otherwise data.frame inputs from earlier pipeline stages
  # would trigger spurious "not numeric" errors.
  va <- getFromNamespace("validate_abundance", "ggpicrust2")

  df <- data.frame(
    feature = paste0("K", 1:3),
    S1 = c(10, 20, 30),
    S2 = c(5, 15, 25),
    stringsAsFactors = FALSE
  )
  expect_true(va(df))
})

# Unit tests for summarize_abundance_by_group(), the per-feature/per-group
# mean/sd helper that pathway_errorbar() and calculate_abundance_stats()
# both route through. Locks down contract so future refactors do not
# silently break one caller while leaving the other intact.
test_that("summarize_abundance_by_group matches manual mean/sd per group", {
  sbg <- getFromNamespace("summarize_abundance_by_group", "ggpicrust2")

  set.seed(42)
  M <- matrix(runif(12), nrow = 3, ncol = 4,
              dimnames = list(c("f1", "f2", "f3"),
                              c("S1", "S2", "S3", "S4")))
  g <- c("A", "A", "B", "B")

  out <- sbg(M, g)
  # Feature-major order: (f1,A), (f1,B), (f2,A), ... — required for
  # downstream compatibility with pathway_errorbar's match()/order() code.
  expect_equal(out$name, rep(c("f1", "f2", "f3"), each = 2))
  expect_equal(out$group, rep(c("A", "B"), times = 3))

  # Numerics must match manual computation exactly.
  for (feat in rownames(M)) {
    manual_A_mean <- mean(M[feat, g == "A"])
    manual_A_sd   <- stats::sd(M[feat, g == "A"])
    manual_B_mean <- mean(M[feat, g == "B"])
    manual_B_sd   <- stats::sd(M[feat, g == "B"])
    row <- out[out$name == feat, ]
    expect_equal(row$mean[row$group == "A"], manual_A_mean)
    expect_equal(row$sd[row$group == "A"], manual_A_sd)
    expect_equal(row$mean[row$group == "B"], manual_B_mean)
    expect_equal(row$sd[row$group == "B"], manual_B_sd)
  }
})

test_that("summarize_abundance_by_group handles NA in abundance via na.rm = TRUE", {
  sbg <- getFromNamespace("summarize_abundance_by_group", "ggpicrust2")

  M <- matrix(c(1, 2, NA, 4,
                5, 6, 7, 8),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("f1", "f2"), c("S1", "S2", "S3", "S4")))
  g <- c("A", "A", "B", "B")
  out <- sbg(M, g)

  # f1 in group B: only S4 is non-NA, so mean = 4, sd = NA (stats::sd on
  # length-1 vector returns NA). This matches calculate_abundance_stats().
  f1_B <- out[out$name == "f1" & out$group == "B", ]
  expect_equal(f1_B$mean, 4)
  expect_true(is.na(f1_B$sd))

  # f2 across all samples is intact.
  f2_A <- out[out$name == "f2" & out$group == "A", ]
  expect_equal(f2_A$mean, mean(c(5, 6)))
  expect_equal(f2_A$sd, stats::sd(c(5, 6)))
})

test_that("summarize_abundance_by_group errors on length mismatch", {
  sbg <- getFromNamespace("summarize_abundance_by_group", "ggpicrust2")

  M <- matrix(1:6, nrow = 2, ncol = 3,
              dimnames = list(c("f1", "f2"), c("S1", "S2", "S3")))
  expect_error(sbg(M, c("A", "B")), regexp = "length")
})
