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
