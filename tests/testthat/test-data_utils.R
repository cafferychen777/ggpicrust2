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

test_that("read_abundance_file reads gzipped delimited files", {
  raf <- getFromNamespace("read_abundance_file", "ggpicrust2")
  plain_file <- tempfile(fileext = ".tsv")
  gz_file <- tempfile(fileext = ".tsv.gz")
  on.exit(unlink(c(plain_file, gz_file)), add = TRUE)

  writeLines(c("feature\tS1\tS2", "K00001\t1\t2"), plain_file)
  con <- gzfile(gz_file, open = "wt")
  writeLines(c("feature\tS1\tS2", "K00001\t1\t2"), con)
  close(con)

  plain <- raf(plain_file)
  gzipped <- raf(gz_file)
  expect_equal(as.data.frame(gzipped), as.data.frame(plain))
})

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

test_that("normalize_abundance_feature_ids moves a leading ID column to row names", {
  normalize <- getFromNamespace("normalize_abundance_feature_ids", "ggpicrust2")

  df <- data.frame(
    feature = c("K00001", "K00002", "K00003"),
    S1 = c(10, 20, 30),
    S2 = c(5, 15, 25),
    stringsAsFactors = FALSE
  )

  normalized <- normalize(df, context = "unit-test abundance")

  expect_equal(rownames(normalized), df$feature)
  expect_equal(colnames(normalized), c("S1", "S2"))
  expect_true(all(vapply(normalized, is.numeric, logical(1))))
})

test_that("normalize_abundance_feature_ids rejects duplicated leading feature IDs", {
  normalize <- getFromNamespace("normalize_abundance_feature_ids", "ggpicrust2")

  df <- data.frame(
    feature = c("K00001", "K00001", "K00003"),
    S1 = c(10, 20, 30),
    S2 = c(5, 15, 25),
    stringsAsFactors = FALSE
  )

  expect_error(
    normalize(df, context = "unit-test abundance"),
    "duplicated identifiers: K00001"
  )
})

test_that("validate_abundance requires numeric sample values", {
  va <- getFromNamespace("validate_abundance", "ggpicrust2")

  bad_matrix <- matrix(
    as.character(1:6),
    nrow = 2,
    dimnames = list(c("f1", "f2"), c("S1", "S2", "S3"))
  )
  expect_error(
    va(bad_matrix),
    "matrix must be numeric"
  )

  bad_df <- data.frame(
    feature = c("f1", "f2"),
    S1 = c(1, 2),
    S2 = c("3", "4"),
    stringsAsFactors = FALSE
  )
  expect_error(
    va(bad_df),
    "Non-numeric sample column\\(s\\): S2"
  )
})

test_that("validate_abundance counts samples after excluding a leading ID column", {
  va <- getFromNamespace("validate_abundance", "ggpicrust2")

  one_sample <- data.frame(
    feature = c("f1", "f2"),
    S1 = c(1, 2),
    stringsAsFactors = FALSE
  )

  expect_error(
    va(one_sample, min_samples = 2),
    "At least 2 samples are required, found 1"
  )
})

test_that("validate_group requires a single metadata column name", {
  validate_group <- getFromNamespace("validate_group", "ggpicrust2")
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    batch = c("X", "X", "Y", "Y"),
    stringsAsFactors = FALSE
  )

  expect_true(validate_group(metadata, "group"))
  expect_error(
    validate_group(metadata, c("group", "batch")),
    "'group' must be a single non-empty character string"
  )
  expect_error(
    validate_group(metadata, NA_character_),
    "'group' must be a single non-empty character string"
  )
  expect_error(
    validate_group(metadata, ""),
    "'group' must be a single non-empty character string"
  )
})

test_that("DAA method name aliases canonicalize the legacy ALDEx2 spelling", {
  canonicalize <- getFromNamespace("canonicalize_daa_method_names", "ggpicrust2")
  validate_results <- getFromNamespace("validate_daa_results", "ggpicrust2")

  old_name <- "ALDEx2_Kruskal-Wallace test"
  new_name <- "ALDEx2_Kruskal-Wallis test"

  expect_equal(canonicalize(old_name), new_name)
  expect_equal(canonicalize(factor(old_name)), new_name)

  mixed_alias_df <- data.frame(
    feature = c("pathway1", "pathway2"),
    method = c(old_name, new_name),
    group1 = c("A", "A"),
    group2 = c("B", "B"),
    p_values = c(0.01, 0.02),
    p_adjust = c(0.03, 0.04),
    stringsAsFactors = FALSE
  )

  expect_true(validate_results(mixed_alias_df))
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

test_that("validate_count_parameter rejects invalid counts without coercion warnings", {
  validate_count <- getFromNamespace("validate_count_parameter", "ggpicrust2")

  expect_error(
    expect_warning(validate_count(1e20, "n_pathways"), NA),
    "n_pathways.*single finite integer"
  )
  expect_error(
    validate_count(1.5, "n_pathways"),
    "n_pathways.*single finite integer"
  )
  expect_error(
    validate_count(0, "n_pathways"),
    "n_pathways.*positive"
  )
  expect_true(validate_count(0, "n_pathways", allow_zero = TRUE))
})

test_that("align_samples is idempotent on already-aligned inputs", {
  # Contract: ggpicrust2() pre-aligns inputs at its top, then calls
  # pathway_daa(), which independently re-aligns. Both sites use the
  # same helper with identical arguments, so the duplicate is bounded
  # by the invariant that running align_samples() twice on the result
  # of the first call is a no-op. This test locks that invariant --
  # any future change to align_samples() that breaks idempotency
  # would silently reintroduce drift between the two paths.
  set.seed(11)
  abund <- matrix(runif(20), nrow = 5, ncol = 4,
                  dimnames = list(paste0("K", 1:5),
                                  paste0("S", 1:4)))
  meta <- data.frame(
    sample = paste0("S", 1:4),
    Env    = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  first  <- align_samples(abund, meta, verbose = FALSE)
  second <- align_samples(first$abundance, first$metadata, verbose = FALSE)

  expect_identical(first$abundance, second$abundance)
  expect_identical(first$metadata,  second$metadata)
  expect_identical(first$sample_col, second$sample_col)
  expect_identical(first$n_samples,  second$n_samples)
})
