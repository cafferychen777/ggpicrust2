# Tests for find_sample_column / samples_match heuristics in data_utils.R

test_that("find_sample_column uses standard names when available", {
  fsc <- getFromNamespace("find_sample_column", "ggpicrust2")
  meta <- data.frame(sample = paste0("S", 1:4),
                     group = c("a", "a", "b", "b"),
                     stringsAsFactors = FALSE)
  expect_equal(fsc(meta, paste0("S", 1:4)), "sample")
})

test_that("find_sample_column rejects a standard-named column with duplicate values", {
  # Regression: Priority 1 used to pick a standard-named column (e.g.
  # "samples") even when its values were duplicated, simply because the
  # column name was on the standard list. That silently misaligned the
  # downstream analysis. Require uniqueness here too, falling through
  # to Priority 2 when the standard-named column is not a valid sample id.
  fsc <- getFromNamespace("find_sample_column", "ggpicrust2")

  meta <- data.frame(
    sample  = c("S1", "S1", "S2", "S2"),   # duplicates -> not a valid id
    real_id = c("S1", "S2", "S3", "S4"),   # correct sample id
    stringsAsFactors = FALSE
  )
  expect_equal(fsc(meta, paste0("S", 1:4)), "real_id")
})

test_that("find_sample_column rejects categorical columns at Priority 2", {
  # Regression: the scan-every-column fallback used a 50% overlap threshold
  # with no uniqueness check, so a categorical column whose levels happened
  # to share a few strings with the sample names would be mistakenly picked
  # as the sample identifier. Require uniqueness + >= 90% overlap.
  fsc <- getFromNamespace("find_sample_column", "ggpicrust2")

  meta <- data.frame(
    subject_label = c("S1", "S1", "S2", "S2"),  # duplicates -> not a sample id
    my_samples    = c("S1", "S2", "S3", "S4"),  # unique, correct
    stringsAsFactors = FALSE
  )
  # `my_samples` is the only real sample column; `subject_label` must be
  # rejected despite sharing values with abundance_samples.
  expect_equal(fsc(meta, paste0("S", 1:4)), "my_samples")
})

test_that("find_sample_column rejects low-overlap coincidental columns", {
  fsc <- getFromNamespace("find_sample_column", "ggpicrust2")

  meta <- data.frame(
    note = c("S1", "other1", "other2", "other3"),  # only 25% overlap
    real_sample = c("S1", "S2", "S3", "S4"),
    stringsAsFactors = FALSE
  )
  expect_equal(fsc(meta, paste0("S", 1:4)), "real_sample")
})

test_that("find_sample_column uses rownames only with near-complete overlap", {
  fsc <- getFromNamespace("find_sample_column", "ggpicrust2")

  meta <- data.frame(group = c("a", "a", "b", "b"),
                     stringsAsFactors = FALSE)
  rownames(meta) <- paste0("S", 1:4)
  expect_equal(fsc(meta, paste0("S", 1:4)), ".rownames")

  # Default integer rownames should not match arbitrary sample names even
  # though intersect() might find a stray match.
  meta2 <- data.frame(group = c("a", "a", "b", "b"),
                      stringsAsFactors = FALSE)
  expect_null(fsc(meta2, paste0("S", 1:4)))
})
