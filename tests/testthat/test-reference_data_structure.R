# Tests for reference data structure

test_that("ko_to_kegg_reference has correct structure and coverage", {
  ko_to_kegg_reference <- ggpicrust2:::load_reference_data("ko_to_kegg")

  # Structure
  expect_s3_class(ko_to_kegg_reference, "data.frame")
  required_cols <- c("pathway_id", "ko_id", "pathway_name")
  expect_true(all(required_cols %in% names(ko_to_kegg_reference)))

  # Coverage (thresholds based on actual KEGG pathway maps,
  # excluding BRITE classification nodes)
  expect_gte(nrow(ko_to_kegg_reference), 50000)
  expect_gte(length(unique(ko_to_kegg_reference$pathway_id)), 500)
  expect_gte(length(unique(ko_to_kegg_reference$ko_id)), 20000)

  # No NA in key columns
  expect_equal(sum(is.na(ko_to_kegg_reference$pathway_id)), 0)
  expect_equal(sum(is.na(ko_to_kegg_reference$ko_id)), 0)

  # ID format (5-digit pathway IDs only, no BRITE classification nodes)
  expect_true(all(grepl("^ko\\d{5}$", ko_to_kegg_reference$pathway_id)))
  expect_true(all(grepl("^K\\d", trimws(ko_to_kegg_reference$ko_id))))
})

create_kegg_reference_test_data <- function(n_kos, n_samples, seed = 123) {
  set.seed(seed)
  ko_to_kegg_reference <- ggpicrust2:::load_reference_data("ko_to_kegg")
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), n_kos)

  data <- data.frame(function. = real_kos, stringsAsFactors = FALSE)
  for (i in seq_len(n_samples)) {
    data[[paste0("S", i)]] <- rpois(n_kos, lambda = 100)
  }
  data
}

test_that("ko2kegg_abundance works with real reference data", {
  test_data <- create_kegg_reference_test_data(100, 5)
  result <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_s3_class(result, "data.frame")
  expect_gt(nrow(result), 0)
  expect_equal(ncol(result), 5)
  expect_true(all(result >= 0))
  expect_false(any(is.na(result)))
  expect_equal(colnames(result), colnames(test_data)[-1])

  # Row names should be valid pathway IDs
  expect_true(all(grepl("^ko\\d{5}$", rownames(result))))
})

test_that("ko2kegg_abundance is deterministic", {
  test_data <- create_kegg_reference_test_data(50, 3, seed = 456)

  result1 <- suppressMessages(ko2kegg_abundance(data = test_data))
  result2 <- suppressMessages(ko2kegg_abundance(data = test_data))

  expect_identical(result1, result2)
})

test_that("ko2kegg_abundance fails fast on total KO/KEGG mismatch", {
  # Since 2.5.13 (commit 72eb140, "Fail fast on ko_to_kegg/pathway misuse"),
  # ko2kegg_abundance() errors out instead of silently returning an empty
  # matrix when no input KO maps into any KEGG pathway. Silent empties
  # cascade into cryptic downstream DAA failures ("No features in abundance
  # data") that are hard to diagnose, so the strict error is intentional.

  # Well-formed KO IDs that are not present in the KEGG reference.
  fake_data <- data.frame(function. = c("K99999", "K88888"),
                          S1 = c(100, 200),
                          stringsAsFactors = FALSE)
  expect_error(
    suppressMessages(ko2kegg_abundance(data = fake_data)),
    "No KO IDs in the input matched any KEGG pathway"
  )

  # Real KO IDs, but all-zero abundances collapse every pathway row to zero,
  # which the function treats as the same "nothing to analyze" condition.
  ko_to_kegg_reference <- ggpicrust2:::load_reference_data("ko_to_kegg")
  real_kos <- head(unique(ko_to_kegg_reference$ko_id), 10)
  zero_data <- data.frame(function. = real_kos,
                          S1 = rep(0, 10),
                          stringsAsFactors = FALSE)
  expect_error(
    suppressMessages(ko2kegg_abundance(data = zero_data)),
    "No KO IDs in the input matched any KEGG pathway"
  )
})

test_that("prepare_gene_sets creates valid gene sets from reference data", {
  gene_sets <- prepare_gene_sets("KEGG")

  expect_type(gene_sets, "list")
  expect_gt(length(gene_sets), 400)
  expect_false("ko01001" %in% names(gene_sets))
  expect_false("ko99980" %in% names(gene_sets))

  set_sizes <- sapply(gene_sets, length)
  expect_gt(min(set_sizes), 0)
  expect_lt(max(set_sizes), 5000)
})
