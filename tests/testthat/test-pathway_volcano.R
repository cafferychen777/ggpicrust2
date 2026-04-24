# Helper: create standard volcano test data
create_volcano_test_data <- function() {
  set.seed(123)
  data.frame(
    feature = paste0("ko", sprintf("%05d", 1:20)),
    pathway_name = paste0("Pathway ", 1:20),
    log2FoldChange = rnorm(20, 0, 2),
    p_adjust = runif(20, 0, 0.1),
    stringsAsFactors = FALSE
  )
}

test_that("pathway_volcano creates a ggplot object", {
  skip_if_not_installed("ggrepel")

  p <- pathway_volcano(create_volcano_test_data())
  expect_s3_class(p, "ggplot")
})

test_that("pathway_volcano handles custom column names", {
  skip_if_not_installed("ggrepel")

  set.seed(123)
  daa_results <- data.frame(
    id = paste0("ko", sprintf("%05d", 1:20)),
    name = paste0("Pathway ", 1:20),
    fc = rnorm(20, 0, 2),
    pval = runif(20, 0, 0.1),
    stringsAsFactors = FALSE
  )

  p <- pathway_volcano(daa_results, fc_col = "fc", p_col = "pval", label_col = "name")
  expect_s3_class(p, "ggplot")
})

test_that("pathway_volcano handles custom thresholds", {
  skip_if_not_installed("ggrepel")

  p <- pathway_volcano(
    create_volcano_test_data(),
    fc_threshold = 0.5,
    p_threshold = 0.01,
    label_top_n = 5
  )
  expect_s3_class(p, "ggplot")
})

test_that("pathway_volcano handles custom colors", {
  skip_if_not_installed("ggrepel")

  p <- pathway_volcano(
    create_volcano_test_data(),
    colors = c("Down" = "blue", "Not Significant" = "gray", "Up" = "red")
  )
  expect_s3_class(p, "ggplot")
})

test_that("pathway_volcano handles no labels", {
  skip_if_not_installed("ggrepel")

  p <- pathway_volcano(create_volcano_test_data(), label_top_n = 0)
  expect_s3_class(p, "ggplot")
})

test_that("pathway_volcano errors on missing columns", {
  set.seed(123)
  daa_results <- data.frame(
    feature = paste0("ko", sprintf("%05d", 1:20)),
    some_value = rnorm(20),
    stringsAsFactors = FALSE
  )

  expect_error(pathway_volcano(daa_results), "not found")
})

test_that("pathway_volcano errors on non-data.frame input", {
  expect_error(pathway_volcano("not a data frame"), "must be a data frame")
})

test_that("pathway_volcano keeps all-zero p-values finite", {
  daa_results <- data.frame(
    feature = c("ko00001", "ko00002"),
    pathway_name = c("Pathway 1", "Pathway 2"),
    log2_fold_change = c(2, -2),
    p_adjust = c(0, 0),
    stringsAsFactors = FALSE
  )

  expect_warning(
    p <- pathway_volcano(daa_results),
    regexp = NA
  )
  expect_true(all(is.finite(p$data$neg_log10_p)))
  expect_true(all(p$data$neg_log10_p > 0))
})

test_that("pathway_volcano rejects negative p-values", {
  daa_results <- data.frame(
    feature = "ko00001",
    pathway_name = "Pathway 1",
    log2_fold_change = 2,
    p_adjust = -0.01,
    stringsAsFactors = FALSE
  )

  expect_error(pathway_volcano(daa_results), "negative")
})

test_that("pathway_volcano works with real DAA workflow", {
  skip_if_not_installed("ggrepel")
  skip_on_cran()

  data(ko_abundance, package = "ggpicrust2")
  data(metadata, package = "ggpicrust2")

  kegg_abundance <- ko2kegg_abundance(data = ko_abundance)

  daa_results <- suppressMessages(pathway_daa(
    abundance = kegg_abundance,
    metadata = metadata,
    group = "Environment",
    daa_method = "LinDA"
  ))

  daa_annotated <- suppressMessages(pathway_annotation(
    pathway = "KO",
    ko_to_kegg = TRUE,
    daa_results_df = daa_results
  ))

  p <- pathway_volcano(daa_annotated, label_top_n = 5)
  expect_s3_class(p, "ggplot")
})
