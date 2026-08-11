test_that("compare_gsea_daa aligns DAA log2 fold changes to GSEA direction", {
  gsea_results <- data.frame(
    pathway_id = c("ko00010", "ko00020"),
    NES = c(1.5, -1.2),
    p.adjust = c(0.01, 0.02),
    group1 = c("Treatment", "Treatment"),
    group2 = c("Control", "Control"),
    stringsAsFactors = FALSE
  )
  daa_results <- data.frame(
    feature = c("ko00010", "ko00020"),
    log2_fold_change = c(-2, 1),
    p_adjust = c(0.01, 0.02),
    group1 = c("Treatment", "Treatment"),
    group2 = c("Control", "Control"),
    stringsAsFactors = FALSE
  )

  result <- compare_gsea_daa(
    gsea_results = gsea_results,
    daa_results = daa_results,
    plot_type = "scatter"
  )

  scatter_data <- result$results$scatter_data
  expect_equal(scatter_data$log2_fold_change, c(-2, 1))
  expect_equal(scatter_data$daa_log2_fold_change_aligned, c(2, -1))
  expect_equal(result$plot$labels$y, "Log2 Fold Change (DAA, aligned to GSEA direction)")
})

test_that("compare_gsea_daa keeps DAA direction when group2 matches GSEA-positive group", {
  gsea_results <- data.frame(
    pathway_id = "ko00010",
    NES = 1.5,
    p.adjust = 0.01,
    group1 = "Treatment",
    group2 = "Control",
    stringsAsFactors = FALSE
  )
  daa_results <- data.frame(
    feature = "ko00010",
    log2_fold_change = 2,
    p_adjust = 0.01,
    group1 = "Control",
    group2 = "Treatment",
    stringsAsFactors = FALSE
  )

  result <- compare_gsea_daa(
    gsea_results = gsea_results,
    daa_results = daa_results,
    plot_type = "scatter"
  )

  expect_equal(result$results$scatter_data$daa_log2_fold_change_aligned, 2)
})

test_that("compare_gsea_daa scatter data preserves GSEA input order", {
  gsea_results <- data.frame(
    pathway_id = c("z_pathway", "a_pathway"),
    NES = c(1.5, -1.2),
    p.adjust = c(0.01, 0.02),
    group1 = "Treatment",
    group2 = "Control",
    stringsAsFactors = FALSE
  )
  daa_results <- data.frame(
    feature = c("a_pathway", "z_pathway"),
    log2_fold_change = c(1, -2),
    p_adjust = c(0.02, 0.01),
    group1 = "Treatment",
    group2 = "Control",
    stringsAsFactors = FALSE
  )

  result <- compare_gsea_daa(gsea_results, daa_results, plot_type = "scatter")

  expect_equal(result$results$scatter_data$pathway_id,
               gsea_results$pathway_id)
  expect_equal(result$results$scatter_data$log2_fold_change, c(-2, 1))
})

test_that("compare_gsea_daa requires explicit directions for scatter effect sizes", {
  gsea_results <- data.frame(
    pathway_id = "ko00010",
    NES = 1.5,
    p.adjust = 0.01,
    stringsAsFactors = FALSE
  )
  daa_results <- data.frame(
    feature = "ko00010",
    log2_fold_change = 2,
    p_adjust = 0.01,
    group1 = "Control",
    group2 = "Treatment",
    stringsAsFactors = FALSE
  )

  expect_error(
    compare_gsea_daa(
      gsea_results = gsea_results,
      daa_results = daa_results,
      plot_type = "scatter"
    ),
    "requires explicit direction columns"
  )
})

test_that("compare_gsea_daa validates plot_type as a single supported choice", {
  gsea_results <- data.frame(
    pathway_id = "ko00010",
    p.adjust = 0.01,
    stringsAsFactors = FALSE
  )
  daa_results <- data.frame(
    feature = "ko00010",
    p_adjust = 0.01,
    stringsAsFactors = FALSE
  )

  expect_error(
    compare_gsea_daa(gsea_results, daa_results, plot_type = c("venn", "scatter")),
    "'plot_type' must be one of"
  )
  expect_error(
    compare_gsea_daa(gsea_results, daa_results, plot_type = NA_character_),
    "'plot_type' must be one of"
  )
})

test_that("compare_gsea_daa handles an empty significant UpSet universe", {
  gsea_results <- data.frame(
    pathway_id = c("ko00010", "ko00020"),
    p.adjust = c(0.5, 0.6),
    stringsAsFactors = FALSE
  )
  daa_results <- data.frame(
    feature = c("ko00010", "ko00030"),
    p_adjust = c(0.7, 0.8),
    stringsAsFactors = FALSE
  )

  result <- compare_gsea_daa(
    gsea_results,
    daa_results,
    plot_type = "upset"
  )

  expect_s3_class(result$plot, "ggplot")
  expect_equal(result$results$n_gsea_total, 0)
  expect_equal(result$results$n_daa_total, 0)
  expect_equal(result$plot$data$count, c(0, 0, 0))
})

test_that("compare_gsea_daa handles an empty significant Venn universe", {
  gsea_results <- data.frame(
    pathway_id = c("ko00010", "ko00020"),
    p.adjust = c(0.5, 0.6),
    stringsAsFactors = FALSE
  )
  daa_results <- data.frame(
    feature = c("ko00010", "ko00030"),
    p_adjust = c(0.7, 0.8),
    stringsAsFactors = FALSE
  )

  result <- compare_gsea_daa(
    gsea_results,
    daa_results,
    plot_type = "venn"
  )

  expect_s3_class(result$plot, "ggplot")
  expect_equal(result$plot$data$count, c(0, 0, 0))
})

test_that("compare_gsea_daa scatter rejects missing p-values", {
  gsea_results <- data.frame(
    pathway_id = "ko00010",
    NES = 1.5,
    p.adjust = NA_real_,
    group1 = "Treatment",
    group2 = "Control",
    stringsAsFactors = FALSE
  )
  daa_results <- data.frame(
    feature = "ko00010",
    log2_fold_change = 2,
    p_adjust = 0.01,
    group1 = "Control",
    group2 = "Treatment",
    stringsAsFactors = FALSE
  )

  expect_error(
    compare_gsea_daa(gsea_results, daa_results, plot_type = "scatter"),
    "scatter plot.*between 0 and 1"
  )
})

test_that("compare_gsea_daa preserves subnormal adjusted p-values", {
  gsea_results <- data.frame(
    pathway_id = "ko00010",
    NES = 1.5,
    p.adjust = 1e-320,
    group1 = "Treatment",
    group2 = "Control",
    stringsAsFactors = FALSE
  )
  daa_results <- data.frame(
    feature = "ko00010",
    log2_fold_change = 2,
    p_adjust = 0.01,
    group1 = "Control",
    group2 = "Treatment",
    stringsAsFactors = FALSE
  )

  result <- compare_gsea_daa(
    gsea_results,
    daa_results,
    plot_type = "scatter"
  )

  expect_equal(
    result$results$scatter_data$gsea_neg_log10_p_adjust,
    320,
    tolerance = 1e-6
  )
})

test_that("compare_gsea_daa set plots require comparable group pairs when available", {
  gsea_results <- data.frame(
    pathway_id = "ko00010",
    p.adjust = 0.01,
    group1 = "Treatment",
    group2 = "Control",
    stringsAsFactors = FALSE
  )
  daa_results <- data.frame(
    feature = "ko00010",
    p_adjust = 0.01,
    group1 = "Control",
    group2 = "Treatment",
    stringsAsFactors = FALSE
  )

  result <- compare_gsea_daa(
    gsea_results = gsea_results,
    daa_results = daa_results,
    plot_type = "venn"
  )
  expect_equal(result$results$n_overlap, 1)

  gsea_multi_pair <- rbind(
    gsea_results,
    data.frame(
      pathway_id = "ko00020",
      p.adjust = 0.02,
      group1 = "Treatment",
      group2 = "Placebo",
      stringsAsFactors = FALSE
    )
  )
  expect_error(
    compare_gsea_daa(
      gsea_results = gsea_multi_pair,
      daa_results = daa_results,
      plot_type = "venn"
    ),
    "requires one comparable group pair"
  )

  daa_incompatible <- daa_results
  daa_incompatible$group1 <- "Placebo"
  expect_error(
    compare_gsea_daa(
      gsea_results = gsea_results,
      daa_results = daa_incompatible,
      plot_type = "upset"
    ),
    "not comparable"
  )
})

test_that("compare_gsea_daa rejects partial direction schemas", {
  gsea_results <- data.frame(
    pathway_id = "ko00010",
    p.adjust = 0.01,
    group1 = "Treatment",
    stringsAsFactors = FALSE
  )
  daa_results <- data.frame(
    feature = "ko00010",
    p_adjust = 0.01,
    stringsAsFactors = FALSE
  )

  expect_error(
    compare_gsea_daa(gsea_results, daa_results, plot_type = "venn"),
    "must provide both 'group1' and 'group2'"
  )

  gsea_results$group2 <- "Control"
  daa_results$group2 <- "Treatment"
  expect_error(
    compare_gsea_daa(gsea_results, daa_results, plot_type = "upset"),
    "daa_results contains only 'group2'"
  )
})

test_that("compare_gsea_daa rejects incompatible scatter directions", {
  gsea_results <- data.frame(
    pathway_id = "ko00010",
    NES = 1.5,
    p.adjust = 0.01,
    group1 = "Treatment",
    group2 = "Control",
    stringsAsFactors = FALSE
  )
  daa_results <- data.frame(
    feature = "ko00010",
    log2_fold_change = 2,
    p_adjust = 0.01,
    group1 = "Control",
    group2 = "Placebo",
    stringsAsFactors = FALSE
  )

  expect_error(
    compare_gsea_daa(
      gsea_results = gsea_results,
      daa_results = daa_results,
      plot_type = "scatter"
    ),
    "not comparable"
  )
})
