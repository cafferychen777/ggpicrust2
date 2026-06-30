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
