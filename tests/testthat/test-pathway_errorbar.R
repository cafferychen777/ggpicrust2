# Helper: create standard errorbar test data
create_errorbar_test_data <- function(n_features = 5, p_adjust = NULL) {
  set.seed(123)
  n_samples <- 10

  abundance <- matrix(runif(n_features * n_samples), nrow = n_features, ncol = n_samples)
  rownames(abundance) <- paste0("pathway", 1:n_features)
  colnames(abundance) <- paste0("sample", 1:n_samples)

  if (is.null(p_adjust)) p_adjust <- rep(0.01, n_features)

  daa_results_df <- data.frame(
    feature = paste0("pathway", 1:n_features),
    pathway_name = paste0("Pathway ", 1:n_features),
    p_adjust = p_adjust,
    method = rep("ALDEx2_Welch's t test", n_features),
    group1 = rep("GroupA", n_features),
    group2 = rep("GroupB", n_features),
    stringsAsFactors = FALSE
  )

  Group <- factor(rep(c("GroupA", "GroupB"), each = n_samples / 2))
  names(Group) <- paste0("sample", 1:n_samples)

  list(abundance = abundance, daa_results_df = daa_results_df, Group = Group)
}

test_that("pathway_errorbar basic functionality works", {
  td <- create_errorbar_test_data(
    n_features = 10,
    p_adjust = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)
  )

  p <- pathway_errorbar(
    abundance = td$abundance,
    daa_results_df = td$daa_results_df,
    Group = td$Group,
    p_values_threshold = 0.05,
    select = paste0("pathway", 1:5),
    x_lab = "pathway_name"
  )

  expect_s3_class(p, "patchwork")
})

test_that("pathway_errorbar pathway_names_text_size parameter works", {
  td <- create_errorbar_test_data(
    n_features = 10,
    p_adjust = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)
  )

  # Test with auto text size
  p1 <- pathway_errorbar(
    abundance = td$abundance,
    daa_results_df = td$daa_results_df,
    Group = td$Group,
    p_values_threshold = 0.05,
    select = paste0("pathway", 1:5),
    x_lab = "pathway_name",
    pathway_names_text_size = "auto"
  )

  expect_s3_class(p1, "patchwork")

  # Test with custom text size
  p2 <- pathway_errorbar(
    abundance = td$abundance,
    daa_results_df = td$daa_results_df,
    Group = td$Group,
    p_values_threshold = 0.05,
    select = paste0("pathway", 1:5),
    x_lab = "pathway_name",
    pathway_names_text_size = 12
  )

  expect_s3_class(p2, "patchwork")
})

test_that("pathway_errorbar handles missing annotations", {
  td <- create_errorbar_test_data(n_features = 2, p_adjust = c(0.01, 0.02))
  td$daa_results_df$pathway_name <- c(NA, "Pathway 2")

  p <- pathway_errorbar(
    abundance = td$abundance,
    daa_results_df = td$daa_results_df,
    Group = td$Group,
    x_lab = "pathway_name"
  )
  expect_s3_class(p, "patchwork")
})

test_that("pathway_errorbar handles too many features", {
  td <- create_errorbar_test_data(n_features = 31)

  expect_warning(
    pathway_errorbar(
      abundance = td$abundance,
      daa_results_df = td$daa_results_df,
      Group = td$Group,
      x_lab = "pathway_name"
    ),
    regexp = "Found \\d+ significant features"
  )
})

test_that("pathway_errorbar handles custom colors correctly", {
  td <- create_errorbar_test_data()

  p <- pathway_errorbar(
    abundance = td$abundance,
    daa_results_df = td$daa_results_df,
    Group = td$Group,
    colors = c("#FF0000", "#0000FF"),
    x_lab = "pathway_name"
  )

  expect_s3_class(p, "patchwork")
})

test_that("pathway_errorbar handles different ordering options", {
  td <- create_errorbar_test_data(p_adjust = c(0.04, 0.01, 0.03, 0.02, 0.05))
  td$daa_results_df$pathway_class <- c("Class1", "Class1", "Class2", "Class2", "Class3")

  for (order_type in c("p_values", "name", "group", "pathway_class")) {
    p <- pathway_errorbar(
      abundance = td$abundance,
      daa_results_df = td$daa_results_df,
      Group = td$Group,
      order = order_type,
      x_lab = "pathway_name"
    )
    expect_s3_class(p, "patchwork")
  }

  # Invalid order type (function lacks upfront validation; crashes downstream)
  expect_error(
    pathway_errorbar(
      abundance = td$abundance,
      daa_results_df = td$daa_results_df,
      Group = td$Group,
      order = "invalid_order",
      x_lab = "pathway_name"
    )
  )
})

test_that("pathway_errorbar regression: ko_to_kegg TRUE with pathway_class order", {
  td <- create_errorbar_test_data(
    n_features = 6,
    p_adjust = c(0.001, 0.002, 0.003, 0.004, 0.005, 0.006)
  )
  td$daa_results_df$pathway_class <- c("Class2", "Class1", "Class2", "Class1", "Class3", "Class3")

  # Regression guard for Issue #177/#196 path:
  # ko_to_kegg=TRUE + order='pathway_class' used to fail with
  # `tibble::column_to_rownames()` / "Can't find column `.`".
  expect_error(
    pathway_errorbar(
      abundance = td$abundance,
      daa_results_df = td$daa_results_df,
      Group = td$Group,
      ko_to_kegg = TRUE,
      order = "pathway_class",
      p_values_threshold = 0.05,
      x_lab = "pathway_name"
    ),
    NA
  )
})

test_that("pathway_errorbar aligns Group by names when provided", {
  td <- create_errorbar_test_data(
    n_features = 4,
    p_adjust = c(0.001, 0.002, 0.003, 0.004)
  )
  td$daa_results_df$pathway_class <- c("Class1", "Class1", "Class2", "Class2")

  # Intentionally shuffle Group order but keep sample names.
  shuffled_group <- td$Group[c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1)]

  expect_error(
    pathway_errorbar(
      abundance = td$abundance,
      daa_results_df = td$daa_results_df,
      Group = shuffled_group,
      ko_to_kegg = TRUE,
      order = "pathway_class",
      p_values_threshold = 0.05,
      x_lab = "pathway_name"
    ),
    NA
  )
})

test_that("pathway_errorbar handles p_value_bar parameter correctly", {
  td <- create_errorbar_test_data()

  p1 <- pathway_errorbar(
    abundance = td$abundance,
    daa_results_df = td$daa_results_df,
    Group = td$Group,
    p_value_bar = TRUE,
    x_lab = "pathway_name"
  )

  p2 <- pathway_errorbar(
    abundance = td$abundance,
    daa_results_df = td$daa_results_df,
    Group = td$Group,
    p_value_bar = FALSE,
    x_lab = "pathway_name"
  )

  expect_s3_class(p1, "patchwork")
  expect_s3_class(p2, "patchwork")
  expect_false(identical(p1, p2))
})

test_that("pathway_errorbar_table function works correctly", {
  td <- create_errorbar_test_data(n_features = 3, p_adjust = c(0.01, 0.02, 0.03))

  metadata <- data.frame(
    sample = colnames(td$abundance),
    group = td$Group,
    stringsAsFactors = FALSE
  )
  daa_results <- pathway_daa(
    abundance = td$abundance,
    metadata = metadata,
    group = "group",
    daa_method = "ALDEx2"
  )
  daa_single_method <- daa_results[daa_results$method == "ALDEx2_Welch's t test", ]

  result <- pathway_errorbar_table(
    abundance = td$abundance,
    daa_results_df = daa_single_method,
    Group = td$Group,
    p_values_threshold = 1.0
  )

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)

  expected_cols <- c("feature", "group1", "group2",
                    "mean_rel_abundance_group1", "sd_rel_abundance_group1",
                    "mean_rel_abundance_group2", "sd_rel_abundance_group2",
                    "log2_fold_change", "p_adjust")
  expect_true(all(expected_cols %in% colnames(result)))
})

test_that("pathway_errorbar and pathway_errorbar_table share the same mean/sd source", {
  # Regression guard for the hand-rolled `pivot_longer %>% group_by %>%
  # summarise(mean(value), sd(value))` path that used to live inside
  # pathway_errorbar(). That path diverged from pathway_errorbar_table()
  # (which went through calculate_abundance_stats()) on NA handling and
  # was a latent bug: fix one path, miss the other. Both entry points now
  # route through summarize_abundance_by_group(), so the per-feature /
  # per-group mean and sd must be numerically identical.
  td <- create_errorbar_test_data(
    n_features = 5,
    p_adjust = c(0.01, 0.02, 0.03, 0.04, 0.05)
  )

  stats_tbl <- pathway_errorbar_table(
    abundance = td$abundance,
    daa_results_df = td$daa_results_df,
    Group = td$Group,
    p_values_threshold = 0.05
  )

  # Reach into the plot to grab the exact data frame ggplot consumed.
  plot <- pathway_errorbar(
    abundance = td$abundance,
    daa_results_df = td$daa_results_df,
    Group = td$Group,
    p_values_threshold = 0.05,
    x_lab = "pathway_name"
  )
  # The first patchwork panel is the errorbar plot itself; its $data holds
  # the long-format summary the new code feeds to ggplot.
  errorbar_data <- plot[[1]]$data

  # For each (feature, group) pair the plot used, the table must agree.
  group1_name <- unique(td$daa_results_df$group1)[1]
  group2_name <- unique(td$daa_results_df$group2)[1]
  for (feat in unique(as.character(errorbar_data$name))) {
    tbl_row <- stats_tbl[stats_tbl$feature == feat, , drop = FALSE]
    if (nrow(tbl_row) == 0) next

    plot_g1 <- errorbar_data[
      as.character(errorbar_data$name) == feat &
        as.character(errorbar_data$group) == group1_name,
    ]
    plot_g2 <- errorbar_data[
      as.character(errorbar_data$name) == feat &
        as.character(errorbar_data$group) == group2_name,
    ]

    expect_equal(plot_g1$mean, tbl_row$mean_rel_abundance_group1,
                 tolerance = 1e-12, info = feat)
    expect_equal(plot_g1$sd, tbl_row$sd_rel_abundance_group1,
                 tolerance = 1e-12, info = feat)
    expect_equal(plot_g2$mean, tbl_row$mean_rel_abundance_group2,
                 tolerance = 1e-12, info = feat)
    expect_equal(plot_g2$sd, tbl_row$sd_rel_abundance_group2,
                 tolerance = 1e-12, info = feat)
  }
})

test_that("pathway_errorbar rejects samples with zero total abundance", {
  # Regression: the relative-abundance conversion used to be
  # `apply(t(mat), 1, function(x) x / sum(x))`, which produced silent NaN
  # for any sample column whose total was 0. Those NaN propagated into
  # the plotted error bars without any warning. The zero-sum sample must
  # now surface as an actionable error that names the offending column.
  td <- create_errorbar_test_data(
    n_features = 6,
    p_adjust = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)
  )
  zero_sample <- colnames(td$abundance)[1]
  td$abundance[, zero_sample] <- 0

  expect_error(
    pathway_errorbar(
      abundance = td$abundance,
      daa_results_df = td$daa_results_df,
      Group = td$Group,
      p_values_threshold = 0.05,
      x_lab = "pathway_name"
    ),
    regexp = zero_sample
  )
})

test_that("pathway_errorbar preserves a method-native log2_fold_change instead of overwriting it with a mean-ratio", {
  # Regression: the function defensively added log2_fold_change as NA
  # only when missing, then unconditionally OVERWROTE every row with a
  # log2 of a relative-abundance mean ratio. For methods that ship
  # their own model-based effect size (DESeq2, edgeR, limma voom, LinDA,
  # Maaslin2, metagenomeSeq, and ALDEx2 with include_effect_size), the
  # bar in the side panel then disagreed with the p_adjust in the next
  # panel -- two outputs of the same fit telling different stories.
  # The fix only runs the mean-ratio fallback when no log2_fold_change
  # column is supplied.
  td <- create_errorbar_test_data(n_features = 5, p_adjust = rep(0.01, 5))

  # Distinctive sentinel values that the mean-ratio fallback cannot
  # coincidentally produce.
  sentinel <- c(99, -99, 77, -77, 42)
  td$daa_results_df$log2_fold_change <- sentinel

  p <- pathway_errorbar(
    abundance      = td$abundance,
    daa_results_df = td$daa_results_df,
    Group          = td$Group,
    p_values_threshold = 0.05,
    x_lab          = "pathway_name"
  )

  # The returned patchwork exposes the log2-fold-change panel data via
  # the main plot and one of the sub-patches; both must carry the user-
  # supplied values exactly. We compare sets because the plot reorders
  # features for display.
  plot_log2fc <- p$data$log2_fold_change
  expect_setequal(plot_log2fc, sentinel)
})

test_that("pathway_errorbar falls back to mean-ratio log2_fold_change when the column is absent", {
  # Guard: methods without a model-based effect size (ALDEx2 with
  # include_effect_size = FALSE, Lefser, user-supplied frames without
  # a log2_fold_change column) should still get a bar in the side
  # panel. The mean-ratio fallback must remain in place for that case.
  td <- create_errorbar_test_data(n_features = 5, p_adjust = rep(0.01, 5))
  # No log2_fold_change column on input.
  expect_false("log2_fold_change" %in% colnames(td$daa_results_df))

  p <- pathway_errorbar(
    abundance      = td$abundance,
    daa_results_df = td$daa_results_df,
    Group          = td$Group,
    p_values_threshold = 0.05,
    x_lab          = "pathway_name"
  )

  # Every feature should get a non-NA fallback log2_fold_change value.
  expect_true("log2_fold_change" %in% colnames(p$data))
  expect_false(any(is.na(p$data$log2_fold_change)))
})

# Regression: several theme() calls used `legend.position = "non"` (typo).
# In ggplot2 4.x any unrecognized string silently falls back to hiding the
# legend, so the plot looked correct -- but a future ggplot2 release that
# tightens this check would turn the typo into an error. Lock the spelling.
test_that("pathway_errorbar uses legend.position = 'none' (not 'non')", {
  body_src <- paste(deparse(body(ggpicrust2::pathway_errorbar)), collapse = "\n")
  expect_false(grepl("legend\\.position\\s*=\\s*\"non\"", body_src))
  expect_true(grepl("legend\\.position\\s*=\\s*\"none\"", body_src))
})
