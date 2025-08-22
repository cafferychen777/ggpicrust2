# Comprehensive tests for pathway_errorbar()

library(testthat)
library(ggplot2)
library(patchwork)
library(dplyr)

skip_if_not_installed("ggprism")

make_dummy_data <- function(n_features = 8, n_samples = 10,
                            groups = c("GroupA", "GroupB"),
                            add_pathway_class = TRUE,
                            all_significant = TRUE) {
  set.seed(123)
  abundance <- matrix(runif(n_features * n_samples, min = 0, max = 1),
                      nrow = n_features, ncol = n_samples)
  rownames(abundance) <- paste0("pathway", seq_len(n_features))
  colnames(abundance) <- paste0("sample", seq_len(n_samples))

  group_df <- data.frame(
    sample = colnames(abundance),
    group  = rep(groups, length.out = n_samples),
    stringsAsFactors = FALSE
  )
  Group <- factor(group_df$group, levels = groups)
  names(Group) <- group_df$sample

  p_vals <- if (all_significant) rep(0.01, n_features) else seq(0.01, 0.2, length.out = n_features)

  daa <- data.frame(
    feature = rownames(abundance),
    pathway_name = paste0("Pathway ", seq_len(n_features)),
    description = paste0("Description ", seq_len(n_features)),
    p_adjust = p_vals,
    method = rep("ALDEx2_Welch's t test", n_features),
    group1 = rep(groups[1], n_features),
    group2 = rep(groups[2], n_features),
    stringsAsFactors = FALSE
  )
  if (add_pathway_class) {
    daa$pathway_class <- rep(c("Class1", "Class2"), length.out = n_features)
  }

  list(abundance = abundance, daa = daa, Group = Group)
}


test_that("valid basic usage returns patchwork object (ko_to_kegg = FALSE)", {
  dd <- make_dummy_data(add_pathway_class = FALSE)
  p <- pathway_errorbar(
    abundance = dd$abundance,
    daa_results_df = dd$daa,
    Group = dd$Group,
    ko_to_kegg = FALSE,
    select = rownames(dd$abundance)[1:5],
    x_lab = "pathway_name"
  )
  expect_s3_class(p, "patchwork")
})


test_that("input validations: type and columns", {
  dd <- make_dummy_data()
  # abundance type check
  expect_error(pathway_errorbar(
    abundance = as.list(dd$abundance),
    daa_results_df = dd$daa,
    Group = dd$Group
  ), "must be a matrix or data frame")

  # daa_results_df type check
  expect_error(pathway_errorbar(
    abundance = dd$abundance,
    daa_results_df = as.list(dd$daa),
    Group = dd$Group
  ), "must be a data frame")

  # missing required columns
  bad_daa <- dd$daa %>% select(-method)
  expect_error(pathway_errorbar(
    abundance = dd$abundance,
    daa_results_df = bad_daa,
    Group = dd$Group
  ), "Missing required columns")
})


test_that("input validations: Group length and multiple methods/groups", {
  dd <- make_dummy_data()
  # Group length mismatch
  bad_Group <- dd$Group[-1]
  expect_error(pathway_errorbar(
    abundance = dd$abundance,
    daa_results_df = dd$daa,
    Group = bad_Group
  ), "Length of Group must match")

  # multiple methods
  mm <- dd$daa
  mm$method[1] <- "Other"
  expect_error(pathway_errorbar(
    abundance = dd$abundance,
    daa_results_df = mm,
    Group = dd$Group
  ), "contains more than one method")

  # multiple group1 / group2 values
  mg <- dd$daa
  mg$group1[1] <- "AnotherGroup"
  expect_error(pathway_errorbar(
    abundance = dd$abundance,
    daa_results_df = mg,
    Group = dd$Group
  ), "contains more than one group")
})


test_that("x_lab handling and messages", {
  dd <- make_dummy_data(add_pathway_class = FALSE)
  # invalid x_lab should error during column selection
  expect_error(pathway_errorbar(
    abundance = dd$abundance,
    daa_results_df = dd$daa,
    Group = dd$Group,
    x_lab = "not_a_col"
  ), regexp = "undefined columns selected")

  # exclude NA annotations and message when missing
  daa_missing <- dd$daa
  daa_missing$pathway_name[1] <- NA
  expect_message(pathway_errorbar(
    abundance = dd$abundance,
    daa_results_df = daa_missing,
    Group = dd$Group,
    x_lab = "pathway_name"
  ), "missing annotations")
})


test_that("threshold, select, and empty after filter", {
  dd <- make_dummy_data(all_significant = FALSE)
  # select subset
  sel <- rownames(dd$abundance)[1:3]
  p <- pathway_errorbar(
    abundance = dd$abundance,
    daa_results_df = dd$daa,
    Group = dd$Group,
    p_values_threshold = 0.2,
    select = sel,
    x_lab = "pathway_name"
  )
  expect_s3_class(p, "patchwork")

  # too strict threshold -> no features
  expect_error(pathway_errorbar(
    abundance = dd$abundance,
    daa_results_df = dd$daa,
    Group = dd$Group,
    p_values_threshold = 0.0001,
    x_lab = "pathway_name"
  ), "cannot be performed because there are no features")
})


test_that("max_features emits warning and can proceed", {
  dd <- make_dummy_data(n_features = 6)
  expect_warning({
    p <- pathway_errorbar(
      abundance = dd$abundance,
      daa_results_df = dd$daa,
      Group = dd$Group,
      x_lab = "pathway_name",
      max_features = 3
    )
    expect_s3_class(p, "patchwork")
  }, "exceeds 3")
})


test_that("p_value_bar toggles when more than two group columns exist", {
  dd <- make_dummy_data()
  daa3 <- dd$daa
  # add a third group column to trigger message and auto-disable
  daa3$group3 <- dd$daa$group1
  expect_message(
    pathway_errorbar(
      abundance = dd$abundance,
      daa_results_df = daa3,
      Group = dd$Group,
      x_lab = "pathway_name",
      p_value_bar = TRUE
    ), "automatically set to FALSE"
  )
})


test_that("legend parameters, p-value formatting and colors work", {
  dd <- make_dummy_data()
  formats <- c("numeric", "scientific", "smart", "stars_only", "combined")
  for (fmt in formats) {
    p <- pathway_errorbar(
      abundance = dd$abundance,
      daa_results_df = dd$daa,
      Group = dd$Group,
      x_lab = "pathway_name",
      legend_position = "bottom",
      legend_direction = "horizontal",
      legend_title = "Groups",
      legend_title_size = 10,
      legend_text_size = 8,
      legend_key_size = 0.6,
      pvalue_format = fmt,
      pvalue_stars = TRUE,
      pvalue_colors = TRUE
    )
    expect_s3_class(p, "patchwork")
  }
})


test_that("smart_colors and accessibility_mode work", {
  dd <- make_dummy_data()
  p <- pathway_errorbar(
    abundance = dd$abundance,
    daa_results_df = dd$daa,
    Group = dd$Group,
    x_lab = "pathway_name",
    smart_colors = TRUE,
    accessibility_mode = TRUE
  )
  expect_s3_class(p, "patchwork")
})


test_that("pathway class annotation and positioning works when ko_to_kegg = TRUE", {
  dd <- make_dummy_data(add_pathway_class = TRUE)
  pos_opts <- c("left", "right", "none")
  for (pos in pos_opts) {
    p <- pathway_errorbar(
      abundance = dd$abundance,
      daa_results_df = dd$daa,
      Group = dd$Group,
      ko_to_kegg = TRUE,
      x_lab = "pathway_name",
      pathway_class_position = pos,
      pathway_class_text_size = 5
    )
    expect_s3_class(p, "patchwork")
  }
})


test_that("pathway_names_text_size works (auto and custom)", {
  dd <- make_dummy_data()
  p1 <- pathway_errorbar(
    abundance = dd$abundance,
    daa_results_df = dd$daa,
    Group = dd$Group,
    x_lab = "pathway_name",
    pathway_names_text_size = "auto"
  )
  p2 <- pathway_errorbar(
    abundance = dd$abundance,
    daa_results_df = dd$daa,
    Group = dd$Group,
    x_lab = "pathway_name",
    pathway_names_text_size = 14
  )
  expect_s3_class(p1, "patchwork")
  expect_s3_class(p2, "patchwork")
})

