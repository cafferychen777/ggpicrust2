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

  expect_message(
    pathway_errorbar(
      abundance = td$abundance,
      daa_results_df = td$daa_results_df,
      Group = td$Group,
      x_lab = "pathway_name"
    ),
    "pathways with missing annotations"
  )
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
