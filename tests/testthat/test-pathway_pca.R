# Tests for pathway_pca function

# Helper to create test data
create_pca_test_data <- function(n_pathways = 3, n_samples = 10, n_groups = 2) {
  set.seed(123)
  test_abundance <- matrix(rnorm(n_pathways * n_samples), nrow = n_pathways, ncol = n_samples)
  colnames(test_abundance) <- paste0("Sample", 1:n_samples)
  rownames(test_abundance) <- paste0("Pathway", LETTERS[1:n_pathways])

  groups <- rep(paste0("Group", 1:n_groups), length.out = n_samples)
  test_metadata <- data.frame(
    sample_name = colnames(test_abundance),
    group = factor(groups)
  )

  list(abundance = test_abundance, metadata = test_metadata)
}

test_that("pathway_pca works with basic inputs", {
  data <- create_pca_test_data()
  result <- pathway_pca(data$abundance, data$metadata, "group")
  expect_s3_class(result, "ggplot")
})

test_that("pathway_pca works with custom colors", {
  data <- create_pca_test_data()
  result <- pathway_pca(data$abundance, data$metadata, "group", colors = c("red", "blue"))
  expect_s3_class(result, "ggplot")
})

test_that("pathway_pca works with multiple groups", {
  data <- create_pca_test_data(n_samples = 12, n_groups = 3)
  result <- pathway_pca(data$abundance, data$metadata, "group")
  expect_s3_class(result, "ggplot")
})

test_that("pathway_pca validates inputs", {
  data <- create_pca_test_data()


  # Invalid input types
  expect_error(pathway_pca(list(1,2,3), data$metadata, "group"), "must be a data frame or matrix")
  expect_error(pathway_pca(data$abundance, list(a=1), "group"), "must be a data frame")

  # Missing group column
  wrong_metadata <- data.frame(sample_name = data$metadata$sample_name, other = 1:10)
  expect_error(pathway_pca(data$abundance, wrong_metadata, "group"), "Group column.*not found")

  # NA values
  data$abundance[1,1] <- NA
  expect_error(pathway_pca(data$abundance, data$metadata, "group"), "NA|missing")
})

test_that("pathway_pca accepts finite zero-sum sample columns", {
  abundance <- matrix(
    c(
      -1, 1, 2, 4,
       0, 2, 3, 5,
       1, 3, 4, 6
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Pathway", 1:3), paste0("Sample", 1:4))
  )
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    p <- pathway_pca(abundance, metadata, "group", show_marginal = FALSE),
    "Skipping PCA confidence ellipse"
  )
  expect_s3_class(p, "ggplot")
})

test_that("pathway_pca keeps zero-variance sample profiles as observations", {
  abundance <- matrix(
    c(
      5, 5, 1, 2, 3,
      5, 5, 2, 3, 4,
      5, 5, 3, 4, 5
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Pathway", 1:3), paste0("Sample", 1:5))
  )
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = c("A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    p <- pathway_pca(abundance, metadata, "group", show_marginal = FALSE),
    "fewer than 4 samples: A=2, B=3"
  )
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), ncol(abundance))
  expect_setequal(as.character(p$data$Group), c("A", "B"))
})

test_that("pathway_pca skips confidence ellipses for groups with fewer than four samples", {
  abundance <- matrix(
    c(
      1, 2, 3, 4, 5, 6,
      2, 4, 6, 8, 10, 13,
      1, 3, 5, 7, 9, 8
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Pathway", 1:3), paste0("Sample", 1:6))
  )
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = c("A", "A", "B", "B", "B", "B"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    p <- pathway_pca(abundance, metadata, "group", show_marginal = FALSE),
    "fewer than 4 samples: A=2"
  )
  expect_s3_class(p, "ggplot")
  ellipse_layers <- vapply(
    p$layers,
    function(layer) inherits(layer$stat, "StatEllipse"),
    logical(1)
  )
  expect_equal(sum(ellipse_layers), 1)
  expect_true(all(as.character(p$layers[[which(ellipse_layers)]]$data$Group) == "B"))
})

test_that("pathway_pca rejects missing group labels after sample alignment", {
  abundance <- matrix(
    c(
      1, 2, 3, 4,
      2, 3, 4, 5,
      3, 4, 5, 6
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(paste0("Pathway", 1:3), paste0("Sample", 1:4))
  )
  metadata <- data.frame(
    sample_name = colnames(abundance),
    group = c("A", "A", "B", NA),
    stringsAsFactors = FALSE
  )

  expect_error(
    pathway_pca(abundance, metadata, "group", show_marginal = FALSE),
    "non-missing, non-empty group labels.*Sample4"
  )
})

test_that("pathway_pca throws error with wrong color count", {
  data <- create_pca_test_data()
  expect_error(
    pathway_pca(data$abundance, data$metadata, "group", colors = c("red", "blue", "green")),
    "Number of colors"
  )
})

test_that("pathway_pca show_marginal parameter works", {
  data <- create_pca_test_data()

  result_with <- pathway_pca(data$abundance, data$metadata, "group", show_marginal = TRUE)
  result_without <- pathway_pca(data$abundance, data$metadata, "group", show_marginal = FALSE)

  expect_s3_class(result_with, "ggplot")
  expect_s3_class(result_without, "ggplot")
})

# Regression: the marginal density panels used to attach
# scale_y_discrete() to geom_density(), which produces a continuous y.
# ggplot2 silently tolerated that mismatch but it was a latent type bug
# and mis-interpreted `expand = c(0, 0.001)` in discrete-category units.
# Assert the source no longer attaches a discrete scale to the density
# y aesthetic.
test_that("pathway_pca marginal density uses a continuous y scale", {
  body_src <- paste(deparse(body(ggpicrust2::pathway_pca)), collapse = "\n")
  # The continuous scale must be present on the density panels.
  expect_true(grepl("scale_y_continuous\\(", body_src))
  # And the discrete scale must be gone from the density construction.
  expect_false(grepl("scale_y_discrete\\(", body_src))
})
