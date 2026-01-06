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
  data <- create_pca_test_data(n_groups = 3)
  result <- pathway_pca(data$abundance, data$metadata, "group")
  expect_s3_class(result, "ggplot")
})

test_that("pathway_pca validates inputs", {
  data <- create_pca_test_data()

  # Missing arguments
  expect_error(pathway_pca(), "Abundance matrix is required")
  expect_error(pathway_pca(data$abundance), "Metadata is required")
  expect_error(pathway_pca(data$abundance, data$metadata), "Group variable name is required")

  # Invalid input types
  expect_error(pathway_pca(list(1,2,3), data$metadata, "group"), "Abundance must be a matrix")
  expect_error(pathway_pca(data$abundance, list(a=1), "group"), "Metadata must be a data frame")

  # Missing group column
  wrong_metadata <- data.frame(sample_name = data$metadata$sample_name, other = 1:10)
  expect_error(pathway_pca(data$abundance, wrong_metadata, "group"), "Group column.*not found")

  # NA values
  data$abundance[1,1] <- NA
  expect_error(pathway_pca(data$abundance, data$metadata, "group"))
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
