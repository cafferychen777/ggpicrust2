# Tests for ggpicrust2() return structure
# Tests the new fields: abundance, metadata, group, daa_results_df, ko_to_kegg

test_that("ggpicrust2 returns all expected fields with ko_to_kegg=TRUE", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")

  # Load example data
  data(ko_abundance, package = "ggpicrust2")
  data(metadata, package = "ggpicrust2")

  # Run ggpicrust2 with ko_to_kegg=TRUE
  result <- suppressMessages(ggpicrust2(
    data = ko_abundance,
    metadata = metadata,
    group = "Environment",
    pathway = "KO",
    daa_method = "ALDEx2",
    ko_to_kegg = TRUE,
    order = "pathway_class",
    p_values_bar = TRUE,
    x_lab = "pathway_name"
  ))

  # Check that result is a list
  expect_type(result, "list")

  # Check for new fields
  expect_true("abundance" %in% names(result),
              info = "Result should contain 'abundance' field")
  expect_true("metadata" %in% names(result),
              info = "Result should contain 'metadata' field")
  expect_true("group" %in% names(result),
              info = "Result should contain 'group' field")
  expect_true("daa_results_df" %in% names(result),
              info = "Result should contain 'daa_results_df' field")
  expect_true("ko_to_kegg" %in% names(result),
              info = "Result should contain 'ko_to_kegg' field")

  # Check field types
  expect_true(is.data.frame(result$abundance) || is.matrix(result$abundance),
              info = "abundance should be a data.frame or matrix")
  expect_true(is.data.frame(result$metadata),
              info = "metadata should be a data.frame")
  expect_equal(result$group, "Environment",
               info = "group should match the input")
  expect_true(is.data.frame(result$daa_results_df),
              info = "daa_results_df should be a data.frame")
  expect_true(result$ko_to_kegg,
              info = "ko_to_kegg should be TRUE")

  # Check backward compatibility - numbered elements should still exist
  expect_true(!is.null(result[[1]]),
              info = "Result should have numbered element [[1]]")
  expect_true("plot" %in% names(result[[1]]),
              info = "Result[[1]] should contain 'plot'")
  expect_true("results" %in% names(result[[1]]),
              info = "Result[[1]] should contain 'results'")
})

test_that("ggpicrust2 returns all expected fields with ko_to_kegg=FALSE", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")

  # Load example data
  data(metacyc_abundance, package = "ggpicrust2")
  data(metadata, package = "ggpicrust2")

  # Run ggpicrust2 with ko_to_kegg=FALSE
  result <- suppressMessages(ggpicrust2(
    data = metacyc_abundance,
    metadata = metadata,
    group = "Environment",
    pathway = "MetaCyc",
    daa_method = "ALDEx2",
    ko_to_kegg = FALSE,
    order = "group",
    p_values_bar = TRUE,
    x_lab = "description"
  ))

  # Check that result is a list
  expect_type(result, "list")

  # Check for new fields
  expect_true("abundance" %in% names(result))
  expect_true("metadata" %in% names(result))
  expect_true("group" %in% names(result))
  expect_true("daa_results_df" %in% names(result))
  expect_true("ko_to_kegg" %in% names(result))

  # Check ko_to_kegg value
  expect_false(result$ko_to_kegg,
               info = "ko_to_kegg should be FALSE")

  # Check that abundance has correct structure
  expect_true(nrow(result$abundance) > 0,
              info = "abundance should have rows")
  expect_true(ncol(result$abundance) > 0,
              info = "abundance should have columns")
})

test_that("ggpicrust2 returned data works with pathway_pca", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")
  skip_if_not_installed("FactoMineR")
  skip_if_not_installed("factoextra")

  # Load example data
  data(ko_abundance, package = "ggpicrust2")
  data(metadata, package = "ggpicrust2")

  # Run ggpicrust2
  result <- suppressMessages(ggpicrust2(
    data = ko_abundance,
    metadata = metadata,
    group = "Environment",
    pathway = "KO",
    daa_method = "ALDEx2",
    ko_to_kegg = TRUE
  ))

  # Use returned data with pathway_pca
  pca_plot <- suppressMessages(pathway_pca(
    abundance = result$abundance,
    metadata = result$metadata,
    group = result$group
  ))

  # Check that PCA plot was created successfully
  expect_s3_class(pca_plot, "ggplot")
})

test_that("ggpicrust2 returned data works with pathway_heatmap for significant features", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")
  skip_if_not_installed("pheatmap")

  # Load example data
  data(ko_abundance, package = "ggpicrust2")
  data(metadata, package = "ggpicrust2")

  # Run ggpicrust2
  result <- suppressMessages(ggpicrust2(
    data = ko_abundance,
    metadata = metadata,
    group = "Environment",
    pathway = "KO",
    daa_method = "ALDEx2",
    ko_to_kegg = TRUE
  ))

  # Filter significant features
  sig_features <- result$daa_results_df %>%
    dplyr::filter(p_adjust < 0.05) %>%
    dplyr::pull(feature)

  # Only test if there are significant features
  if (length(sig_features) > 0) {
    # Use returned data with pathway_heatmap
    heatmap_result <- suppressMessages(pathway_heatmap(
      abundance = result$abundance[sig_features, , drop = FALSE],
      metadata = result$metadata,
      group = result$group
    ))

    # Check that heatmap was created (returns list with plot element)
    expect_true(!is.null(heatmap_result))
  } else {
    # If no significant features, skip heatmap test
    skip("No significant features found for heatmap test")
  }
})

test_that("ggpicrust2 daa_results_df contains required columns", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")

  # Load example data
  data(ko_abundance, package = "ggpicrust2")
  data(metadata, package = "ggpicrust2")

  # Run ggpicrust2
  result <- suppressMessages(ggpicrust2(
    data = ko_abundance,
    metadata = metadata,
    group = "Environment",
    pathway = "KO",
    daa_method = "ALDEx2",
    ko_to_kegg = TRUE
  ))

  # Check that daa_results_df has expected columns
  expect_true("feature" %in% names(result$daa_results_df))
  expect_true("p_adjust" %in% names(result$daa_results_df))
  expect_true("method" %in% names(result$daa_results_df))
})

test_that("ggpicrust2 abundance matches sample columns to metadata rows", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")

  # Load example data
  data(ko_abundance, package = "ggpicrust2")
  data(metadata, package = "ggpicrust2")

  # Run ggpicrust2
  result <- suppressMessages(ggpicrust2(
    data = ko_abundance,
    metadata = metadata,
    group = "Environment",
    pathway = "KO",
    daa_method = "ALDEx2",
    ko_to_kegg = TRUE
  ))

  # Check that abundance columns match metadata rows
  abundance_samples <- colnames(result$abundance)

  # Metadata should have sample identifiers
  expect_true(nrow(result$metadata) > 0)

  # Number of samples in abundance should match
  expect_equal(length(abundance_samples), nrow(result$metadata))
})

test_that("ggpicrust2 metadata includes group column", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")

  # Load example data
  data(ko_abundance, package = "ggpicrust2")
  data(metadata, package = "ggpicrust2")

  # Run ggpicrust2
  result <- suppressMessages(ggpicrust2(
    data = ko_abundance,
    metadata = metadata,
    group = "Environment",
    pathway = "KO",
    daa_method = "ALDEx2",
    ko_to_kegg = TRUE
  ))

  # Check that metadata contains the group column
  expect_true(result$group %in% names(result$metadata),
              info = "metadata should contain the group column")
})
