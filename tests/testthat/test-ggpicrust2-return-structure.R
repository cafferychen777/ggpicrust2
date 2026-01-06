# Tests for ggpicrust2() return structure

test_that("ggpicrust2 returns all expected fields", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")

  data(ko_abundance, package = "ggpicrust2")
  data(metadata, package = "ggpicrust2")

  result <- suppressMessages(ggpicrust2(
    data = ko_abundance,
    metadata = metadata,
    group = "Environment",
    pathway = "KO",
    daa_method = "ALDEx2",
    ko_to_kegg = TRUE
  ))

  # Check structure
  expect_type(result, "list")
  expect_true(all(c("abundance", "metadata", "group", "daa_results_df", "ko_to_kegg") %in% names(result)))

  # Check field types
  expect_true(is.data.frame(result$abundance) || is.matrix(result$abundance))
  expect_true(is.data.frame(result$metadata))
  expect_true(is.data.frame(result$daa_results_df))
  expect_equal(result$group, "Environment")
  expect_true(result$ko_to_kegg)

  # Check backward compatibility
  expect_true(!is.null(result[[1]]))
  expect_true(all(c("plot", "results") %in% names(result[[1]])))
})

test_that("ggpicrust2 returned data works with downstream functions", {
  skip_on_cran()
  skip_if_not_installed("ALDEx2")
  skip_if_not_installed("FactoMineR")
  skip_if_not_installed("factoextra")

  data(ko_abundance, package = "ggpicrust2")
  data(metadata, package = "ggpicrust2")

  result <- suppressMessages(ggpicrust2(
    data = ko_abundance,
    metadata = metadata,
    group = "Environment",
    pathway = "KO",
    daa_method = "ALDEx2",
    ko_to_kegg = TRUE
  ))

  # Test with pathway_pca
  pca_plot <- suppressMessages(pathway_pca(
    abundance = result$abundance,
    metadata = result$metadata,
    group = result$group
  ))
  expect_s3_class(pca_plot, "ggplot")
})
