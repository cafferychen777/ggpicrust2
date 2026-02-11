# Tests for ggpicrust2() return structure

test_that("ggpicrust2 returns all expected fields", {
  skip_if(
    Sys.getenv("GGPICRUST2_RUN_E2E_TESTS", "false") != "true",
    "Set GGPICRUST2_RUN_E2E_TESTS=true to run full ggpicrust2 end-to-end tests."
  )
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

  # Backward compatibility: old API returned unnamed list elements accessible
  # by position. result[[1]] is the first DA method's plot+results.
  expect_true(!is.null(result[[1]]))
  expect_true(all(c("plot", "results") %in% names(result[[1]])))
})

test_that("ggpicrust2 aligns Group vector to abundance sample order before plotting", {
  captured_group <- NULL

  mock_abundance <- data.frame(
    S_B = c(10, 20),
    S_A = c(30, 40),
    row.names = c("K00001", "K00002"),
    check.names = FALSE
  )

  metadata <- data.frame(
    sample = c("S_A", "S_B"),
    Environment = c("A", "B"),
    stringsAsFactors = FALSE
  )

  mock_pathway_daa <- function(abundance, metadata, group, daa_method, select, p.adjust, reference) {
    data.frame(
      feature = rownames(abundance)[1],
      method = "mock_method",
      p_values = 0.01,
      adj_method = "BH",
      p_adjust = 0.01,
      group1 = "A",
      group2 = "B",
      stringsAsFactors = FALSE
    )
  }

  mock_pathway_annotation <- function(pathway, ko_to_kegg, daa_results_df, p_adjust_threshold) {
    daa_results_df$pathway_name <- daa_results_df$feature
    daa_results_df$pathway_class <- "Class1"
    daa_results_df
  }

  mock_pathway_errorbar <- function(abundance, daa_results_df, Group, ko_to_kegg, p_value_bar, order, colors, select, x_lab, p_values_threshold) {
    captured_group <<- Group
    NULL
  }

  result <- testthat::with_mocked_bindings(
    ggpicrust2(
      data = data.frame(function. = rownames(mock_abundance), mock_abundance, check.names = FALSE),
      metadata = metadata,
      group = "Environment",
      pathway = "KO",
      daa_method = "ALDEx2",
      ko_to_kegg = FALSE
    ),
    pathway_daa = mock_pathway_daa,
    pathway_annotation = mock_pathway_annotation,
    pathway_errorbar = mock_pathway_errorbar
  )

  expect_type(result, "list")
  expect_equal(names(captured_group), c("S_B", "S_A"))
  expect_equal(unname(captured_group), c("B", "A"))
})
