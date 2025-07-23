# Test for ggpicrust2_extended function
library(testthat)

# Helper function to create test data
create_test_data <- function() {
  # Create abundance data
  abundance <- matrix(rpois(300, lambda = 20), nrow = 30, ncol = 10)
  rownames(abundance) <- paste0("K", sprintf("%05d", 1:30))
  colnames(abundance) <- paste0("Sample", 1:10)
  abundance <- as.data.frame(abundance)
  abundance <- cbind(data.frame("#NAME" = rownames(abundance)), abundance)

  # Create metadata
  metadata <- data.frame(
    sample = paste0("Sample", 1:10),
    Environment = factor(rep(c("Forest", "Desert"), each = 5))
  )

  return(list(abundance = abundance, metadata = metadata))
}

test_that("ggpicrust2_extended runs without GSEA", {
  # Skip if required packages are not available
  skip_if_not_installed("dplyr")

  # Create test data
  test_data <- create_test_data()

  # Mock ggpicrust2 function
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    list(
      list(
        plot = ggplot2::ggplot(),
        results = data.frame(
          feature = paste0("path:ko", sprintf("%05d", 1:10)),
          p_adjust = runif(10, 0, 0.1),
          log_2_fold_change = rnorm(10)
        )
      )
    )
  })

  # Test without GSEA
  result <- ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Environment",
    pathway = "KO",
    daa_method = "ALDEx2",
    run_gsea = FALSE
  )

  # Check the result
  expect_type(result, "list")
  expect_true("daa_results" %in% names(result))
  expect_false("gsea_results" %in% names(result))
})

test_that("ggpicrust2_extended runs with GSEA", {
  # Skip if required packages are not available
  skip_if_not_installed("dplyr")
  skip_if_not_installed("fgsea")

  # Create test data
  test_data <- create_test_data()

  # Mock ggpicrust2 function
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    list(
      list(
        plot = ggplot2::ggplot(),
        results = data.frame(
          feature = paste0("path:ko", sprintf("%05d", 1:10)),
          p_adjust = runif(10, 0, 0.1),
          log_2_fold_change = rnorm(10)
        )
      )
    )
  })

  # Mock pathway_gsea function
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    data.frame(
      pathway_id = paste0("path:ko", sprintf("%05d", 1:10)),
      pathway_name = paste0("path:ko", sprintf("%05d", 1:10)),
      size = sample(10:100, 10, replace = TRUE),
      ES = runif(10, -0.8, 0.8),
      NES = runif(10, -2, 2),
      pvalue = runif(10, 0, 0.1),
      p.adjust = runif(10, 0, 0.2),
      leading_edge = replicate(10, paste(paste0("K", sprintf("%05d", sample(1:1000, 5))), collapse = ";")),
      method = rep("fgsea", 10),
      stringsAsFactors = FALSE
    )
  })

  # Mock gsea_pathway_annotation function
  mockery::stub(ggpicrust2_extended, "gsea_pathway_annotation", function(...) {
    data.frame(
      pathway_id = paste0("path:ko", sprintf("%05d", 1:10)),
      pathway_name = paste("KEGG Pathway", 1:10),
      pathway_class = rep(c("Metabolism", "Genetic Information Processing"), each = 5),
      size = sample(10:100, 10, replace = TRUE),
      ES = runif(10, -0.8, 0.8),
      NES = runif(10, -2, 2),
      pvalue = runif(10, 0, 0.1),
      p.adjust = runif(10, 0, 0.2),
      leading_edge = replicate(10, paste(paste0("K", sprintf("%05d", sample(1:1000, 5))), collapse = ";")),
      method = rep("fgsea", 10),
      stringsAsFactors = FALSE
    )
  })

  # Mock visualize_gsea function
  mockery::stub(ggpicrust2_extended, "visualize_gsea", function(...) {
    ggplot2::ggplot()
  })

  # Mock compare_gsea_daa function
  mockery::stub(ggpicrust2_extended, "compare_gsea_daa", function(...) {
    list(
      plot = ggplot2::ggplot(),
      results = list(
        overlap = c("path:ko00001", "path:ko00002"),
        gsea_only = c("path:ko00003", "path:ko00004"),
        daa_only = c("path:ko00005", "path:ko00006"),
        n_overlap = 2,
        n_gsea_only = 2,
        n_daa_only = 2
      )
    )
  })

  # Test with GSEA
  result <- ggpicrust2_extended(
    data = test_data$abundance,
    metadata = test_data$metadata,
    group = "Environment",
    pathway = "KO",
    daa_method = "ALDEx2",
    run_gsea = TRUE,
    gsea_params = list(
      method = "fgsea",
      rank_method = "signal2noise"
    )
  )

  # Check the result
  expect_type(result, "list")
  expect_true("daa_results" %in% names(result))
  expect_true("gsea_results" %in% names(result))
  expect_true("gsea_plot" %in% names(result))
  expect_true("comparison" %in% names(result))
})

test_that("ggpicrust2_extended handles missing packages", {
  # Create test data
  test_data <- create_test_data()

  # Mock ggpicrust2 function
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    list(
      list(
        plot = ggplot2::ggplot(),
        results = data.frame(
          feature = paste0("path:ko", sprintf("%05d", 1:10)),
          p_adjust = runif(10, 0, 0.1),
          log_2_fold_change = rnorm(10)
        )
      )
    )
  })

  # Mock requireNamespace to simulate missing package
  mockery::stub(ggpicrust2_extended, "requireNamespace", function(pkg, ...) {
    return(FALSE)
  })

  # Test with GSEA but missing package
  expect_warning(
    result <- ggpicrust2_extended(
      data = test_data$abundance,
      metadata = test_data$metadata,
      group = "Environment",
      pathway = "KO",
      daa_method = "ALDEx2",
      run_gsea = TRUE
    ),
    "Package 'fgsea' is required for GSEA analysis. Skipping GSEA."
  )

  # Check the result
  expect_type(result, "list")
  expect_true("daa_results" %in% names(result))
  expect_false("gsea_results" %in% names(result))
})

test_that("ggpicrust2_extended validates inputs correctly", {
  # Create test data
  test_data <- create_test_data()

  # Mock ggpicrust2 function
  mockery::stub(ggpicrust2_extended, "ggpicrust2", function(...) {
    list(
      list(
        plot = ggplot2::ggplot(),
        results = data.frame(
          feature = paste0("path:ko", sprintf("%05d", 1:10)),
          p_adjust = runif(10, 0, 0.1),
          log_2_fold_change = rnorm(10)
        )
      )
    )
  })

  # Test with missing abundance data
  # We need to mock the function to avoid the warning about missing fgsea package
  mockery::stub(ggpicrust2_extended, "requireNamespace", function(pkg, ...) {
    return(TRUE)  # Pretend all packages are available
  })

  # Mock pathway_gsea to avoid actual computation
  mockery::stub(ggpicrust2_extended, "pathway_gsea", function(...) {
    data.frame()
  })

  # Now test with missing abundance data
  expect_error(
    ggpicrust2_extended(
      metadata = test_data$metadata,
      group = "Environment",
      pathway = "KO",
      daa_method = "ALDEx2",
      run_gsea = TRUE
    ),
    "No abundance data provided for GSEA analysis"
  )

  # Test with invalid gsea_params
  # We need to check if modifyList is called with invalid parameters
  # Instead of testing the error directly, we'll check that the function validates gsea_params
  expect_error(
    ggpicrust2_extended(
      data = test_data$abundance,
      metadata = test_data$metadata,
      group = "Environment",
      pathway = "KO",
      daa_method = "ALDEx2",
      run_gsea = TRUE,
      gsea_params = "not a list"
    )
  )
})
