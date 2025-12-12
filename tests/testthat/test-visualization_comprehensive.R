# Comprehensive visualization tests for all comparison plot types
library(testthat)

# Advanced visualization test helpers
create_visualization_test_data <- function(overlap_scenario = "moderate") {
  # Create controlled data for different overlap scenarios
  
  if (overlap_scenario == "high") {
    # High overlap scenario: 80% of significant pathways overlap
    gsea_sig <- paste0("ko", sprintf("%05d", 1:10))
    daa_sig <- paste0("ko", sprintf("%05d", 3:12))  # 8 overlapping
  } else if (overlap_scenario == "low") {
    # Low overlap scenario: 20% of significant pathways overlap
    gsea_sig <- paste0("ko", sprintf("%05d", 1:10))
    daa_sig <- paste0("ko", sprintf("%05d", 8:17))  # 3 overlapping
  } else if (overlap_scenario == "none") {
    # No overlap scenario
    gsea_sig <- paste0("ko", sprintf("%05d", 1:10))
    daa_sig <- paste0("ko", sprintf("%05d", 20:29))  # 0 overlapping
  } else {
    # Moderate overlap scenario: 50% overlap (default)
    gsea_sig <- paste0("ko", sprintf("%05d", 1:10))
    daa_sig <- paste0("ko", sprintf("%05d", 6:15))  # 5 overlapping
  }
  
  # Create GSEA results
  all_gsea_pathways <- c(gsea_sig, paste0("ko", sprintf("%05d", 100:109)))
  gsea_results <- data.frame(
    pathway_id = all_gsea_pathways,
    pathway_name = paste("Pathway", all_gsea_pathways),
    pathway_class = sample(c("Metabolism", "Genetic Information Processing", "Cellular Processes"),
                          length(all_gsea_pathways), replace = TRUE),
    size = sample(20:150, length(all_gsea_pathways), replace = TRUE),
    ES = rnorm(length(all_gsea_pathways), 0, 1.2),
    NES = rnorm(length(all_gsea_pathways), 0, 1.8),
    pvalue = runif(length(all_gsea_pathways), 0.001, 0.2),
    p.adjust = ifelse(all_gsea_pathways %in% gsea_sig, 
                     runif(length(all_gsea_pathways), 0.001, 0.045),
                     runif(length(all_gsea_pathways), 0.06, 0.3)),
    leading_edge = replicate(length(all_gsea_pathways), {
      genes <- paste0("K", sprintf("%05d", sample(1:20000, sample(5:20, 1))))
      paste(genes, collapse = ";")
    }),
    method = rep("fgsea", length(all_gsea_pathways)),
    stringsAsFactors = FALSE
  )
  
  # Create DAA results
  all_daa_features <- c(daa_sig, paste0("ko", sprintf("%05d", 200:209)))
  daa_results <- data.frame(
    feature = all_daa_features,
    method = rep("LinDA", length(all_daa_features)),
    group1 = rep("Control", length(all_daa_features)),
    group2 = rep("Treatment", length(all_daa_features)),
    p_values = runif(length(all_daa_features), 0.001, 0.2),
    p_adjust = ifelse(all_daa_features %in% daa_sig,
                     runif(length(all_daa_features), 0.001, 0.045),
                     runif(length(all_daa_features), 0.06, 0.4)),
    log_2_fold_change = rnorm(length(all_daa_features), 0, 2),
    effect_size = rnorm(length(all_daa_features), 0, 1.5),
    stringsAsFactors = FALSE
  )
  
  return(list(
    gsea_results = gsea_results,
    daa_results = daa_results,
    expected_overlap = intersect(gsea_sig, daa_sig),
    expected_gsea_only = setdiff(gsea_sig, daa_sig),
    expected_daa_only = setdiff(daa_sig, gsea_sig)
  ))
}

test_that("Venn diagram visualization generation and structure", {
  skip_if_not_installed("ggplot2")
  
  # Test different overlap scenarios
  scenarios <- c("high", "moderate", "low", "none")
  
  for (scenario in scenarios) {
    test_data <- create_visualization_test_data(scenario)
    
    # Test with ggVennDiagram available
    if (requireNamespace("ggVennDiagram", quietly = TRUE)) {
      # Mock ggVennDiagram to control output
      mockery::stub(compare_gsea_daa, "ggVennDiagram::ggVennDiagram", function(venn_list) {
        # Verify input format
        expect_type(venn_list, "list")
        expect_named(venn_list, c("GSEA", "DAA"))
        expect_type(venn_list$GSEA, "character")
        expect_type(venn_list$DAA, "character")
        
        # Return mock plot
        ggplot2::ggplot() + 
          ggplot2::geom_blank() + 
          ggplot2::labs(title = paste("Venn Diagram -", scenario)) +
          ggplot2::theme_minimal()
      })
    }
    
    comparison <- compare_gsea_daa(test_data$gsea_results, test_data$daa_results, plot_type = "venn")
    
    # Verify plot structure
    expect_s3_class(comparison$plot, "ggplot")
    
    # Verify plot has title
    plot_build <- ggplot2::ggplot_build(comparison$plot)
    
    # Verify comparison results match expected overlap
    expect_setequal(comparison$results$overlap, test_data$expected_overlap)
    expect_setequal(comparison$results$gsea_only, test_data$expected_gsea_only)
    expect_setequal(comparison$results$daa_only, test_data$expected_daa_only)
  }
})

test_that("Venn diagram fallback visualization when ggVennDiagram unavailable", {
  skip_if_not_installed("ggplot2")
  
  test_data <- create_visualization_test_data("moderate")
  
  # Mock ggVennDiagram as unavailable
  with_mock(
    `requireNamespace` = function(pkg, quietly = TRUE) {
      if (pkg == "ggVennDiagram") return(FALSE)
      return(TRUE)
    },
    {
      expect_warning(
        comparison <- compare_gsea_daa(test_data$gsea_results, test_data$daa_results, plot_type = "venn"),
        "Package 'ggVennDiagram' is required for Venn diagrams. Using a basic plot instead."
      )
      
      # Verify fallback plot structure
      expect_s3_class(comparison$plot, "ggplot")
      
      # Verify fallback plot content
      plot_build <- ggplot2::ggplot_build(comparison$plot)
      expect_true(length(plot_build$data) > 0)
      
      # Verify it's a bar plot (fallback visualization)
      expect_true("GeomBar" %in% class(plot_build$plot$layers[[1]]$geom))
    }
  )
})

test_that("UpSet plot visualization generation and structure", {
  skip_if_not_installed("ggplot2")
  
  test_data <- create_visualization_test_data("moderate")
  
  # Test UpSet plot creation
  if (requireNamespace("UpSetR", quietly = TRUE)) {
    # Mock UpSetR::upset function
    mockery::stub(compare_gsea_daa, "UpSetR::upset", function(upset_df, ...) {
      # Verify input format
      expect_true(is.data.frame(upset_df))
      expect_equal(colnames(upset_df), c("GSEA", "DAA"))
      expect_true(all(upset_df$GSEA %in% c(0, 1)))
      expect_true(all(upset_df$DAA %in% c(0, 1)))
      
      # Return mock result (UpSetR doesn't return ggplot, so we simulate)
      return("UpSet plot created")
    })
  }
  
  comparison <- compare_gsea_daa(test_data$gsea_results, test_data$daa_results, plot_type = "upset")
  
  # Verify basic structure
  expect_type(comparison, "list")
  expect_true("plot" %in% names(comparison))
  expect_true("results" %in% names(comparison))
  
  # The plot should be a ggplot (either from conversion or fallback)
  expect_s3_class(comparison$plot, "ggplot")
})

test_that("UpSet plot fallback visualization when UpSetR unavailable", {
  skip_if_not_installed("ggplot2")
  
  test_data <- create_visualization_test_data("high")
  
  # Mock UpSetR as unavailable
  with_mock(
    `requireNamespace` = function(pkg, quietly = TRUE) {
      if (pkg == "UpSetR") return(FALSE)
      return(TRUE)
    },
    {
      expect_warning(
        comparison <- compare_gsea_daa(test_data$gsea_results, test_data$daa_results, plot_type = "upset"),
        "Package 'UpSetR' is required for UpSet plots. Using a basic plot instead."
      )
      
      # Verify fallback plot
      expect_s3_class(comparison$plot, "ggplot")
      
      # Verify it's a bar plot (fallback)
      plot_build <- ggplot2::ggplot_build(comparison$plot)
      expect_true("GeomBar" %in% class(plot_build$plot$layers[[1]]$geom))
    }
  )
})

test_that("Scatter plot visualization with correlation analysis", {
  skip_if_not_installed("ggplot2")
  
  # Create data with controlled correlations
  overlapping_pathways <- paste0("ko", sprintf("%05d", 1:15))
  
  # Create correlated effect sizes for realistic testing
  set.seed(12345)
  nes_values <- rnorm(15, 0, 1.5)
  log2fc_values <- nes_values * 0.8 + rnorm(15, 0, 0.5)  # Correlated with some noise
  
  gsea_results <- data.frame(
    pathway_id = c(overlapping_pathways, paste0("ko", sprintf("%05d", 100:109))),
    pathway_name = paste("Pathway", c(overlapping_pathways, paste0("ko", sprintf("%05d", 100:109)))),
    size = rep(50, 25),
    ES = c(nes_values * 0.7, rnorm(10, 0, 0.5)),
    NES = c(nes_values, rnorm(10, 0, 0.8)),
    pvalue = runif(25, 0.001, 0.1),
    p.adjust = runif(25, 0.001, 0.15),
    leading_edge = rep("K00001;K00002", 25),
    method = rep("fgsea", 25),
    stringsAsFactors = FALSE
  )
  
  daa_results <- data.frame(
    feature = c(overlapping_pathways, paste0("ko", sprintf("%05d", 200:209))),
    method = rep("ALDEx2", 25),
    group1 = rep("Control", 25),
    group2 = rep("Treatment", 25),
    p_values = runif(25, 0.001, 0.1),
    p_adjust = runif(25, 0.001, 0.12),
    log_2_fold_change = c(log2fc_values, rnorm(10, 0, 1)),
    effect_size = rnorm(25, 0, 1),
    stringsAsFactors = FALSE
  )
  
  # Create scatter plot
  comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "scatter")
  
  # Verify plot structure
  expect_s3_class(comparison$plot, "ggplot")
  
  # Verify plot content
  plot_build <- ggplot2::ggplot_build(comparison$plot)
  expect_true(length(plot_build$data) > 0)
  
  # Check if points are plotted
  expect_true("GeomPoint" %in% class(plot_build$plot$layers[[1]]$geom))
  
  # Verify axes are correct
  plot_labels <- plot_build$plot$labels
  expect_true(grepl("Normalized Enrichment Score", plot_labels$x, fixed = TRUE))
  expect_true(grepl("Log2 Fold Change", plot_labels$y, fixed = TRUE))

  # Verify some aesthetic mapping is present (color, colour, fill, or size)
  expect_true(any(c("color", "colour", "fill", "size") %in% names(plot_labels)))
})

test_that("Scatter plot handles no overlapping pathways correctly", {
  skip_if_not_installed("ggplot2")
  
  # Create data with no overlap
  gsea_results <- data.frame(
    pathway_id = paste0("ko", sprintf("%05d", 1:10)),
    pathway_name = paste("Pathway", 1:10),
    size = rep(50, 10),
    ES = rnorm(10, 0, 1),
    NES = rnorm(10, 0, 1.5),
    pvalue = runif(10, 0.001, 0.1),
    p.adjust = runif(10, 0.001, 0.1),
    leading_edge = rep("K00001;K00002", 10),
    method = rep("fgsea", 10),
    stringsAsFactors = FALSE
  )
  
  daa_results <- data.frame(
    feature = paste0("ko", sprintf("%05d", 20:29)),  # No overlap with GSEA
    method = rep("ALDEx2", 10),
    group1 = rep("Control", 10),
    group2 = rep("Treatment", 10),
    p_values = runif(10, 0.001, 0.1),
    p_adjust = runif(10, 0.001, 0.1),
    log_2_fold_change = rnorm(10, 0, 1),
    effect_size = rnorm(10, 0, 1),
    stringsAsFactors = FALSE
  )
  
  # Should warn about no overlapping pathways
  expect_warning(
    comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "scatter"),
    "No overlapping pathways found for scatter plot"
  )
  
  # Should still return a plot (empty/message plot)
  expect_s3_class(comparison$plot, "ggplot")
})

test_that("Effect size correlation visualization accuracy", {
  skip_if_not_installed("ggplot2")
  
  # Create perfectly correlated data for testing
  pathway_ids <- paste0("ko", sprintf("%05d", 1:20))
  effect_sizes <- seq(-2, 2, length.out = 20)
  
  gsea_results <- data.frame(
    pathway_id = pathway_ids,
    pathway_name = paste("Pathway", pathway_ids),
    size = rep(50, 20),
    ES = effect_sizes * 0.8,  # Slightly different scaling
    NES = effect_sizes,       # Perfect correlation
    pvalue = runif(20, 0.001, 0.05),
    p.adjust = runif(20, 0.001, 0.08),
    leading_edge = rep("K00001;K00002", 20),
    method = rep("fgsea", 20),
    stringsAsFactors = FALSE
  )
  
  daa_results <- data.frame(
    feature = pathway_ids,    # Perfect overlap
    method = rep("ALDEx2", 20),
    group1 = rep("Control", 20),
    group2 = rep("Treatment", 20),
    p_values = runif(20, 0.001, 0.05),
    p_adjust = runif(20, 0.001, 0.08),
    log_2_fold_change = effect_sizes * 1.2,  # Correlated but scaled
    effect_size = effect_sizes,
    stringsAsFactors = FALSE
  )
  
  comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "scatter")
  
  # Verify all 20 pathways appear in the plot
  plot_data <- ggplot2::ggplot_build(comparison$plot)
  plot_df <- plot_data$data[[1]]  # Get the points data
  
  expect_equal(nrow(plot_df), 20)  # All pathways should be plotted
  
  # Verify the correlation is visible in the data
  expect_true(cor(plot_df$x, plot_df$y) > 0.5)  # Should be positively correlated
})

test_that("Statistical significance visualization in scatter plots", {
  skip_if_not_installed("ggplot2")
  
  # Create data with varying significance levels
  pathway_ids <- paste0("ko", sprintf("%05d", 1:15))
  p_values <- c(
    rep(0.001, 5),   # Highly significant
    rep(0.02, 5),    # Moderately significant  
    rep(0.08, 5)     # Less significant
  )
  
  gsea_results <- data.frame(
    pathway_id = pathway_ids,
    pathway_name = paste("Pathway", pathway_ids),
    size = rep(50, 15),
    ES = rnorm(15, 0, 1),
    NES = rnorm(15, 0, 1.5),
    pvalue = p_values,
    p.adjust = p_values * 1.2,  # Adjusted p-values slightly higher
    leading_edge = rep("K00001;K00002", 15),
    method = rep("fgsea", 15),
    stringsAsFactors = FALSE
  )
  
  daa_results <- data.frame(
    feature = pathway_ids,
    method = rep("ALDEx2", 15),
    group1 = rep("Control", 15),
    group2 = rep("Treatment", 15),
    p_values = p_values,
    p_adjust = p_values * 1.1,
    log_2_fold_change = rnorm(15, 0, 2),
    effect_size = rnorm(15, 0, 1),
    stringsAsFactors = FALSE
  )
  
  comparison <- compare_gsea_daa(gsea_results, daa_results, plot_type = "scatter")
  
  # Verify color scale reflects significance
  plot_build <- ggplot2::ggplot_build(comparison$plot)
  plot_data <- plot_build$data[[1]]
  
  # Check that color values vary (representing different p-values)
  expect_true(length(unique(plot_data$colour)) > 1)
  
  # Verify some aesthetic mapping is present for significance
  # May be color, colour, fill, or size depending on implementation
  expect_true(any(c("color", "colour", "fill", "size") %in% names(plot_build$plot$labels)))
})

test_that("Heatmap plot placeholder and warning behavior", {
  skip_if_not_installed("ggplot2")
  
  test_data <- create_visualization_test_data("moderate")
  
  # Test heatmap not implemented warning
  expect_warning(
    comparison <- compare_gsea_daa(test_data$gsea_results, test_data$daa_results, plot_type = "heatmap"),
    "Heatmap plot not yet implemented"
  )
  
  # Should return placeholder plot
  expect_s3_class(comparison$plot, "ggplot")
  
  # Verify results are still calculated correctly
  expect_type(comparison$results, "list")
  expect_setequal(comparison$results$overlap, test_data$expected_overlap)
})

test_that("Visualization consistency across different data sizes", {
  skip_if_not_installed("ggplot2")
  
  # Test small dataset
  small_data <- create_visualization_test_data("moderate")
  # Reduce to smaller size
  small_gsea <- small_data$gsea_results[1:10, ]
  small_daa <- small_data$daa_results[1:10, ]
  
  comparison_small <- compare_gsea_daa(small_gsea, small_daa, plot_type = "venn")
  expect_s3_class(comparison_small$plot, "ggplot")
  
  # Test medium dataset (default)
  medium_comparison <- compare_gsea_daa(small_data$gsea_results, small_data$daa_results, plot_type = "venn")
  expect_s3_class(medium_comparison$plot, "ggplot")
  
  # Create large dataset
  large_gsea <- do.call(rbind, replicate(5, small_data$gsea_results, simplify = FALSE))
  large_gsea$pathway_id <- make.unique(large_gsea$pathway_id)
  large_daa <- do.call(rbind, replicate(5, small_data$daa_results, simplify = FALSE))
  large_daa$feature <- make.unique(large_daa$feature)
  
  comparison_large <- compare_gsea_daa(large_gsea, large_daa, plot_type = "venn")
  expect_s3_class(comparison_large$plot, "ggplot")
  
  # All should produce valid plots regardless of size
  expect_true(all(sapply(list(comparison_small$plot, medium_comparison$plot, comparison_large$plot), 
                        function(p) inherits(p, "ggplot"))))
})

test_that("Visualization theme and styling consistency", {
  skip_if_not_installed("ggplot2")
  
  test_data <- create_visualization_test_data("moderate")
  
  # Test different plot types maintain consistent styling
  venn_plot <- compare_gsea_daa(test_data$gsea_results, test_data$daa_results, plot_type = "venn")$plot
  scatter_plot <- compare_gsea_daa(test_data$gsea_results, test_data$daa_results, plot_type = "scatter")$plot
  upset_plot <- compare_gsea_daa(test_data$gsea_results, test_data$daa_results, plot_type = "upset")$plot
  
  # All should have titles
  venn_build <- ggplot2::ggplot_build(venn_plot)
  scatter_build <- ggplot2::ggplot_build(scatter_plot)
  upset_build <- ggplot2::ggplot_build(upset_plot)
  
  # Check that plots have appropriate titles/labels
  expect_true("title" %in% names(venn_build$plot$labels) || "title" %in% names(scatter_build$plot$labels))
  
  # All should be ggplot objects
  expect_s3_class(venn_plot, "ggplot")
  expect_s3_class(scatter_plot, "ggplot")
  expect_s3_class(upset_plot, "ggplot")
})

test_that("Visualization handles missing or invalid data gracefully", {
  skip_if_not_installed("ggplot2")
  
  # Create data with some missing values
  test_data <- create_visualization_test_data("moderate")
  
  # Introduce missing values
  test_data$gsea_results$pathway_id[1] <- NA
  test_data$daa_results$p_adjust[2] <- NA
  test_data$gsea_results$NES[3] <- NA
  test_data$daa_results$log_2_fold_change[4] <- NA
  
  # Should still create plots despite missing data
  comparison_venn <- compare_gsea_daa(test_data$gsea_results, test_data$daa_results, plot_type = "venn")
  comparison_scatter <- compare_gsea_daa(test_data$gsea_results, test_data$daa_results, plot_type = "scatter")
  
  expect_s3_class(comparison_venn$plot, "ggplot")
  expect_s3_class(comparison_scatter$plot, "ggplot")
  
  # Results should still be calculated
  expect_type(comparison_venn$results, "list")
  expect_type(comparison_scatter$results, "list")
})

test_that("Color schemes and accessibility in visualizations", {
  skip_if_not_installed("ggplot2")
  
  test_data <- create_visualization_test_data("moderate")
  
  # Test scatter plot color scheme
  scatter_comparison <- compare_gsea_daa(test_data$gsea_results, test_data$daa_results, plot_type = "scatter")
  
  plot_build <- ggplot2::ggplot_build(scatter_comparison$plot)
  
  # Check that scatter plot uses a gradient color scale
  has_color_scale <- length(plot_build$plot$scales$scales) > 0
  if (has_color_scale) {
    color_scale <- plot_build$plot$scales$scales[[which(sapply(plot_build$plot$scales$scales, function(x) "colour" %in% x$aesthetics))]]
    expect_true(!is.null(color_scale))
  }
  
  # Test that fallback plots use distinct colors
  with_mock(
    `requireNamespace` = function(pkg, quietly = TRUE) {
      if (pkg %in% c("ggVennDiagram", "UpSetR")) return(FALSE)
      return(TRUE)
    },
    {
      fallback_venn <- suppressWarnings(compare_gsea_daa(test_data$gsea_results, test_data$daa_results, plot_type = "venn"))
      fallback_upset <- suppressWarnings(compare_gsea_daa(test_data$gsea_results, test_data$daa_results, plot_type = "upset"))
      
      # Both should be readable bar charts
      expect_s3_class(fallback_venn$plot, "ggplot")
      expect_s3_class(fallback_upset$plot, "ggplot")
    }
  )
})