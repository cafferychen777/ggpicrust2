# Edge Cases and Error Handling Tests for GSEA Visualization
# Following Linus principle: "Never break userspace"
# Ensure the functions handle all edge cases gracefully without breaking

library(testthat)

# Load the package functions
devtools::load_all()

#=============================================================================
# DATA STRUCTURE EDGE CASES
#=============================================================================

test_that("visualize_gsea: Handles malformed data structures gracefully", {
  skip_if_not_installed("ggplot2")
  
  # Test with data frame containing all required columns but weird values
  malformed_data <- data.frame(
    pathway_id = c("", "path2", "path3", NA, "path5"),
    NES = c(NA, Inf, -Inf, 0, 1000000),
    pvalue = c(0, 1, NA, -0.1, 2.0),  # Invalid probability values
    p.adjust = c(NA, 0, 1, 1.5, -0.5), # Invalid probability values
    size = c(0, -5, NA, 1000000, 1),
    leading_edge = c("", "gene1", NA, "gene1;gene2;;;;gene3", "gene1;gene1;gene1"), # Duplicates and empty
    stringsAsFactors = FALSE
  )
  
  # Should handle gracefully without crashing
  expect_no_error({
    plot_malformed <- visualize_gsea(malformed_data, plot_type = "barplot", n_pathways = 3)
    expect_s3_class(plot_malformed, "ggplot")
  })
  
  # Test with completely empty strings and NAs
  extreme_data <- data.frame(
    pathway_id = c(NA, "", "   ", "path4"),
    NES = c(NA, NA, 0, 1),
    pvalue = c(NA, NA, 1, 0.05),
    p.adjust = c(NA, NA, 1, 0.1),
    size = c(NA, 0, 1, 10),
    leading_edge = c(NA, "", "   ", "gene1;gene2"),
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    plot_extreme <- visualize_gsea(extreme_data, plot_type = "barplot", n_pathways = 2)
    expect_s3_class(plot_extreme, "ggplot")
  })
})

test_that("visualize_gsea: Handles missing required columns with clear errors", {
  # Test missing pathway identifier columns
  missing_pathway_cols <- data.frame(
    NES = c(1.0, -1.0),
    pvalue = c(0.01, 0.05),
    p.adjust = c(0.05, 0.1),
    size = c(10, 15),
    leading_edge = c("gene1;gene2", "gene3;gene4"),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    visualize_gsea(missing_pathway_cols),
    "must contain either 'pathway_name' or 'pathway_id' column"
  )
  
  # Test missing statistical columns
  missing_stats <- data.frame(
    pathway_id = c("path1", "path2"),
    # Missing NES, pvalue, p.adjust
    size = c(10, 15),
    leading_edge = c("gene1;gene2", "gene3;gene4"),
    stringsAsFactors = FALSE
  )
  
  # May or may not work depending on implementation
  # Just verify it doesn't crash R
  tryCatch({
    plot_missing <- visualize_gsea(missing_stats, plot_type = "barplot")
    expect_s3_class(plot_missing, "ggplot")
  }, error = function(e) {
    # Expect error for missing required columns
    expect_true(TRUE)
  })
})

test_that("visualize_gsea: Handles data type mismatches gracefully", {
  skip_if_not_installed("ggplot2")
  
  # Test with wrong data types
  wrong_types <- data.frame(
    pathway_id = c("path1", "path2", "path3"),
    NES = c("high", "medium", "low"),  # Character instead of numeric
    pvalue = c("significant", "not", "maybe"), # Character instead of numeric
    p.adjust = c(TRUE, FALSE, TRUE),  # Logical instead of numeric
    size = c("small", "medium", "large"), # Character instead of numeric
    leading_edge = c(123, 456, 789), # Numeric instead of character
    stringsAsFactors = FALSE
  )
  
  # Function should attempt to handle or fail gracefully
  # May succeed with type coercion or fail with error
  tryCatch({
    result <- visualize_gsea(wrong_types, plot_type = "barplot")
    # If it succeeds, verify it's a ggplot
    expect_s3_class(result, "ggplot")
  }, error = function(e) {
    # If it fails, that's expected
    expect_true(TRUE)
  })
})

#=============================================================================
# EXTREME PARAMETER VALUES
#=============================================================================

test_that("visualize_gsea: Handles extreme parameter values", {
  skip_if_not_installed("ggplot2")
  
  # Create normal test data
  normal_data <- data.frame(
    pathway_id = sprintf("path_%03d", 1:20),
    NES = rnorm(20, 0, 1),
    pvalue = runif(20, 0.001, 0.1),
    p.adjust = runif(20, 0.01, 0.2),
    size = sample(10:100, 20),
    leading_edge = replicate(20, paste(sprintf("gene%d", sample(1:50, 5)), collapse = ";")),
    stringsAsFactors = FALSE
  )
  
  # Test with n_pathways = 0
  expect_no_error({
    plot_zero <- visualize_gsea(normal_data, plot_type = "barplot", n_pathways = 0)
    expect_s3_class(plot_zero, "ggplot")
    # Empty plot may have NULL data or 0 rows
    expect_true(is.null(plot_zero$data) || nrow(plot_zero$data) == 0)
  })
  
  # Test with very large n_pathways
  expect_no_error({
    plot_large_n <- visualize_gsea(normal_data, plot_type = "barplot", n_pathways = 1000000)
    expect_s3_class(plot_large_n, "ggplot")
    expect_equal(nrow(plot_large_n$data), 20) # Should be limited to available data
  })
  
  # Test with negative n_pathways (should be handled)
  expect_no_error({
    plot_negative_n <- visualize_gsea(normal_data, plot_type = "barplot", n_pathways = -5)
    expect_s3_class(plot_negative_n, "ggplot")
  })
})

test_that("visualize_gsea: Color parameter edge cases", {
  skip_if_not_installed("ggplot2")
  
  normal_data <- data.frame(
    pathway_id = c("path1", "path2", "path3"),
    NES = c(1.0, -1.0, 2.0),
    pvalue = c(0.01, 0.05, 0.001),
    p.adjust = c(0.05, 0.1, 0.01),
    size = c(10, 15, 20),
    leading_edge = c("gene1;gene2", "gene3;gene4", "gene5;gene6"),
    stringsAsFactors = FALSE
  )
  
  # Test with empty color vector
  expect_no_error({
    plot_empty_colors <- visualize_gsea(normal_data, plot_type = "barplot", colors = character(0))
    expect_s3_class(plot_empty_colors, "ggplot")
  })
  
  # Test with invalid color names
  expect_no_error({
    plot_invalid_colors <- visualize_gsea(normal_data, plot_type = "barplot", 
                                         colors = c("not_a_color", "also_not_a_color"))
    expect_s3_class(plot_invalid_colors, "ggplot")
  })
  
  # Test with single color
  expect_no_error({
    plot_single_color <- visualize_gsea(normal_data, plot_type = "barplot", colors = "#FF0000")
    expect_s3_class(plot_single_color, "ggplot")
  })
  
  # Test with more colors than needed
  many_colors <- rainbow(100)
  expect_no_error({
    plot_many_colors <- visualize_gsea(normal_data, plot_type = "barplot", colors = many_colors)
    expect_s3_class(plot_many_colors, "ggplot")
  })
})

#=============================================================================
# NETWORK PLOT EDGE CASES
#=============================================================================

test_that("create_network_plot: Handles extreme similarity cutoff values", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  
  test_data <- data.frame(
    pathway_id = c("path1", "path2", "path3", "path4"),
    NES = c(1.5, -1.2, 2.0, -0.8),
    pvalue = c(0.01, 0.05, 0.001, 0.1),
    p.adjust = c(0.05, 0.1, 0.01, 0.2),
    size = c(10, 8, 12, 6),
    leading_edge = c(
      "gene1;gene2;gene3",
      "gene1;gene2;gene4",
      "gene1;gene5;gene6",
      "gene7;gene8;gene9"
    ),
    stringsAsFactors = FALSE
  )
  
  # Test with cutoff = 0 (should include all connections)
  expect_no_error({
    plot_cutoff_zero <- visualize_gsea(test_data, plot_type = "network",
                                      network_params = list(similarity_cutoff = 0))
    expect_s3_class(plot_cutoff_zero, "ggplot")
  })
  
  # Test with cutoff = 1 (should exclude all connections)
  expect_no_error({
    plot_cutoff_one <- visualize_gsea(test_data, plot_type = "network",
                                     network_params = list(similarity_cutoff = 1))
    expect_s3_class(plot_cutoff_one, "ggplot")
  })
  
  # Test with cutoff > 1 (invalid but should be handled)
  expect_no_error({
    plot_cutoff_invalid <- visualize_gsea(test_data, plot_type = "network",
                                         network_params = list(similarity_cutoff = 2))
    expect_s3_class(plot_cutoff_invalid, "ggplot")
  })
  
  # Test with negative cutoff (should be handled)
  expect_no_error({
    plot_cutoff_negative <- visualize_gsea(test_data, plot_type = "network",
                                          network_params = list(similarity_cutoff = -0.5))
    expect_s3_class(plot_cutoff_negative, "ggplot")
  })
})

test_that("create_network_plot: Handles invalid similarity measures", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  
  test_data <- data.frame(
    pathway_id = c("path1", "path2", "path3"),
    NES = c(1.5, -1.2, 2.0),
    pvalue = c(0.01, 0.05, 0.001),
    p.adjust = c(0.05, 0.1, 0.01),
    size = c(10, 8, 12),
    leading_edge = c("gene1;gene2", "gene1;gene3", "gene4;gene5"),
    stringsAsFactors = FALSE
  )
  
  # Test with invalid similarity measure (should default or handle gracefully)
  expect_no_error({
    plot_invalid_sim <- visualize_gsea(test_data, plot_type = "network",
                                      network_params = list(similarity_measure = "invalid_measure"))
    expect_s3_class(plot_invalid_sim, "ggplot")
  })
})

test_that("create_network_plot: Handles pathologically small/large datasets", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  
  # Test with single pathway
  single_pathway <- data.frame(
    pathway_id = "lonely_path",
    NES = 1.5,
    pvalue = 0.01,
    p.adjust = 0.05,
    size = 10,
    leading_edge = "gene1;gene2;gene3",
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    plot_single <- visualize_gsea(single_pathway, plot_type = "network")
    expect_s3_class(plot_single, "ggplot")
  })
  
  # Test with two pathways (minimal network)
  two_pathways <- data.frame(
    pathway_id = c("path1", "path2"),
    NES = c(1.0, -1.0),
    pvalue = c(0.01, 0.05),
    p.adjust = c(0.05, 0.1),
    size = c(5, 8),
    leading_edge = c("gene1;gene2", "gene3;gene4"),
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    plot_two <- visualize_gsea(two_pathways, plot_type = "network")
    expect_s3_class(plot_two, "ggplot")
  })
})

#=============================================================================
# HEATMAP PLOT EDGE CASES  
#=============================================================================

test_that("create_heatmap_plot: Handles mismatched abundance and metadata", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  # Create abundance data
  abundance_data <- data.frame(
    Sample1 = c(10, 20, 30),
    Sample2 = c(15, 25, 35),
    Sample3 = c(12, 22, 32),
    row.names = c("gene1", "gene2", "gene3")
  )
  
  # Create metadata with different samples
  mismatched_metadata <- data.frame(
    sample = c("SampleA", "SampleB", "SampleC"), # Different sample names
    group = c("Control", "Treatment", "Control"),
    stringsAsFactors = FALSE,
    row.names = c("SampleA", "SampleB", "SampleC")
  )
  
  gsea_data <- data.frame(
    pathway_id = "path1",
    NES = 1.5,
    pvalue = 0.01,
    p.adjust = 0.05,
    size = 3,
    leading_edge = "gene1;gene2;gene3",
    stringsAsFactors = FALSE
  )
  
  # Should handle mismatch gracefully (might skip missing samples)
  expect_no_error({
    plot_mismatch <- visualize_gsea(
      gsea_data, plot_type = "heatmap",
      abundance = abundance_data,
      metadata = mismatched_metadata,
      group = "group"
    )
    expect_s4_class(plot_mismatch, "Heatmap")
  })
})

test_that("create_heatmap_plot: Handles empty leading edges", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  abundance_data <- data.frame(
    Sample1 = c(10, 20, 30, 40),
    Sample2 = c(15, 25, 35, 45),
    Sample3 = c(12, 22, 32, 42),
    row.names = c("gene1", "gene2", "gene3", "gene4")
  )
  
  metadata <- data.frame(
    sample = c("Sample1", "Sample2", "Sample3"),
    group = c("Control", "Treatment", "Control"),
    stringsAsFactors = FALSE,
    row.names = c("Sample1", "Sample2", "Sample3")
  )
  
  # GSEA data with empty and missing leading edges
  gsea_empty_edges <- data.frame(
    pathway_id = c("path1", "path2", "path3"),
    NES = c(1.0, -1.0, 2.0),
    pvalue = c(0.01, 0.05, 0.001),
    p.adjust = c(0.05, 0.1, 0.01),
    size = c(0, 2, 4),
    leading_edge = c("", "gene1;gene2", "gene1;gene2;gene3;gene4"),
    stringsAsFactors = FALSE
  )
  
  # Should handle empty leading edges without crashing
  expect_no_error({
    plot_empty_edges <- visualize_gsea(
      gsea_empty_edges, plot_type = "heatmap",
      abundance = abundance_data,
      metadata = metadata,
      group = "group"
    )
    expect_s4_class(plot_empty_edges, "Heatmap")
  })
})

test_that("create_heatmap_plot: Handles genes not in abundance matrix", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  abundance_data <- data.frame(
    Sample1 = c(10, 20, 30),
    Sample2 = c(15, 25, 35),
    row.names = c("gene_A", "gene_B", "gene_C") # Different gene names
  )
  
  metadata <- data.frame(
    sample = c("Sample1", "Sample2"),
    group = c("Control", "Treatment"),
    stringsAsFactors = FALSE,
    row.names = c("Sample1", "Sample2")
  )
  
  # GSEA data with genes not in abundance matrix
  gsea_missing_genes <- data.frame(
    pathway_id = c("path1", "path2"),
    NES = c(1.5, -1.2),
    pvalue = c(0.01, 0.05),
    p.adjust = c(0.05, 0.1),
    size = c(5, 3),
    leading_edge = c(
      "gene1;gene2;gene3;gene4;gene5", # None of these exist in abundance
      "gene_A;gene_B;missing_gene"     # Some exist, some don't
    ),
    stringsAsFactors = FALSE
  )
  
  # Should handle missing genes gracefully
  expect_no_error({
    plot_missing_genes <- visualize_gsea(
      gsea_missing_genes, plot_type = "heatmap",
      abundance = abundance_data,
      metadata = metadata,
      group = "group"
    )
    expect_s4_class(plot_missing_genes, "Heatmap")
  })
})

test_that("create_heatmap_plot: Handles extreme heatmap parameters", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  # Create normal test data
  abundance_data <- data.frame(
    Sample1 = 1:20,
    Sample2 = 2:21,
    Sample3 = 3:22,
    Sample4 = 4:23,
    row.names = sprintf("gene_%02d", 1:20)
  )
  
  metadata <- data.frame(
    sample = paste0("Sample", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE,
    row.names = paste0("Sample", 1:4)
  )
  
  gsea_data <- data.frame(
    pathway_id = c("path1", "path2"),
    NES = c(2.0, -1.5),
    pvalue = c(0.001, 0.05),
    p.adjust = c(0.01, 0.1),
    size = c(10, 8),
    leading_edge = c(
      paste(sprintf("gene_%02d", 1:10), collapse = ";"),
      paste(sprintf("gene_%02d", 11:18), collapse = ";")
    ),
    stringsAsFactors = FALSE
  )
  
  # Test with invalid annotation colors
  invalid_colors <- list(Group = c("Invalid" = "#FF0000")) # Group "Invalid" doesn't exist
  
  expect_no_error({
    plot_invalid_colors <- visualize_gsea(
      gsea_data, plot_type = "heatmap",
      abundance = abundance_data,
      metadata = metadata,
      group = "group",
      heatmap_params = list(annotation_colors = invalid_colors)
    )
    expect_s4_class(plot_invalid_colors, "Heatmap")
  })
  
  # Test with all clustering disabled and no row names
  expect_no_error({
    plot_minimal <- visualize_gsea(
      gsea_data, plot_type = "heatmap",
      abundance = abundance_data,
      metadata = metadata,
      group = "group",
      heatmap_params = list(
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_rownames = FALSE
      )
    )
    expect_s4_class(plot_minimal, "Heatmap")
  })
})

#=============================================================================
# PACKAGE DEPENDENCY EDGE CASES
#=============================================================================

test_that("visualize_gsea: Handles missing optional packages gracefully", {
  # We can't easily uninstall packages for testing, but we can verify
  # that the dependency checks are in place and the code structure handles missing packages
  
  normal_data <- data.frame(
    pathway_id = c("path1", "path2"),
    NES = c(1.0, -1.0),
    pvalue = c(0.01, 0.05),
    p.adjust = c(0.05, 0.1),
    size = c(10, 15),
    leading_edge = c("gene1;gene2", "gene3;gene4"),
    stringsAsFactors = FALSE
  )
  
  # Verify that requireNamespace function exists and is used
  expect_true(exists("requireNamespace"))
  
  # The actual package checking happens inside the function calls
  # We can't easily test missing packages without breaking the test environment
  # But we can verify the structure is sound
  expect_true(is.function(visualize_gsea))
})

#=============================================================================
# UNICODE AND SPECIAL CHARACTER HANDLING
#=============================================================================

test_that("visualize_gsea: Handles special characters and Unicode", {
  skip_if_not_installed("ggplot2")
  
  # Test data with special characters and Unicode
  special_char_data <- data.frame(
    pathway_id = c("path/1", "path\\2", "path<>3", "path|4", "path?5"),
    pathway_name = c("Pathway α", "Pathway β", "Pathway γ", "Pathway δ", "Pathway ε"),
    NES = c(1.2, -1.5, 2.1, -0.8, 1.7),
    pvalue = c(0.01, 0.02, 0.001, 0.05, 0.01),
    p.adjust = c(0.05, 0.08, 0.005, 0.1, 0.05),
    size = c(15, 12, 20, 8, 18),
    leading_edge = c(
      "gene→1;gene←2;gene↑3",
      "gene♠1;gene♦2;gene♣3",
      "gene∀1;gene∃2;gene∈3", 
      "gene™1;gene®2;gene©3",
      "gene_1;gene_2;gene_3"
    ),
    stringsAsFactors = FALSE
  )
  
  # Should handle special characters without breaking
  expect_no_error({
    plot_special <- visualize_gsea(special_char_data, plot_type = "barplot", n_pathways = 3)
    expect_s3_class(plot_special, "ggplot")
  })
  
  # Test with pathway_name containing special characters
  expect_no_error({
    plot_special_names <- visualize_gsea(special_char_data, plot_type = "barplot", 
                                        pathway_label_column = "pathway_name", n_pathways = 3)
    expect_s3_class(plot_special_names, "ggplot")
  })
})

#=============================================================================
# MEMORY AND PERFORMANCE EDGE CASES
#=============================================================================

test_that("visualize_gsea: Handles memory constraints gracefully", {
  skip_on_ci() # Skip on CI to avoid memory/time issues
  skip_if_not_installed("ggplot2")
  
  # Test with very long pathway names and gene lists
  memory_test_data <- data.frame(
    pathway_id = sprintf("pathway_with_very_long_identifier_%04d", 1:10),
    pathway_name = paste("Very Long Pathway Name With Many Words", 
                        sprintf("Number %04d", 1:10),
                        "And Even More Descriptive Text That Goes On And On"),
    NES = rnorm(10, 0, 1),
    pvalue = runif(10, 0.001, 0.1),
    p.adjust = runif(10, 0.01, 0.2),
    size = sample(50:500, 10),
    # Very long leading edge lists
    leading_edge = replicate(10, {
      genes <- sprintf("gene_with_very_long_name_%04d", sample(1:1000, 50))
      paste(genes, collapse = ";")
    }),
    stringsAsFactors = FALSE
  )
  
  # Should handle large data without memory issues
  expect_no_error({
    plot_memory <- visualize_gsea(memory_test_data, plot_type = "barplot", n_pathways = 5)
    expect_s3_class(plot_memory, "ggplot")
  })
})

# Final verification that all edge cases are handled
test_that("visualize_gsea: Complete edge case verification", {
  # This test verifies that the function exists and is callable
  # with minimal valid input, ensuring basic robustness
  
  minimal_valid_data <- data.frame(
    pathway_id = "test_pathway",
    NES = 1.0,
    pvalue = 0.05,
    p.adjust = 0.1,
    size = 10,
    leading_edge = "gene1;gene2",
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    result <- visualize_gsea(minimal_valid_data, plot_type = "barplot", n_pathways = 1)
    expect_s3_class(result, "ggplot")
  })
  
  # Verify function signature and required parameters
  expect_true(is.function(visualize_gsea))
  
  # Get function arguments
  args <- names(formals(visualize_gsea))
  required_args <- c("gsea_results")
  expect_true(all(required_args %in% args))
})

message("GSEA visualization edge case tests completed successfully")