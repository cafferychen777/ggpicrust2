test_that("visualize_gsea handles pathway labels correctly", {
  # Create mock GSEA results without pathway_name
  gsea_results_no_names <- data.frame(
    pathway_id = c("ko00010", "ko00020", "ko00030"),
    NES = c(2.5, -1.8, 3.2),
    pvalue = c(0.001, 0.05, 0.0001),
    p.adjust = c(0.01, 0.1, 0.001),
    size = c(50, 30, 80),
    leading_edge = c("gene1;gene2", "gene3;gene4", "gene5;gene6"),
    stringsAsFactors = FALSE
  )
  
  # Create mock GSEA results with pathway_name
  gsea_results_with_names <- gsea_results_no_names
  gsea_results_with_names$pathway_name <- c(
    "Glycolysis / Gluconeogenesis",
    "Citrate cycle (TCA cycle)",
    "Pentose phosphate pathway"
  )
  
  # Test 1: Function should use pathway_id when pathway_name is not available
  expect_no_error({
    p1 <- visualize_gsea(gsea_results_no_names, plot_type = "barplot", n_pathways = 3)
  })
  
  # Test 2: Function should automatically use pathway_name when available
  expect_no_error({
    p2 <- visualize_gsea(gsea_results_with_names, plot_type = "barplot", n_pathways = 3)
  })
  
  # Test 3: Function should use custom pathway_label_column when specified
  expect_no_error({
    p3 <- visualize_gsea(gsea_results_with_names, plot_type = "barplot", 
                        pathway_label_column = "pathway_name", n_pathways = 3)
  })
  
  # Test 4: Function should use pathway_id when custom column is specified as pathway_id
  expect_no_error({
    p4 <- visualize_gsea(gsea_results_with_names, plot_type = "barplot", 
                        pathway_label_column = "pathway_id", n_pathways = 3)
  })
  
  # Test 5: Function should error when custom column doesn't exist
  expect_error({
    visualize_gsea(gsea_results_with_names, plot_type = "barplot", 
                  pathway_label_column = "nonexistent_column", n_pathways = 3)
  }, "not found in gsea_results")
  
  # Test 6: Function should error when neither pathway_name nor pathway_id exist
  gsea_results_invalid <- gsea_results_no_names
  names(gsea_results_invalid)[names(gsea_results_invalid) == "pathway_id"] <- "invalid_column"
  
  expect_error({
    visualize_gsea(gsea_results_invalid, plot_type = "barplot", n_pathways = 3)
  }, "must contain either 'pathway_name' or 'pathway_id' column")
})

test_that("visualize_gsea pathway labels work across different plot types", {
  # Create mock GSEA results with pathway_name
  gsea_results <- data.frame(
    pathway_id = c("ko00010", "ko00020", "ko00030", "ko00040", "ko00050"),
    pathway_name = c(
      "Glycolysis / Gluconeogenesis",
      "Citrate cycle (TCA cycle)",
      "Pentose phosphate pathway",
      "Pentose and glucuronate interconversions",
      "Fructose and mannose metabolism"
    ),
    NES = c(2.5, -1.8, 3.2, -2.1, 1.9),
    pvalue = c(0.001, 0.05, 0.0001, 0.02, 0.03),
    p.adjust = c(0.01, 0.1, 0.001, 0.05, 0.06),
    size = c(50, 30, 80, 40, 60),
    leading_edge = c("gene1;gene2", "gene3;gene4", "gene5;gene6", "gene7;gene8", "gene9;gene10"),
    stringsAsFactors = FALSE
  )
  
  # Test different plot types with pathway names
  expect_no_error({
    p_enrichment <- visualize_gsea(gsea_results, plot_type = "enrichment_plot", n_pathways = 3)
  })
  
  expect_no_error({
    p_dotplot <- visualize_gsea(gsea_results, plot_type = "dotplot", n_pathways = 3)
  })
  
  expect_no_error({
    p_barplot <- visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 3)
  })
  
  expect_no_error({
    p_network <- visualize_gsea(gsea_results, plot_type = "network", n_pathways = 3)
  })
})

test_that("pathway_label_column parameter validation works correctly", {
  gsea_results <- data.frame(
    pathway_id = c("ko00010", "ko00020"),
    pathway_name = c("Pathway 1", "Pathway 2"),
    NES = c(2.5, -1.8),
    pvalue = c(0.001, 0.05),
    p.adjust = c(0.01, 0.1),
    size = c(50, 30),
    leading_edge = c("gene1;gene2", "gene3;gene4"),
    stringsAsFactors = FALSE
  )
  
  # Test invalid pathway_label_column type
  expect_error({
    visualize_gsea(gsea_results, plot_type = "barplot", 
                  pathway_label_column = 123, n_pathways = 2)
  }, "pathway_label_column must be NULL or a character string")
  
  # Test valid pathway_label_column
  expect_no_error({
    visualize_gsea(gsea_results, plot_type = "barplot", 
                  pathway_label_column = "pathway_name", n_pathways = 2)
  })
})
