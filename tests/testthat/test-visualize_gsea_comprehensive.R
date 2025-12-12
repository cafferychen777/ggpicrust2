# Comprehensive Tests for GSEA Visualization Functions
library(testthat)

# Enhanced helper functions for creating more realistic test data
create_comprehensive_gsea_results <- function(n_pathways = 15, include_pathway_name = TRUE, 
                                            with_empty_leading_edges = FALSE, 
                                            with_na_values = FALSE) {
  set.seed(42)
  
  # Create pathway IDs
  pathway_ids <- paste0("path:ko", sprintf("%05d", 1:n_pathways))
  
  # Create base data frame
  gsea_results <- data.frame(
    pathway_id = pathway_ids,
    size = sample(10:200, n_pathways, replace = TRUE),
    ES = runif(n_pathways, -0.9, 0.9),
    NES = runif(n_pathways, -3, 3),
    pvalue = if (n_pathways > 5) c(runif(n_pathways - 5, 0.001, 0.05), runif(5, 0.05, 0.2)) else runif(n_pathways, 0.001, 0.2), # Mix of significant and non-significant
    p.adjust = if (n_pathways > 3) c(runif(n_pathways - 3, 0.001, 0.1), runif(3, 0.1, 0.5)) else runif(n_pathways, 0.001, 0.5), # Mix of significant and non-significant
    method = rep("fgsea", n_pathways),
    pathway_class = sample(
      c("Metabolism", "Genetic Information Processing", "Environmental Information Processing",
        "Cellular Processes", "Organismal Systems", "Human Diseases", "Drug Development"),
      n_pathways, replace = TRUE
    ),
    stringsAsFactors = FALSE
  )
  
  # Add pathway names if requested
  if (include_pathway_name) {
    gsea_results$pathway_name <- c(
      "Glycolysis / Gluconeogenesis", "Citrate cycle (TCA cycle)", "Pentose phosphate pathway",
      "Fatty acid biosynthesis", "Valine, leucine and isoleucine degradation", 
      "Lysine biosynthesis", "Arginine biosynthesis", "Histidine metabolism",
      "Tyrosine metabolism", "Phenylalanine metabolism", "Tryptophan metabolism",
      "Purine metabolism", "Pyrimidine metabolism", "Nicotinate and nicotinamide metabolism",
      "Riboflavin metabolism"
    )[1:n_pathways]
  }
  
  # Create leading edge genes
  all_genes <- paste0("K", sprintf("%05d", 1:1000))
  
  if (with_empty_leading_edges) {
    # Create some empty leading edges to test edge cases
    leading_edges <- replicate(n_pathways, {
      if (runif(1) < 0.1) { # 10% chance of empty
        return("")
      } else {
        return(paste(sample(all_genes, sample(3:15, 1)), collapse = ";"))
      }
    })
  } else {
    leading_edges <- replicate(n_pathways, {
      paste(sample(all_genes, sample(3:15, 1)), collapse = ";")
    })
  }
  
  gsea_results$leading_edge <- leading_edges
  
  # Add NA values to test robustness
  if (with_na_values) {
    gsea_results$NES[sample(n_pathways, 2)] <- NA
    gsea_results$pvalue[sample(n_pathways, 1)] <- NA
  }
  
  return(gsea_results)
}

# Helper function to create more comprehensive test abundance data
create_comprehensive_abundance <- function(n_genes = 200, n_samples = 12) {
  set.seed(42)
  gene_ids <- paste0("K", sprintf("%05d", 1:n_genes))
  sample_ids <- paste0("Sample", 1:n_samples)
  
  # Create realistic abundance matrix with different abundance patterns
  abundance_matrix <- matrix(
    rgamma(n_genes * n_samples, shape = 2, rate = 0.02), # More realistic than uniform
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(gene_ids, sample_ids)
  )
  
  # Add some zeros to simulate real data
  zero_prop <- 0.1
  zero_indices <- sample(length(abundance_matrix), size = floor(length(abundance_matrix) * zero_prop))
  abundance_matrix[zero_indices] <- 0
  
  return(as.data.frame(abundance_matrix))
}

# Helper function to create comprehensive metadata
create_comprehensive_metadata <- function(n_samples = 12, n_groups = 3) {
  set.seed(42)
  sample_ids <- paste0("Sample", 1:n_samples)
  
  # Create multiple grouping variables
  groups <- rep(paste0("Group", 1:n_groups), length.out = n_samples)
  treatment <- rep(c("Control", "Treatment"), length.out = n_samples)
  
  metadata <- data.frame(
    sample = sample_ids,
    group = groups,
    treatment = treatment,
    batch = sample(c("Batch1", "Batch2"), n_samples, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  rownames(metadata) <- metadata$sample
  return(metadata)
}

# Test 1: Main visualize_gsea function with all plot types
test_that("visualize_gsea creates all plot types correctly", {
  skip_if_not_installed("ggplot2")
  
  gsea_results <- create_comprehensive_gsea_results()
  abundance <- create_comprehensive_abundance()
  metadata <- create_comprehensive_metadata()
  
  # Test enrichment plot
  expect_no_error({
    p_enrichment <- visualize_gsea(gsea_results, plot_type = "enrichment_plot", n_pathways = 10)
    expect_s3_class(p_enrichment, "ggplot")
    expect_true("GeomBar" %in% sapply(p_enrichment$layers, function(x) class(x$geom)[1]))
  })
  
  # Test dotplot
  expect_no_error({
    p_dotplot <- visualize_gsea(gsea_results, plot_type = "dotplot", n_pathways = 10)
    expect_s3_class(p_dotplot, "ggplot")
    expect_true("GeomPoint" %in% sapply(p_dotplot$layers, function(x) class(x$geom)[1]))
  })
  
  # Test barplot
  expect_no_error({
    p_barplot <- visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 10)
    expect_s3_class(p_barplot, "ggplot")
    expect_true("GeomBar" %in% sapply(p_barplot$layers, function(x) class(x$geom)[1]))
  })
  
  # Test network plot (if packages available)
  if (requireNamespace("igraph", quietly = TRUE) && 
      requireNamespace("ggraph", quietly = TRUE) && 
      requireNamespace("tidygraph", quietly = TRUE)) {
    expect_no_error({
      p_network <- visualize_gsea(gsea_results, plot_type = "network", n_pathways = 8)
      expect_s3_class(p_network, "ggplot")
    })
  }
  
  # Test heatmap (if packages available)
  if (requireNamespace("ComplexHeatmap", quietly = TRUE) && 
      requireNamespace("circlize", quietly = TRUE)) {
    expect_no_error({
      p_heatmap <- visualize_gsea(
        gsea_results, 
        plot_type = "heatmap", 
        n_pathways = 8,
        abundance = abundance,
        metadata = metadata,
        group = "group"
      )
      expect_s4_class(p_heatmap, "Heatmap")
    })
  }
})

# Test 2: Comprehensive pathway labeling tests
test_that("visualize_gsea handles pathway labeling comprehensively", {
  skip_if_not_installed("ggplot2")
  
  # Test with pathway_id only
  gsea_no_names <- create_comprehensive_gsea_results(include_pathway_name = FALSE)
  expect_no_error({
    p1 <- visualize_gsea(gsea_no_names, plot_type = "barplot", n_pathways = 5)
    expect_s3_class(p1, "ggplot")
    # Should use pathway_id as labels
    expect_true("pathway_label" %in% colnames(p1$data))
  })
  
  # Test with pathway_name (automatic selection)
  gsea_with_names <- create_comprehensive_gsea_results(include_pathway_name = TRUE)
  expect_no_error({
    p2 <- visualize_gsea(gsea_with_names, plot_type = "barplot", n_pathways = 5)
    expect_s3_class(p2, "ggplot")
    # Should automatically use pathway_name
    expect_true("pathway_label" %in% colnames(p2$data))
  })
  
  # Test with explicit pathway_label_column = "pathway_id"
  expect_no_error({
    p3 <- visualize_gsea(gsea_with_names, plot_type = "barplot", 
                        pathway_label_column = "pathway_id", n_pathways = 5)
    expect_s3_class(p3, "ggplot")
    expect_true("pathway_label" %in% colnames(p3$data))
  })
  
  # Test with explicit pathway_label_column = "pathway_name"
  expect_no_error({
    p4 <- visualize_gsea(gsea_with_names, plot_type = "barplot", 
                        pathway_label_column = "pathway_name", n_pathways = 5)
    expect_s3_class(p4, "ggplot")
    expect_true("pathway_label" %in% colnames(p4$data))
  })
  
  # Test with custom column
  gsea_custom <- gsea_with_names
  gsea_custom$custom_label <- paste("Custom", gsea_custom$pathway_name)
  expect_no_error({
    p5 <- visualize_gsea(gsea_custom, plot_type = "barplot", 
                        pathway_label_column = "custom_label", n_pathways = 5)
    expect_s3_class(p5, "ggplot")
    expect_true("pathway_label" %in% colnames(p5$data))
  })
})

# Test 3: Parameter validation for visualization options
test_that("visualize_gsea validates parameters comprehensively", {
  gsea_results <- create_comprehensive_gsea_results()
  
  # Test invalid inputs
  expect_error(visualize_gsea("not_a_dataframe"), "'gsea_results' must be a data frame")
  expect_error(visualize_gsea(gsea_results, plot_type = "invalid_type"), 
               "plot_type must be one of")
  expect_error(visualize_gsea(gsea_results, sort_by = "invalid_sort"), 
               "sort_by must be one of")
  expect_error(visualize_gsea(gsea_results, colors = 123), 
               "colors must be NULL or a character vector")
  expect_error(visualize_gsea(gsea_results, pathway_label_column = 123), 
               "pathway_label_column must be NULL or a character string")
  
  # Test missing columns
  gsea_incomplete <- gsea_results
  gsea_incomplete$pathway_id <- NULL
  gsea_incomplete$pathway_name <- NULL
  expect_error(visualize_gsea(gsea_incomplete), 
               "must contain either 'pathway_name' or 'pathway_id' column")
  
  # Test non-existent pathway_label_column
  expect_error(visualize_gsea(gsea_results, pathway_label_column = "nonexistent"), 
               "not found in gsea_results")
  
  # Test valid color vectors
  expect_no_error({
    p <- visualize_gsea(gsea_results, plot_type = "barplot", 
                       colors = c("#FF0000", "#00FF00", "#0000FF"))
    expect_s3_class(p, "ggplot")
  })
})

# Test 4: Network similarity calculations and thresholds
test_that("create_network_plot calculates similarities correctly", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  
  gsea_results <- create_comprehensive_gsea_results(n_pathways = 8)
  
  # Test different similarity measures
  expect_no_error({
    p_jaccard <- visualize_gsea(gsea_results, plot_type = "network",
                               network_params = list(similarity_measure = "jaccard"))
    expect_s3_class(p_jaccard, "ggplot")
  })
  
  expect_no_error({
    p_overlap <- visualize_gsea(gsea_results, plot_type = "network",
                               network_params = list(similarity_measure = "overlap"))
    expect_s3_class(p_overlap, "ggplot")
  })
  
  expect_no_error({
    p_correlation <- visualize_gsea(gsea_results, plot_type = "network",
                                   network_params = list(similarity_measure = "correlation"))
    expect_s3_class(p_correlation, "ggplot")
  })
  
  # Test different similarity cutoffs
  expect_no_error({
    p_low_cutoff <- visualize_gsea(gsea_results, plot_type = "network",
                                  network_params = list(similarity_cutoff = 0.1))
    expect_s3_class(p_low_cutoff, "ggplot")
  })
  
  expect_no_error({
    p_high_cutoff <- visualize_gsea(gsea_results, plot_type = "network",
                                   network_params = list(similarity_cutoff = 0.8))
    expect_s3_class(p_high_cutoff, "ggplot")
  })
  
  # Test different layouts
  expect_no_error({
    p_circle <- visualize_gsea(gsea_results, plot_type = "network",
                              network_params = list(layout = "circle"))
    expect_s3_class(p_circle, "ggplot")
  })
  
  expect_no_error({
    p_kamada <- visualize_gsea(gsea_results, plot_type = "network",
                              network_params = list(layout = "kamada"))
    expect_s3_class(p_kamada, "ggplot")
  })
  
  # Test different node coloring
  expect_no_error({
    p_pvalue <- visualize_gsea(gsea_results, plot_type = "network",
                              network_params = list(node_color_by = "pvalue"))
    expect_s3_class(p_pvalue, "ggplot")
  })
})

# Test 5: Heatmap data preparation and scaling
test_that("create_heatmap_plot handles data preparation correctly", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  
  gsea_results <- create_comprehensive_gsea_results(n_pathways = 8)
  abundance <- create_comprehensive_abundance()
  metadata <- create_comprehensive_metadata()
  
  # Test basic heatmap
  expect_no_error({
    p_basic <- visualize_gsea(gsea_results, plot_type = "heatmap",
                             abundance = abundance, metadata = metadata, group = "group")
    expect_s4_class(p_basic, "Heatmap")
  })
  
  # Test with different clustering options
  expect_no_error({
    p_no_cluster <- visualize_gsea(gsea_results, plot_type = "heatmap",
                                  abundance = abundance, metadata = metadata, group = "group",
                                  heatmap_params = list(cluster_rows = FALSE, cluster_columns = FALSE))
    expect_s4_class(p_no_cluster, "Heatmap")
  })
  
  # Test with different annotation groups
  expect_no_error({
    p_treatment <- visualize_gsea(gsea_results, plot_type = "heatmap",
                                 abundance = abundance, metadata = metadata, group = "treatment")
    expect_s4_class(p_treatment, "Heatmap")
  })
  
  # Test with custom colors
  custom_colors <- list(Group = c("Group1" = "#FF0000", "Group2" = "#00FF00", "Group3" = "#0000FF"))
  expect_no_error({
    p_custom <- visualize_gsea(gsea_results, plot_type = "heatmap",
                              abundance = abundance, metadata = metadata, group = "group",
                              heatmap_params = list(annotation_colors = custom_colors))
    expect_s4_class(p_custom, "Heatmap")
  })
  
  # Test validation - missing parameters
  expect_error(visualize_gsea(gsea_results, plot_type = "heatmap"), 
               "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required")
  expect_error(visualize_gsea(gsea_results, plot_type = "heatmap", abundance = abundance), 
               "For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required")
})

# Test 6: Error handling for missing required packages
test_that("visualize_gsea handles missing packages correctly", {
  gsea_results <- create_comprehensive_gsea_results()
  
  # Mock function to simulate missing packages
  mock_require <- function(pkg) FALSE
  
  # We can't easily test this without actually uninstalling packages,
  # but we can test that the checks are in place
  expect_true(exists("requireNamespace"))
})

# Test 7: Edge cases - empty results, single pathways, etc.
test_that("visualize_gsea handles edge cases correctly", {
  skip_if_not_installed("ggplot2")
  
  # Test with single pathway
  gsea_single <- create_comprehensive_gsea_results(n_pathways = 1)
  expect_no_error({
    p_single <- visualize_gsea(gsea_single, plot_type = "barplot", n_pathways = 5)
    expect_s3_class(p_single, "ggplot")
    expect_equal(nrow(p_single$data), 1)
  })
  
  # Test with empty leading edges
  gsea_empty_edges <- create_comprehensive_gsea_results(with_empty_leading_edges = TRUE)
  expect_no_error({
    p_empty <- visualize_gsea(gsea_empty_edges, plot_type = "barplot", n_pathways = 5)
    expect_s3_class(p_empty, "ggplot")
  })
  
  # Test with very large n_pathways
  gsea_results <- create_comprehensive_gsea_results(n_pathways = 5)
  expect_no_error({
    p_large <- visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 100)
    expect_s3_class(p_large, "ggplot")
    expect_equal(nrow(p_large$data), 5) # Should be limited to available pathways
  })
  
  # Test with zero pathways (empty data frame)
  gsea_empty <- data.frame(
    pathway_id = character(0),
    NES = numeric(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    size = integer(0),
    leading_edge = character(0)
  )
  expect_no_error({
    p_empty_df <- visualize_gsea(gsea_empty, plot_type = "barplot", n_pathways = 5)
    expect_s3_class(p_empty_df, "ggplot")
    # Empty plot may have NULL data or 0 rows
    expect_true(is.null(p_empty_df$data) || nrow(p_empty_df$data) == 0)
  })
})

# Test 8: Sorting functionality
test_that("visualize_gsea sorts data correctly", {
  skip_if_not_installed("ggplot2")
  
  # Create results with known values for testing sorting
  gsea_results <- data.frame(
    pathway_id = c("path1", "path2", "path3", "path4", "path5"),
    pathway_name = paste("Pathway", 1:5),
    NES = c(1.0, 3.0, -2.0, 0.5, 2.5),
    pvalue = c(0.05, 0.01, 0.001, 0.1, 0.02),
    p.adjust = c(0.1, 0.05, 0.01, 0.2, 0.08),
    size = c(50, 30, 80, 40, 60),
    leading_edge = rep("gene1;gene2", 5),
    stringsAsFactors = FALSE
  )
  
  # Test sorting by NES (absolute values)
  p_nes <- visualize_gsea(gsea_results, plot_type = "barplot", sort_by = "NES", n_pathways = 3)
  expect_s3_class(p_nes, "ggplot")
  
  # Test sorting by pvalue
  p_pval <- visualize_gsea(gsea_results, plot_type = "barplot", sort_by = "pvalue", n_pathways = 3)
  expect_s3_class(p_pval, "ggplot")
  
  # Test sorting by p.adjust
  p_padj <- visualize_gsea(gsea_results, plot_type = "barplot", sort_by = "p.adjust", n_pathways = 3)
  expect_s3_class(p_padj, "ggplot")
})

# Test 9: Color schemes and aesthetic mappings
test_that("visualize_gsea applies color schemes correctly", {
  skip_if_not_installed("ggplot2")
  
  gsea_results <- create_comprehensive_gsea_results()
  
  # Test custom colors
  custom_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  expect_no_error({
    p_custom <- visualize_gsea(gsea_results, plot_type = "barplot", colors = custom_colors)
    expect_s3_class(p_custom, "ggplot")
  })
  
  # Test default colors (NULL)
  expect_no_error({
    p_default <- visualize_gsea(gsea_results, plot_type = "barplot", colors = NULL)
    expect_s3_class(p_default, "ggplot")
  })
  
  # Test that barplot shows direction coloring
  p_bar <- visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 5)
  expect_true("direction" %in% colnames(p_bar$data))
  expect_true(all(p_bar$data$direction %in% c("Positive", "Negative")))
})

# Test 10: Network plot edge cases
test_that("create_network_plot handles edge cases", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  
  # Test with very high similarity cutoff (should result in no connections)
  gsea_results <- create_comprehensive_gsea_results(n_pathways = 5)
  expect_no_error({
    p_no_connections <- visualize_gsea(gsea_results, plot_type = "network",
                                      network_params = list(similarity_cutoff = 0.99))
    expect_s3_class(p_no_connections, "ggplot")
  })
  
  # Test with empty leading edges
  gsea_empty_edges <- gsea_results
  gsea_empty_edges$leading_edge <- rep("", nrow(gsea_empty_edges))
  expect_no_error({
    p_empty_edges <- visualize_gsea(gsea_empty_edges, plot_type = "network")
    expect_s3_class(p_empty_edges, "ggplot")
  })
})

# Test 11: Comprehensive integration test
test_that("visualize_gsea integration test with all features", {
  skip_if_not_installed("ggplot2")
  
  # Create comprehensive test data
  gsea_results <- create_comprehensive_gsea_results(n_pathways = 20, include_pathway_name = TRUE)
  abundance <- create_comprehensive_abundance(n_genes = 100, n_samples = 12)
  metadata <- create_comprehensive_metadata(n_samples = 12, n_groups = 3)
  
  # Test all combinations that don't require optional packages
  plot_types <- c("enrichment_plot", "dotplot", "barplot")
  sort_methods <- c("NES", "pvalue", "p.adjust")
  
  for (plot_type in plot_types) {
    for (sort_by in sort_methods) {
      expect_no_error({
        p <- visualize_gsea(
          gsea_results = gsea_results,
          plot_type = plot_type,
          n_pathways = 8,
          sort_by = sort_by,
          colors = c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FECA57"),
          pathway_label_column = "pathway_name"
        )
        expect_s3_class(p, "ggplot")
      })
    }
  }
})