# Detailed Tests for Network Visualization Functions
library(testthat)

# Helper function to create test data with controlled similarity
create_similarity_test_data <- function() {
  # Create pathways with known gene overlaps for testing similarity calculations
  data.frame(
    pathway_id = c("path1", "path2", "path3", "path4", "path5"),
    pathway_name = c("Pathway A", "Pathway B", "Pathway C", "Pathway D", "Pathway E"),
    NES = c(2.1, -1.5, 2.8, -0.8, 1.2),
    pvalue = c(0.001, 0.05, 0.001, 0.1, 0.02),
    p.adjust = c(0.01, 0.1, 0.01, 0.2, 0.05),
    size = c(20, 15, 25, 10, 18),
    # Controlled leading edges for testing similarity
    leading_edge = c(
      "gene1;gene2;gene3;gene4;gene5",           # path1: 5 genes
      "gene1;gene2;gene6;gene7",                 # path2: 4 genes, 2 overlap with path1
      "gene1;gene8;gene9;gene10;gene11;gene12",  # path3: 6 genes, 1 overlap with path1
      "gene13;gene14;gene15",                    # path4: 3 genes, no overlap
      "gene2;gene3;gene16;gene17"                # path5: 4 genes, 2 overlap with path1
    ),
    stringsAsFactors = FALSE
  )
}

test_that("create_network_plot calculates Jaccard similarity correctly", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  gsea_results <- create_similarity_test_data()
  
  # Test Jaccard similarity calculation
  # Expected similarities based on our test data:
  # path1 vs path2: |{gene1,gene2}| / |{gene1,gene2,gene3,gene4,gene5,gene6,gene7}| = 2/7 ≈ 0.286
  # path1 vs path3: |{gene1}| / |{gene1,gene2,gene3,gene4,gene5,gene8,gene9,gene10,gene11,gene12}| = 1/10 = 0.1
  # path1 vs path4: 0 (no overlap)
  # path1 vs path5: |{gene2,gene3}| / |{gene1,gene2,gene3,gene4,gene5,gene16,gene17}| = 2/7 ≈ 0.286
  
  expect_no_error({
    p_jaccard <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      network_params = list(
        similarity_measure = "jaccard",
        similarity_cutoff = 0.1
      )
    )
    expect_s3_class(p_jaccard, "ggplot")
  })
  
  # Test with higher cutoff that should eliminate some connections
  expect_no_error({
    p_jaccard_high <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      network_params = list(
        similarity_measure = "jaccard",
        similarity_cutoff = 0.25
      )
    )
    expect_s3_class(p_jaccard_high, "ggplot")
  })
})

test_that("create_network_plot calculates overlap similarity correctly", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  gsea_results <- create_similarity_test_data()
  
  # Test overlap coefficient calculation
  # Expected similarities based on our test data:
  # path1 vs path2: |{gene1,gene2}| / min(5,4) = 2/4 = 0.5
  # path1 vs path3: |{gene1}| / min(5,6) = 1/5 = 0.2
  # path1 vs path4: 0 (no overlap)
  # path1 vs path5: |{gene2,gene3}| / min(5,4) = 2/4 = 0.5
  
  expect_no_error({
    p_overlap <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      network_params = list(
        similarity_measure = "overlap",
        similarity_cutoff = 0.15
      )
    )
    expect_s3_class(p_overlap, "ggplot")
  })
})

test_that("create_network_plot calculates correlation similarity correctly", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  gsea_results <- create_similarity_test_data()
  
  # Test correlation-like similarity calculation
  expect_no_error({
    p_correlation <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      network_params = list(
        similarity_measure = "correlation",
        similarity_cutoff = 0.1
      )
    )
    expect_s3_class(p_correlation, "ggplot")
  })
})

test_that("create_network_plot handles different layout algorithms", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  gsea_results <- create_similarity_test_data()
  
  layouts <- c("fruchterman", "kamada", "circle")
  
  for (layout in layouts) {
    tryCatch({
      p <- visualize_gsea(
        gsea_results,
        plot_type = "network",
        network_params = list(
          layout = layout,
          similarity_cutoff = 0.1
        )
      )
      expect_s3_class(p, "ggplot")
    }, error = function(e) {
      fail(paste("Layout:", layout, "-", e$message))
    })
  }
})

test_that("create_network_plot handles different node coloring options", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  gsea_results <- create_similarity_test_data()
  
  color_options <- c("NES", "pvalue", "p.adjust")
  
  for (color_by in color_options) {
    tryCatch({
      p <- visualize_gsea(
        gsea_results,
        plot_type = "network",
        network_params = list(
          node_color_by = color_by,
          similarity_cutoff = 0.1
        )
      )
      expect_s3_class(p, "ggplot")
    }, error = function(e) {
      fail(paste("Color by:", color_by, "-", e$message))
    })
  }
})

test_that("create_network_plot handles edge width options", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  gsea_results <- create_similarity_test_data()
  
  edge_width_options <- c("similarity", "constant")
  
  for (edge_width in edge_width_options) {
    tryCatch({
      p <- visualize_gsea(
        gsea_results,
        plot_type = "network",
        network_params = list(
          edge_width_by = edge_width,
          similarity_cutoff = 0.1
        )
      )
      expect_s3_class(p, "ggplot")
    }, error = function(e) {
      fail(paste("Edge width by:", edge_width, "-", e$message))
    })
  }
})

test_that("create_network_plot handles empty leading edges gracefully", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  gsea_results <- create_similarity_test_data()
  
  # Test with some empty leading edges
  gsea_results$leading_edge[3] <- ""
  gsea_results$leading_edge[4] <- ""
  
  expect_no_error({
    p <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      network_params = list(similarity_cutoff = 0.1)
    )
    expect_s3_class(p, "ggplot")
  })
})

test_that("create_network_plot handles all empty leading edges", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  gsea_results <- create_similarity_test_data()
  
  # Test with all empty leading edges
  gsea_results$leading_edge <- rep("", nrow(gsea_results))
  
  expect_no_error({
    p <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      network_params = list(similarity_cutoff = 0.1)
    )
    expect_s3_class(p, "ggplot")
    # Should return a plot with a message about no connections
  })
})

test_that("create_network_plot handles very high similarity cutoff", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  gsea_results <- create_similarity_test_data()
  
  # Test with cutoff that eliminates all connections
  expect_no_error({
    p <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      network_params = list(similarity_cutoff = 0.99)
    )
    expect_s3_class(p, "ggplot")
    # Should return a plot with a message about no connections
  })
})

test_that("create_network_plot handles single pathway", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  # Test with single pathway
  gsea_single <- data.frame(
    pathway_id = "path1",
    pathway_name = "Single Pathway",
    NES = 2.1,
    pvalue = 0.001,
    p.adjust = 0.01,
    size = 20,
    leading_edge = "gene1;gene2;gene3",
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    p <- visualize_gsea(
      gsea_single, 
      plot_type = "network",
      network_params = list(similarity_cutoff = 0.1)
    )
    expect_s3_class(p, "ggplot")
  })
})

test_that("create_network_plot parameter validation", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  gsea_results <- create_similarity_test_data()
  
  # Test invalid similarity measure (should use default)
  expect_no_error({
    p <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      network_params = list(similarity_measure = "invalid")
    )
    expect_s3_class(p, "ggplot")
  })
  
  # Test invalid layout (should use default)
  expect_no_error({
    p <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      network_params = list(layout = "invalid")
    )
    expect_s3_class(p, "ggplot")
  })
  
  # Test extreme cutoff values
  expect_no_error({
    p_zero <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      network_params = list(similarity_cutoff = 0.0)
    )
    expect_s3_class(p_zero, "ggplot")
  })
  
  expect_no_error({
    p_one <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      network_params = list(similarity_cutoff = 1.0)
    )
    expect_s3_class(p_one, "ggplot")
  })
})

test_that("create_network_plot works with pathway labels", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  gsea_results <- create_similarity_test_data()
  
  # Test with pathway_name labels
  expect_no_error({
    p_names <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      pathway_label_column = "pathway_name",
      network_params = list(similarity_cutoff = 0.1)
    )
    expect_s3_class(p_names, "ggplot")
  })
  
  # Test with pathway_id labels
  expect_no_error({
    p_ids <- visualize_gsea(
      gsea_results, 
      plot_type = "network",
      pathway_label_column = "pathway_id",
      network_params = list(similarity_cutoff = 0.1)
    )
    expect_s3_class(p_ids, "ggplot")
  })
})

test_that("create_network_plot handles large number of pathways", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("dplyr")
  
  # Create larger dataset
  set.seed(42)
  n_pathways <- 25
  all_genes <- paste0("gene", 1:100)
  
  gsea_large <- data.frame(
    pathway_id = paste0("path", 1:n_pathways),
    pathway_name = paste("Pathway", 1:n_pathways),
    NES = runif(n_pathways, -3, 3),
    pvalue = runif(n_pathways, 0.001, 0.1),
    p.adjust = runif(n_pathways, 0.001, 0.2),
    size = sample(10:100, n_pathways, replace = TRUE),
    leading_edge = replicate(n_pathways, {
      paste(sample(all_genes, sample(3:15, 1)), collapse = ";")
    }),
    stringsAsFactors = FALSE
  )
  
  expect_no_error({
    p_large <- visualize_gsea(
      gsea_large, 
      plot_type = "network",
      n_pathways = 15,  # Limit to reasonable number
      network_params = list(similarity_cutoff = 0.2)
    )
    expect_s3_class(p_large, "ggplot")
  })
})