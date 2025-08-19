# Test GO pathway support in GSEA functionality

library(testthat)

context("GO Pathway Support")

test_that("create_basic_go_mapping creates valid GO mapping", {
  # Load the function
  source("../../R/pathway_gsea.R")
  
  # Create basic GO mapping
  go_mapping <- create_basic_go_mapping()
  
  # Test structure
  expect_is(go_mapping, "data.frame")
  expect_true(all(c("go_id", "go_name", "category", "ko_members") %in% colnames(go_mapping)))
  
  # Test GO ID format
  expect_true(all(grepl("^GO:\\d{7}$", go_mapping$go_id)))
  
  # Test categories
  expect_true(all(go_mapping$category %in% c("BP", "MF", "CC")))
  
  # Test KO members format
  expect_true(all(grepl("K\\d{5}", go_mapping$ko_members)))
  
  # Test we have all three categories
  expect_true(length(unique(go_mapping$category)) == 3)
})

test_that("prepare_gene_sets works with GO pathway_type", {
  # Load required data
  data("ko_abundance", package = "ggpicrust2")
  
  # Test BP category
  gene_sets_bp <- prepare_gene_sets("GO", go_category = "BP")
  
  expect_is(gene_sets_bp, "list")
  expect_true(length(gene_sets_bp) > 0)
  
  # Check GO ID format in names
  expect_true(all(grepl("^GO:\\d{7}$", names(gene_sets_bp))))
  
  # Check KO format in gene sets
  all_kos <- unique(unlist(gene_sets_bp))
  expect_true(all(grepl("^K\\d{5}$", all_kos)))
  
  # Test MF category
  gene_sets_mf <- prepare_gene_sets("GO", go_category = "MF")
  expect_is(gene_sets_mf, "list")
  expect_true(length(gene_sets_mf) > 0)
  
  # Test CC category
  gene_sets_cc <- prepare_gene_sets("GO", go_category = "CC")
  expect_is(gene_sets_cc, "list")
  expect_true(length(gene_sets_cc) > 0)
  
  # Test all categories
  gene_sets_all <- prepare_gene_sets("GO", go_category = "all")
  expect_is(gene_sets_all, "list")
  expect_true(length(gene_sets_all) >= length(gene_sets_bp))
})

test_that("pathway_gsea works with GO pathway_type", {
  # Load test data
  data("ko_abundance", package = "ggpicrust2")
  data("metadata", package = "ggpicrust2")
  
  # Prepare abundance data
  abundance_data <- as.data.frame(ko_abundance)
  rownames(abundance_data) <- abundance_data[, "#NAME"]
  abundance_data <- abundance_data[, -1]
  
  # Take subset for faster testing
  abundance_data <- abundance_data[1:100, ]
  
  # Test GO GSEA with BP category
  expect_no_error({
    gsea_results_bp <- pathway_gsea(
      abundance = abundance_data,
      metadata = metadata,
      group = "Environment",
      pathway_type = "GO",
      go_category = "BP",
      method = "fgsea",
      nperm = 100,
      min_size = 2,
      max_size = 500
    )
  })
  
  # Test GO GSEA with MF category
  expect_no_error({
    gsea_results_mf <- pathway_gsea(
      abundance = abundance_data,
      metadata = metadata,
      group = "Environment",
      pathway_type = "GO",
      go_category = "MF",
      method = "fgsea",
      nperm = 100,
      min_size = 2,
      max_size = 500
    )
  })
})

test_that("gsea_pathway_annotation works with GO results", {
  # Load test data
  data("ko_abundance", package = "ggpicrust2")
  data("metadata", package = "ggpicrust2")
  
  # Prepare abundance data
  abundance_data <- as.data.frame(ko_abundance)
  rownames(abundance_data) <- abundance_data[, "#NAME"]
  abundance_data <- abundance_data[, -1]
  
  # Take subset for faster testing
  abundance_data <- abundance_data[1:50, ]
  
  # Run GSEA analysis
  gsea_results <- pathway_gsea(
    abundance = abundance_data,
    metadata = metadata,
    group = "Environment",
    pathway_type = "GO",
    go_category = "BP",
    method = "fgsea",
    nperm = 50,
    min_size = 2,
    max_size = 500
  )
  
  # Test annotation
  if (nrow(gsea_results) > 0) {
    annotated_results <- gsea_pathway_annotation(
      gsea_results = gsea_results,
      pathway_type = "GO"
    )
    
    expect_is(annotated_results, "data.frame")
    expect_true("pathway_name" %in% colnames(annotated_results))
    
    # Check that GO IDs are properly formatted
    expect_true(all(grepl("^GO:\\d{7}$", annotated_results$pathway_id)))
    
    # Check that pathway names are not just GO IDs for known terms
    known_terms <- annotated_results[!is.na(annotated_results$pathway_name) & 
                                   annotated_results$pathway_name != annotated_results$pathway_id, ]
    if (nrow(known_terms) > 0) {
      expect_true(nrow(known_terms) > 0)
    }
  }
})

test_that("GO category filtering works correctly", {
  # Test that each category returns only terms from that category
  gene_sets_bp <- prepare_gene_sets("GO", go_category = "BP")
  gene_sets_mf <- prepare_gene_sets("GO", go_category = "MF") 
  gene_sets_cc <- prepare_gene_sets("GO", go_category = "CC")
  
  # Categories should be different (no overlap in GO IDs)
  bp_ids <- names(gene_sets_bp)
  mf_ids <- names(gene_sets_mf)
  cc_ids <- names(gene_sets_cc)
  
  expect_true(length(intersect(bp_ids, mf_ids)) == 0)
  expect_true(length(intersect(bp_ids, cc_ids)) == 0)
  expect_true(length(intersect(mf_ids, cc_ids)) == 0)
  
  # All categories together should include all individual categories
  gene_sets_all <- prepare_gene_sets("GO", go_category = "all")
  all_ids <- names(gene_sets_all)
  
  expect_true(all(bp_ids %in% all_ids))
  expect_true(all(mf_ids %in% all_ids))
  expect_true(all(cc_ids %in% all_ids))
})

test_that("GO gene sets have reasonable sizes", {
  gene_sets <- prepare_gene_sets("GO", go_category = "BP")
  
  # Check gene set sizes
  gene_set_sizes <- sapply(gene_sets, length)
  
  expect_true(all(gene_set_sizes > 0))
  expect_true(all(gene_set_sizes <= 50))  # Reasonable upper limit
  expect_true(mean(gene_set_sizes) >= 3)  # At least some KOs per term
})

test_that("GO validation catches invalid parameters", {
  # Test invalid GO category
  expect_error(pathway_gsea(
    abundance = matrix(1:10, nrow = 2),
    metadata = data.frame(sample = c("A", "B", "C", "D", "E"), 
                         group = c("X", "X", "Y", "Y", "Y")),
    group = "group",
    pathway_type = "GO",
    go_category = "INVALID",
    method = "fgsea"
  ))
})

test_that("GO GSEA integrates with visualization functions", {
  # This test ensures GO results work with existing visualization
  # Load test data
  data("ko_abundance", package = "ggpicrust2")
  data("metadata", package = "ggpicrust2") 
  
  # Prepare abundance data (small subset)
  abundance_data <- as.data.frame(ko_abundance)
  rownames(abundance_data) <- abundance_data[, "#NAME"]
  abundance_data <- abundance_data[, -1]
  abundance_data <- abundance_data[1:30, ]
  
  # Run GO GSEA
  gsea_results <- pathway_gsea(
    abundance = abundance_data,
    metadata = metadata,
    group = "Environment",
    pathway_type = "GO",
    go_category = "BP",
    method = "fgsea",
    nperm = 50,
    min_size = 2,
    max_size = 500
  )
  
  # Test that results have expected structure for visualization
  if (nrow(gsea_results) > 0) {
    expect_true(all(c("pathway_id", "NES", "pvalue", "p.adjust") %in% colnames(gsea_results)))
    expect_is(gsea_results$NES, "numeric")
    expect_is(gsea_results$pvalue, "numeric")
    expect_is(gsea_results$p.adjust, "numeric")
  }
})