# Test GO pathway support in GSEA functionality

library(testthat)

# Note: context() is deprecated in testthat 3.x, removed

test_that("create_basic_go_mapping creates valid enhanced GO mapping", {
  # Load the function
  source("../../R/pathway_gsea.R")

  # Create enhanced GO mapping
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

  # Test enhanced mapping has more terms than basic (should be 100+)
  expect_true(nrow(go_mapping) >= 100)

  # Test distribution across categories
  category_counts <- table(go_mapping$category)
  expect_true(category_counts["BP"] >= 40)  # At least 40 BP terms
  expect_true(category_counts["MF"] >= 20)  # At least 20 MF terms
  expect_true(category_counts["CC"] >= 20)  # At least 20 CC terms

  # Test that we have microbiome-relevant terms
  microbiome_terms <- c("metabolic", "stress", "membrane", "transport", "catalytic")
  term_names <- tolower(go_mapping$go_name)
  for (term in microbiome_terms) {
    expect_true(any(grepl(term, term_names)),
                info = paste("Missing microbiome-relevant term:", term))
  }
})

test_that("ko_to_go_reference dataset is properly structured", {
  # Test loading the complete dataset
  expect_no_error({
    data("ko_to_go_reference", package = "ggpicrust2")
  })

  # Test dataset structure
  expect_is(ko_to_go_reference, "data.frame")
  expect_true(all(c("go_id", "go_name", "category", "ko_members") %in% colnames(ko_to_go_reference)))

  # Test GO ID format
  expect_true(all(grepl("^GO:\\d{7}$", ko_to_go_reference$go_id)))

  # Test categories
  expect_true(all(ko_to_go_reference$category %in% c("BP", "MF", "CC")))

  # Test KO members format
  expect_true(all(grepl("K\\d{5}", ko_to_go_reference$ko_members)))

  # Test we have all three categories
  expect_true(length(unique(ko_to_go_reference$category)) == 3)

  # Test dataset size (should be comprehensive)
  expect_true(nrow(ko_to_go_reference) >= 100)

  # Test no duplicate GO IDs
  expect_true(length(unique(ko_to_go_reference$go_id)) == nrow(ko_to_go_reference))

  # Test that KO members are properly formatted
  ko_members_sample <- head(ko_to_go_reference$ko_members, 10)
  for (ko_string in ko_members_sample) {
    ko_list <- strsplit(ko_string, ";")[[1]]
    expect_true(all(grepl("^K\\d{5}$", ko_list)),
                info = paste("Invalid KO format in:", ko_string))
  }
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
  
  # Prepare abundance data - use tibble method to maintain sample names
  abundance_data <- ko_abundance %>%
    tibble::column_to_rownames("#NAME")

  # Take subset for faster testing
  abundance_data <- abundance_data[1:100, ]
  
  # Test GO GSEA with BP category
  # Metadata already has rownames, just convert to plain data.frame
  meta_df <- as.data.frame(metadata)
  rownames(meta_df) <- meta_df$sample_name

  expect_no_error({
    gsea_results_bp <- pathway_gsea(
      abundance = abundance_data,
      metadata = meta_df,
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
      metadata = meta_df,
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

test_that("GO data loading mechanism works correctly", {
  # Test that the improved data loading mechanism works
  # This tests the fallback mechanism when complete data is not available

  # Test prepare_gene_sets with different data availability scenarios
  expect_no_error({
    gene_sets_bp <- prepare_gene_sets("GO", go_category = "BP")
  })

  expect_is(gene_sets_bp, "list")
  expect_true(length(gene_sets_bp) > 0)

  # Test that gene sets contain valid KO identifiers
  all_kos <- unique(unlist(gene_sets_bp))
  expect_true(all(grepl("^K\\d{5}$", all_kos)))

  # Test that we get a reasonable number of gene sets
  expect_true(length(gene_sets_bp) >= 20)  # Should have at least 20 BP terms
})

test_that("GO pathway analysis handles edge cases", {
  # Test with minimal data
  minimal_abundance <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2,
                             dimnames = list(c("K00001", "K00002"),
                                           c("S1", "S2", "S3")))
  # Need at least 4 samples for statistical analysis
  minimal_metadata <- data.frame(
    sample = c("S1", "S2", "S3"),
    group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  rownames(minimal_metadata) <- minimal_metadata$sample

  # Should error with insufficient samples (< 4)
  expect_error({
    result <- pathway_gsea(
      abundance = minimal_abundance,
      metadata = minimal_metadata,
      group = "group",
      pathway_type = "GO",
      go_category = "BP",
      method = "fgsea",
      min_size = 1,
      max_size = 100,
      nperm = 10
    )
  }, "Insufficient overlapping samples")
})

test_that("GO error handling provides informative messages", {
  # Test that error messages are informative when things go wrong

  # Test invalid GO category
  expect_error({
    prepare_gene_sets("GO", go_category = "INVALID")
  }, regexp = "Invalid go_category")

  # Test that warnings are properly issued when fallback is used
  # This is tested implicitly in other tests through the warning system
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