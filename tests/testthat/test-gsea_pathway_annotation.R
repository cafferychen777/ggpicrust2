# Test for gsea_pathway_annotation function
library(testthat)

# Helper function to create test GSEA results
create_test_gsea_results <- function(n_pathways = 10) {
  set.seed(123)
  data.frame(
    pathway_id = paste0("path:ko", sprintf("%05d", 1:n_pathways)),
    pathway_name = paste0("path:ko", sprintf("%05d", 1:n_pathways)),  # Same as pathway_id, will be replaced
    size = sample(10:100, n_pathways, replace = TRUE),
    ES = runif(n_pathways, -0.8, 0.8),
    NES = runif(n_pathways, -2, 2),
    pvalue = runif(n_pathways, 0, 0.1),
    p.adjust = runif(n_pathways, 0, 0.2),
    leading_edge = replicate(n_pathways, paste(paste0("K", sprintf("%05d", sample(1:1000, 5))), collapse = ";")),
    method = rep("fgsea", n_pathways),
    stringsAsFactors = FALSE
  )
}

test_that("gsea_pathway_annotation annotates KEGG pathways correctly", {
  # Create test data using real KEGG pathway IDs from the reference
  gsea_results <- create_test_gsea_results()

  # Use real pathway IDs that exist in the reference data
  kegg_ref <- ggpicrust2:::load_reference_data("KEGG")
  real_pathway_ids <- head(kegg_ref$pathway, 10)
  gsea_results$pathway_id <- real_pathway_ids

  # Test annotation
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")

  # Check the result
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))
  expect_true("pathway_name" %in% colnames(annotated))

  # Check that pathway names were updated (should not just be the pathway_id)
  # At least some pathways should have different names from their IDs
  expect_true(any(annotated$pathway_name != annotated$pathway_id))
})

test_that("gsea_pathway_annotation annotates MetaCyc pathways correctly", {
  # Create test data with some common MetaCyc pathway IDs
  gsea_results <- create_test_gsea_results()
  gsea_results$pathway_id <- paste0("PWY-", 1000 + 1:10)  # MetaCyc IDs

  # Test annotation - the function should work with real reference data
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "MetaCyc")

  # Check the result structure
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))
  expect_true("pathway_name" %in% colnames(annotated))

  # For unknown pathways, pathway_name should default to pathway_id
  # (since these are not real MetaCyc IDs in the reference)
  expect_true(all(!is.na(annotated$pathway_name)))
})

test_that("gsea_pathway_annotation handles GO pathways", {
  # Create test data
  gsea_results <- create_test_gsea_results()
  gsea_results$pathway_id <- paste0("GO:", sprintf("%07d", 1:10))  # GO IDs

  # GO annotation is now fully implemented using ko_to_go_reference
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "GO")

  # Check the result
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))
  expect_true("pathway_name" %in% colnames(annotated))
})

test_that("gsea_pathway_annotation validates inputs correctly", {
  # Create test data
  gsea_results <- create_test_gsea_results()

  # Test invalid gsea_results
  expect_error(
    gsea_pathway_annotation(gsea_results = "invalid"),
    "'gsea_results' must be a data frame"
  )

  # Test invalid pathway_type
  expect_error(
    gsea_pathway_annotation(gsea_results, pathway_type = "invalid"),
    "pathway_type must be one of: KEGG, MetaCyc, GO"
  )

  # Test missing required columns
  gsea_results_missing <- gsea_results[, !names(gsea_results) %in% c("pathway_id")]
  expect_error(
    gsea_pathway_annotation(gsea_results_missing),
    "GSEA results missing required column: pathway_id"
  )
})

test_that("gsea_pathway_annotation handles missing annotations", {
  # Create test data with a mix of real and fake pathway IDs
  gsea_results <- create_test_gsea_results()

  # Use some real pathway IDs and some fake ones
  kegg_ref <- ggpicrust2:::load_reference_data("KEGG")
  real_pathway_ids <- head(kegg_ref$pathway, 7)

  # Create fake pathway IDs that don't exist in reference
  fake_pathway_ids <- c("path:ko99991", "path:ko99992", "path:ko99993")

  gsea_results$pathway_id <- c(fake_pathway_ids, real_pathway_ids)

  # Test annotation
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")

  # Check the result
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))

  # The missing pathways (fake IDs) should use pathway_id as name
  missing_rows <- annotated[annotated$pathway_id %in% fake_pathway_ids, ]
  expect_true(all(missing_rows$pathway_name == missing_rows$pathway_id))

  # Found pathways should have reference names (not equal to pathway_id)
  found_rows <- annotated[annotated$pathway_id %in% real_pathway_ids, ]
  expect_true(any(found_rows$pathway_name != found_rows$pathway_id))
})
