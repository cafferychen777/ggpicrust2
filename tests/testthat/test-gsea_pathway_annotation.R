# Test for gsea_pathway_annotation function
# create_test_gsea_results() is defined in helper-gsea.R

test_that("gsea_pathway_annotation annotates KEGG pathways correctly", {
  gsea_results <- create_test_gsea_results()

  kegg_ref <- ggpicrust2:::load_reference_data("KEGG")
  gsea_results$pathway_id <- head(kegg_ref$pathway, 10)

  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")

  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))
  expect_true("pathway_name" %in% colnames(annotated))
  expect_true(any(annotated$pathway_name != annotated$pathway_id))
})

test_that("gsea_pathway_annotation annotates MetaCyc pathways correctly", {
  mock_results <- data.frame(
    pathway_id = c("1CMET2-PWY", "ANAGLYCOLYSIS-PWY", "UNKNOWN-PATHWAY-TEST"),
    NES = c(1.5, -1.2, 0.2),
    pvalue = c(0.02, 0.08, 0.45),
    stringsAsFactors = FALSE
  )

  annotated <- gsea_pathway_annotation(mock_results, pathway_type = "MetaCyc")

  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), 3)
  expect_true("pathway_name" %in% colnames(annotated))
  expect_true(all(!is.na(annotated$pathway_name)))

  # Unknown pathway should keep its ID as name
  unknown_row <- annotated[annotated$pathway_id == "UNKNOWN-PATHWAY-TEST", ]
  expect_equal(unknown_row$pathway_name, "UNKNOWN-PATHWAY-TEST")
})

test_that("gsea_pathway_annotation annotates GO pathways", {
  gsea_results <- create_test_gsea_results()
  gsea_results$pathway_id <- paste0("GO:", sprintf("%07d", 1:10))

  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "GO")

  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))
  expect_true("pathway_name" %in% colnames(annotated))
})

test_that("gsea_pathway_annotation validates inputs correctly", {
  gsea_results <- create_test_gsea_results()

  expect_error(gsea_pathway_annotation(gsea_results = "invalid"), "'gsea_results' must be a data frame")
  expect_error(gsea_pathway_annotation(gsea_results, pathway_type = "invalid"), "pathway_type must be one of")

  gsea_results_missing <- gsea_results[, !names(gsea_results) %in% c("pathway_id")]
  expect_error(gsea_pathway_annotation(gsea_results_missing), "missing required column: pathway_id")
})
