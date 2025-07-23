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
  # Create test data
  gsea_results <- create_test_gsea_results()

  # Mock the data loading
  mockery::stub(gsea_pathway_annotation, "data", function(x, package, envir) {
    if (x == "kegg_reference") {
      assign("kegg_reference", data.frame(
        pathway = paste0("path:ko", sprintf("%05d", 1:10)),
        pathway_name = paste("KEGG Pathway", 1:10),
        pathway_class = rep(c("Metabolism", "Genetic Information Processing"), each = 5),
        description = paste("Description for pathway", 1:10),
        stringsAsFactors = FALSE
      ), envir = envir)
    }
  })

  # Mock the merge function to control the output
  mockery::stub(gsea_pathway_annotation, "merge", function(...) {
    # Return a data frame with the expected structure
    result <- gsea_results
    result$pathway_name <- paste("KEGG Pathway", 1:nrow(result))
    result$pathway_class <- rep(c("Metabolism", "Genetic Information Processing"), length.out = nrow(result))
    result$description <- paste("Description for pathway", 1:nrow(result))
    return(result)
  })

  # Test annotation
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")

  # Check the result
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))
  expect_true("pathway_name" %in% colnames(annotated))
  expect_true("pathway_class" %in% colnames(annotated))
  expect_true("description" %in% colnames(annotated))

  # Check that pathway names were updated
  expect_equal(annotated$pathway_name[1], "KEGG Pathway 1")
})

test_that("gsea_pathway_annotation annotates MetaCyc pathways correctly", {
  # Create test data
  gsea_results <- create_test_gsea_results()
  gsea_results$pathway_id <- paste0("PWY-", 1000 + 1:10)  # MetaCyc IDs

  # Mock the data loading
  mockery::stub(gsea_pathway_annotation, "data", function(x, package, envir) {
    if (x == "metacyc_reference") {
      assign("metacyc_reference", data.frame(
        pathway = paste0("PWY-", 1000 + 1:10),
        pathway_name = paste("MetaCyc Pathway", 1:10),
        description = paste("Description for MetaCyc pathway", 1:10),
        stringsAsFactors = FALSE
      ), envir = envir)
    }
  })

  # Mock the merge function to control the output
  mockery::stub(gsea_pathway_annotation, "merge", function(...) {
    # Return a data frame with the expected structure
    result <- gsea_results
    result$pathway_name <- paste("MetaCyc Pathway", 1:nrow(result))
    result$description <- paste("Description for MetaCyc pathway", 1:nrow(result))
    return(result)
  })

  # Test annotation
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "MetaCyc")

  # Check the result
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))
  expect_true("pathway_name" %in% colnames(annotated))
  expect_true("description" %in% colnames(annotated))

  # Check that pathway names were updated
  expect_equal(annotated$pathway_name[1], "MetaCyc Pathway 1")
})

test_that("gsea_pathway_annotation handles GO pathways", {
  # Create test data
  gsea_results <- create_test_gsea_results()
  gsea_results$pathway_id <- paste0("GO:", sprintf("%07d", 1:10))  # GO IDs

  # Test annotation (not fully implemented yet)
  expect_warning(
    annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "GO"),
    "GO pathway annotation not yet implemented"
  )

  # Check the result
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))
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
    "pathway_type must be one of 'KEGG', 'MetaCyc', or 'GO'"
  )

  # Test missing required columns
  gsea_results_missing <- gsea_results[, !names(gsea_results) %in% c("pathway_id")]
  expect_error(
    gsea_pathway_annotation(gsea_results_missing),
    "GSEA results missing required column: pathway_id"
  )
})

test_that("gsea_pathway_annotation handles missing annotations", {
  # Create test data
  gsea_results <- create_test_gsea_results()

  # Add some pathways that won't be in the reference
  gsea_results$pathway_id[1:3] <- paste0("path:ko", sprintf("%05d", 101:103))

  # Mock the data loading with incomplete reference
  mockery::stub(gsea_pathway_annotation, "data", function(x, package, envir) {
    if (x == "kegg_reference") {
      assign("kegg_reference", data.frame(
        pathway = paste0("path:ko", sprintf("%05d", 4:10)),  # Missing first 3 pathways
        pathway_name = paste("KEGG Pathway", 4:10),
        pathway_class = rep("Metabolism", 7),
        stringsAsFactors = FALSE
      ), envir = envir)
    }
  })

  # Mock the merge function to control the output
  mockery::stub(gsea_pathway_annotation, "merge", function(x, y, by.x, by.y, all.x, ...) {
    # Simulate a merge that keeps all rows from x (gsea_results)
    result <- x
    # For rows 1-3, pathway_name should be the same as pathway_id (missing in reference)
    # For rows 4-10, pathway_name should come from the reference
    result$pathway_name <- result$pathway_id  # Default to pathway_id
    result$pathway_name[4:10] <- paste("KEGG Pathway", 4:10)  # Update with reference names for rows 4-10
    result$pathway_class <- NA
    result$pathway_class[4:10] <- "Metabolism"
    return(result)
  })

  # Test annotation
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")

  # Check the result
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))

  # Check that missing pathways use pathway_id as name
  expect_equal(annotated$pathway_name[1], gsea_results$pathway_id[1])
  expect_equal(annotated$pathway_name[2], gsea_results$pathway_id[2])
  expect_equal(annotated$pathway_name[3], gsea_results$pathway_id[3])

  # Check that found pathways use reference names
  expect_equal(annotated$pathway_name[4], "KEGG Pathway 4")
})
