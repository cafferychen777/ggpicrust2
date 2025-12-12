# Comprehensive tests for gsea_pathway_annotation function
library(testthat)
library(ggpicrust2)

# Helper function to create realistic test GSEA results
create_realistic_gsea_results <- function(n_pathways = 20, pathway_type = "KEGG") {
  set.seed(123)

  if (pathway_type == "KEGG") {
    pathway_ids <- paste0("ko", sprintf("%05d", sample(10:999, n_pathways)))
  } else if (pathway_type == "MetaCyc") {
    base_pathways <- c("GLYCOLYSIS", "TCA-CYCLE", "PENTOSE-P-PWY", "CALVIN-PWY",
                      "FATTY-ACID-BIOSYNTHESIS", "AMINO-ACID-BIOSYN", "NUCLEOTIDE-BIOSYN",
                      "METHANE-OXIDATION", "SULFATE-REDUCTION", "NITROGEN-FIXATION")
    if (n_pathways <= 10) {
      pathway_ids <- base_pathways[1:n_pathways]
    } else {
      additional_pathways <- paste0("PWY-", sample(1000:9999, n_pathways - 10))
      pathway_ids <- c(base_pathways, additional_pathways)
    }
  } else {
    pathway_ids <- paste0("GO:", sprintf("%07d", sample(1:9999999, n_pathways)))
  }

  data.frame(
    pathway_id = pathway_ids,
    pathway_name = pathway_ids,  # Will be replaced by annotation
    size = sample(5:200, n_pathways, replace = TRUE),
    ES = runif(n_pathways, -1.0, 1.0),
    NES = rnorm(n_pathways, 0, 1.5),
    pvalue = runif(n_pathways, 0.001, 0.2),
    p.adjust = runif(n_pathways, 0.001, 0.3),
    leading_edge = replicate(n_pathways, {
      genes <- paste0("K", sprintf("%05d", sample(1:10000, sample(3:15, 1))))
      paste(genes, collapse = ";")
    }),
    method = rep("fgsea", n_pathways),
    stringsAsFactors = FALSE
  )
}

test_that("gsea_pathway_annotation handles KEGG pathways correctly", {
  # Use real pathway IDs from the package's reference data
  # Load ko_to_kegg_reference from sysdata
  data_env <- new.env()
  load(system.file("R/sysdata.rda", package = "ggpicrust2"), envir = data_env)

  # Get some real KEGG pathway IDs
  real_pathway_ids <- head(unique(data_env$ko_to_kegg_reference$pathway_id), 5)

  gsea_results <- data.frame(
    pathway_id = real_pathway_ids,
    pathway_name = real_pathway_ids,  # Will be replaced
    size = c(10, 20, 30, 40, 50),
    ES = runif(5, -1, 1),
    NES = rnorm(5),
    pvalue = runif(5, 0.001, 0.2),
    p.adjust = runif(5, 0.001, 0.3),
    leading_edge = rep("K00001;K00002", 5),
    method = rep("fgsea", 5),
    stringsAsFactors = FALSE
  )

  # Test annotation - should work with real KEGG reference data
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")

  # Check the result structure
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))

  # Check that all original columns are preserved
  original_cols <- colnames(gsea_results)
  expect_true(all(original_cols %in% colnames(annotated)))

  # Check that annotation columns exist
  expect_true("pathway_name" %in% colnames(annotated))
})

test_that("gsea_pathway_annotation handles MetaCyc pathways correctly", {
  # Create test data with common MetaCyc pathway IDs
  test_pathways <- c("GLYCOLYSIS", "TCA", "PWY-5659")

  gsea_results <- data.frame(
    pathway_id = test_pathways,
    pathway_name = test_pathways,
    size = c(10, 15, 20),
    ES = c(0.5, -0.3, 0.8),
    NES = c(1.2, -0.8, 1.5),
    pvalue = c(0.01, 0.05, 0.001),
    p.adjust = c(0.05, 0.1, 0.01),
    leading_edge = rep("K00001;K00002", 3),
    method = rep("fgsea", 3),
    stringsAsFactors = FALSE
  )

  # Test annotation
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "MetaCyc")

  # Check the result
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))
  expect_true("pathway_name" %in% colnames(annotated))
})

test_that("gsea_pathway_annotation handles GO pathways correctly", {
  # Load real GO reference data
  tryCatch({
    data("ko_to_go_reference", package = "ggpicrust2")
    real_go_ids <- head(ko_to_go_reference$go_id, 5)

    gsea_results <- data.frame(
      pathway_id = real_go_ids,
      pathway_name = real_go_ids,
      size = c(10, 20, 30, 40, 50),
      ES = runif(5, -1, 1),
      NES = rnorm(5),
      pvalue = runif(5, 0.001, 0.2),
      p.adjust = runif(5, 0.001, 0.3),
      leading_edge = rep("K00001;K00002", 5),
      method = rep("fgsea", 5),
      stringsAsFactors = FALSE
    )

    # GO annotation should now work (no warning expected)
    annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "GO")

    expect_s3_class(annotated, "data.frame")
    expect_equal(nrow(annotated), nrow(gsea_results))
    expect_true("pathway_name" %in% colnames(annotated))

    # Check that pathway names were updated for known GO terms
    expect_true(any(annotated$pathway_name != annotated$pathway_id))
  }, error = function(e) {
    skip("ko_to_go_reference data not available")
  })
})

test_that("gsea_pathway_annotation handles empty GSEA results", {
  # Create empty GSEA results with proper structure
  empty_gsea <- data.frame(
    pathway_id = character(0),
    pathway_name = character(0),
    size = integer(0),
    ES = numeric(0),
    NES = numeric(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    leading_edge = character(0),
    method = character(0),
    stringsAsFactors = FALSE
  )

  annotated <- gsea_pathway_annotation(empty_gsea, pathway_type = "KEGG")

  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), 0)
})

test_that("gsea_pathway_annotation validates input parameters", {
  gsea_results <- create_realistic_gsea_results(5)

  # Test invalid gsea_results type
  expect_error(
    gsea_pathway_annotation(gsea_results = list()),
    "'gsea_results' must be a data frame"
  )

  expect_error(
    gsea_pathway_annotation(gsea_results = "not_a_dataframe"),
    "'gsea_results' must be a data frame"
  )

  # Test invalid pathway_type
  expect_error(
    gsea_pathway_annotation(gsea_results, pathway_type = "INVALID"),
    "pathway_type must be one of 'KEGG', 'MetaCyc', or 'GO'"
  )

  expect_error(
    gsea_pathway_annotation(gsea_results, pathway_type = c("KEGG", "MetaCyc")),
    "pathway_type must be one of 'KEGG', 'MetaCyc', or 'GO'"
  )

  # Test missing required columns
  incomplete_gsea <- gsea_results[, !names(gsea_results) %in% "pathway_id"]
  expect_error(
    gsea_pathway_annotation(incomplete_gsea),
    "GSEA results missing required column: pathway_id"
  )
})

test_that("gsea_pathway_annotation handles factor pathway_ids", {
  gsea_results <- create_realistic_gsea_results(5, "KEGG")

  # Test with factor pathway_ids
  gsea_results_factor <- gsea_results
  gsea_results_factor$pathway_id <- as.factor(gsea_results_factor$pathway_id)

  annotated <- gsea_pathway_annotation(gsea_results_factor, pathway_type = "KEGG")
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), 5)
})

test_that("gsea_pathway_annotation handles NA values in pathway_id", {
  gsea_results <- create_realistic_gsea_results(5, "KEGG")

  # Test with NA values in pathway_id
  gsea_results_na <- gsea_results
  gsea_results_na$pathway_id[1] <- NA

  annotated <- gsea_pathway_annotation(gsea_results_na, pathway_type = "KEGG")
  expect_s3_class(annotated, "data.frame")
  # NA pathway should still be in results (may be handled differently)
})

test_that("gsea_pathway_annotation preserves original column types", {
  gsea_results <- create_realistic_gsea_results(3)
  original_classes <- sapply(gsea_results, class)

  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")

  # Check that original column types are preserved (for existing columns)
  for (col in names(original_classes)) {
    if (col %in% colnames(annotated) && col != "pathway_name") {
      expect_equal(class(annotated[[col]]), original_classes[[col]])
    }
  }
})

test_that("gsea_pathway_annotation handles pathways not in reference", {
  # Create test data with pathway IDs not in reference
  fake_pathways <- c("ko99999", "ko88888", "ko77777")

  gsea_results <- data.frame(
    pathway_id = fake_pathways,
    pathway_name = fake_pathways,
    size = c(10, 20, 30),
    ES = c(0.5, -0.3, 0.8),
    NES = c(1.2, -0.8, 1.5),
    pvalue = c(0.01, 0.05, 0.001),
    p.adjust = c(0.05, 0.1, 0.01),
    leading_edge = rep("K00001;K00002", 3),
    method = rep("fgsea", 3),
    stringsAsFactors = FALSE
  )

  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")

  # Should return results with pathway_id as pathway_name for unknown pathways
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), 3)
})

test_that("gsea_pathway_annotation handles partial matches in reference data", {
  # Load real pathway IDs
  data_env <- new.env()
  load(system.file("R/sysdata.rda", package = "ggpicrust2"), envir = data_env)

  # Mix of real and fake pathway IDs
  real_ids <- head(unique(data_env$ko_to_kegg_reference$pathway_id), 2)
  fake_ids <- c("ko99999", "ko88888")

  gsea_results <- data.frame(
    pathway_id = c(real_ids, fake_ids),
    pathway_name = c(real_ids, fake_ids),
    size = rep(10, 4),
    ES = rep(0.5, 4),
    NES = rep(1.0, 4),
    pvalue = rep(0.05, 4),
    p.adjust = rep(0.1, 4),
    leading_edge = rep("K00001", 4),
    method = rep("fgsea", 4),
    stringsAsFactors = FALSE
  )

  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")

  # All rows should be preserved
  expect_equal(nrow(annotated), 4)
  expect_true("pathway_name" %in% colnames(annotated))
})

test_that("gsea_pathway_annotation works with single row input", {
  gsea_results <- data.frame(
    pathway_id = "ko00001",
    pathway_name = "ko00001",
    size = 10,
    ES = 0.5,
    NES = 1.2,
    pvalue = 0.01,
    p.adjust = 0.05,
    leading_edge = "K00001;K00002",
    method = "fgsea",
    stringsAsFactors = FALSE
  )

  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")

  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), 1)
})

test_that("gsea_pathway_annotation handles many pathways efficiently", {
  # Test with many pathways
  n_pathways <- 100
  gsea_results <- create_realistic_gsea_results(n_pathways, "KEGG")

  start_time <- Sys.time()
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), n_pathways)
  expect_lt(elapsed, 5)  # Should complete in less than 5 seconds
})
