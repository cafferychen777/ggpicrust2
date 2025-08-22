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

# Helper function to create realistic KEGG reference data
create_mock_kegg_reference <- function(pathways = NULL) {
  if (is.null(pathways) || length(pathways) == 0) {
    # Return empty data frame with correct structure
    return(data.frame(
      pathway = character(0),
      pathway_name = character(0),
      pathway_class = character(0),
      description = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  pathway_classes <- c("Metabolism", "Genetic Information Processing", 
                      "Environmental Information Processing", "Cellular Processes",
                      "Human Diseases", "Drug Development")
  
  pathway_names <- c("Glycolysis", "TCA cycle", "Fatty acid biosynthesis", 
                    "Amino acid metabolism", "Nucleotide metabolism",
                    "Cell cycle", "DNA replication", "RNA transcription")
  
  data.frame(
    pathway = pathways,
    pathway_name = paste("KEGG", pathways, "-", 
                        sample(pathway_names, length(pathways), replace = TRUE)),
    pathway_class = sample(pathway_classes, length(pathways), replace = TRUE),
    description = paste("Detailed description for pathway", pathways),
    stringsAsFactors = FALSE
  )
}

# Helper function to create realistic MetaCyc reference data  
create_mock_metacyc_reference <- function(pathways = NULL) {
  if (is.null(pathways) || length(pathways) == 0) {
    return(data.frame(
      X1 = character(0),
      X2 = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # MetaCyc uses X1 and X2 column names based on the data structure we saw
  data.frame(
    X1 = pathways,
    X2 = paste(pathways, "pathway description", sep = " - "),
    stringsAsFactors = FALSE
  )
}

test_that("gsea_pathway_annotation handles KEGG pathways with full annotation", {
  # Create test data with specific pathway IDs
  test_pathways <- paste0("ko", sprintf("%05d", c(10, 20, 30, 40, 50)))
  gsea_results <- create_realistic_gsea_results(5, "KEGG")
  gsea_results$pathway_id <- test_pathways
  
  # Create mock KEGG reference data with all test pathways
  mock_kegg_ref <- create_mock_kegg_reference(test_pathways)
  
  # Use with_mocked_bindings (newer testthat approach)
  local_mocked_bindings(
    `file.exists` = function(...) TRUE,
    `load` = function(file) {
      # Simulate loading kegg_reference into the calling environment
      parent_env <- parent.frame()
      assign("kegg_reference", mock_kegg_ref, envir = parent_env)
    }
  )
  
  # Test annotation
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")
  
  # Check the result structure
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))
  
  # Check that all original columns are preserved
  original_cols <- colnames(gsea_results)
  expect_true(all(original_cols %in% colnames(annotated)))
  
  # Check that annotation columns were added
  expect_true("pathway_name" %in% colnames(annotated))
  expect_true("pathway_class" %in% colnames(annotated))
  expect_true("description" %in% colnames(annotated))
  
  # Check that pathway names were properly annotated
  expect_true(all(grepl("KEGG", annotated$pathway_name)))
  expect_true(all(!is.na(annotated$pathway_class)))
})

test_that("gsea_pathway_annotation handles KEGG pathways with missing reference file", {
  gsea_results <- create_realistic_gsea_results(5, "KEGG")
  
  # Mock file existence check to return FALSE
  local_mocked_bindings(
    `file.exists` = function(...) FALSE
  )
  
  expect_error(
    gsea_pathway_annotation(gsea_results, pathway_type = "KEGG"),
    "kegg_reference data file not found"
  )
})

test_that("gsea_pathway_annotation handles MetaCyc pathways correctly", {
  # Create test data
  test_pathways <- c("GLYCOLYSIS", "TCA-CYCLE", "PENTOSE-P-PWY")
  gsea_results <- create_realistic_gsea_results(3, "MetaCyc")
  gsea_results$pathway_id <- test_pathways
  
  # Create mock MetaCyc reference data
  mock_metacyc_ref <- create_mock_metacyc_reference(test_pathways)
  
  # Mock the data loading - MetaCyc uses the data() function
  local_mocked_bindings(
    `data` = function(x, package, envir) {
      if (x == "metacyc_reference") {
        # Note: the actual data uses "metacyc_reference" but file has "MetaCyc_reference"
        assign("metacyc_reference", mock_metacyc_ref, envir = envir)
      }
    }
  )
  
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "MetaCyc")
  
  # Check the result
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), nrow(gsea_results))
  expect_true("pathway_name" %in% colnames(annotated))
  
  # Check that pathway names were updated from MetaCyc reference
  expect_true(all(grepl("pathway description", annotated$pathway_name)))
})

test_that("gsea_pathway_annotation handles partial matches in reference data", {
  # Create test data with some pathways not in reference
  test_pathways <- paste0("ko", sprintf("%05d", c(10, 20, 30, 40, 50)))
  gsea_results <- create_realistic_gsea_results(5, "KEGG")
  gsea_results$pathway_id <- test_pathways
  
  # Create mock reference with only first 3 pathways
  mock_kegg_ref <- create_mock_kegg_reference(test_pathways[1:3])
  
  local_mocked_bindings(
    `file.exists` = function(...) TRUE,
    `load` = function(file) {
      parent_env <- parent.frame()
      assign("kegg_reference", mock_kegg_ref, envir = parent_env)
    }
  )
  
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")
  
  # Check that all rows are preserved
  expect_equal(nrow(annotated), 5)
  
  # Check that missing pathways use pathway_id as name
  # Pathways 4 and 5 should not be in reference and use pathway_id
  missing_pathways <- test_pathways[4:5]
  for (pathway in missing_pathways) {
    pathway_row <- annotated[annotated$pathway_id == pathway, ]
    expect_equal(pathway_row$pathway_name, pathway)
  }
  
  # Check that found pathways have proper names
  found_pathways <- test_pathways[1:3]
  for (pathway in found_pathways) {
    pathway_row <- annotated[annotated$pathway_id == pathway, ]
    expect_true(grepl("KEGG", pathway_row$pathway_name))
  }
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
  
  mock_kegg_ref <- create_mock_kegg_reference(character(0))
  
  local_mocked_bindings(
    `file.exists` = function(...) TRUE,
    `load` = function(file) {
      parent_env <- parent.frame()
      assign("kegg_reference", mock_kegg_ref, envir = parent_env)
    }
  )
  
  annotated <- gsea_pathway_annotation(empty_gsea, pathway_type = "KEGG")
  
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), 0)
  expect_true("pathway_name" %in% colnames(annotated))
})

test_that("gsea_pathway_annotation handles GO pathways (not implemented)", {
  gsea_results <- create_realistic_gsea_results(5, "GO")
  
  expect_warning(
    annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "GO"),
    "GO pathway annotation not yet implemented"
  )
  
  # Should return original results unchanged
  expect_equal(annotated, gsea_results)
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

test_that("gsea_pathway_annotation handles various data types and edge cases", {
  gsea_results <- create_realistic_gsea_results(5)
  
  # Test with factor pathway_ids
  gsea_results_factor <- gsea_results
  gsea_results_factor$pathway_id <- as.factor(gsea_results_factor$pathway_id)
  
  mock_kegg_ref <- create_mock_kegg_reference(as.character(gsea_results_factor$pathway_id))
  
  local_mocked_bindings(
    `file.exists` = function(...) TRUE,
    `load` = function(file) {
      parent_env <- parent.frame()
      assign("kegg_reference", mock_kegg_ref, envir = parent_env)
    }
  )
  
  annotated <- gsea_pathway_annotation(gsea_results_factor, pathway_type = "KEGG")
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), 5)
  
  # Test with NA values in pathway_id
  gsea_results_na <- gsea_results
  gsea_results_na$pathway_id[1] <- NA
  
  local_mocked_bindings(
    `file.exists` = function(...) TRUE,
    `load` = function(file) {
      parent_env <- parent.frame()
      assign("kegg_reference", mock_kegg_ref, envir = parent_env)
    }
  )
  
  annotated <- gsea_pathway_annotation(gsea_results_na, pathway_type = "KEGG")
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), 5)
})

test_that("gsea_pathway_annotation preserves original column order and types", {
  gsea_results <- create_realistic_gsea_results(3)
  original_names <- colnames(gsea_results)
  original_classes <- sapply(gsea_results, class)
  
  mock_kegg_ref <- create_mock_kegg_reference(gsea_results$pathway_id)
  
  local_mocked_bindings(
    `file.exists` = function(...) TRUE,
    `load` = function(file) {
      parent_env <- parent.frame()
      assign("kegg_reference", mock_kegg_ref, envir = parent_env)
    }
  )
  
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")
  
  # Check that all original columns are preserved
  expect_true(all(original_names %in% colnames(annotated)))
  
  # Check that original column types are preserved (for existing columns)
  for (col in original_names) {
    if (col != "pathway_name") {  # pathway_name gets replaced
      expect_equal(class(annotated[[col]]), original_classes[[col]])
    }
  }
})

test_that("gsea_pathway_annotation handles reference data conversion correctly", {
  # Test that the function properly converts reference data to data frame
  gsea_results <- create_realistic_gsea_results(3, "KEGG")
  
  # Create mock reference as matrix (should be converted to data frame)
  mock_kegg_matrix <- as.matrix(create_mock_kegg_reference(gsea_results$pathway_id))
  
  local_mocked_bindings(
    `file.exists` = function(...) TRUE,
    `load` = function(file) {
      parent_env <- parent.frame()
      assign("kegg_reference", mock_kegg_matrix, envir = parent_env)
    }
  )
  
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")
  
  expect_s3_class(annotated, "data.frame")
  expect_equal(nrow(annotated), 3)
})

test_that("gsea_pathway_annotation handles missing pathway_name column in reference", {
  # Test when reference data doesn't have pathway_name column
  gsea_results <- create_realistic_gsea_results(3, "KEGG")
  
  # Create reference without pathway_name column
  mock_kegg_ref <- data.frame(
    pathway = gsea_results$pathway_id,
    description = paste("Description for", gsea_results$pathway_id),
    stringsAsFactors = FALSE
  )
  
  local_mocked_bindings(
    `file.exists` = function(...) TRUE,
    `load` = function(file) {
      parent_env <- parent.frame()
      assign("kegg_reference", mock_kegg_ref, envir = parent_env)
    }
  )
  
  annotated <- gsea_pathway_annotation(gsea_results, pathway_type = "KEGG")
  
  expect_s3_class(annotated, "data.frame")
  expect_true("pathway_name" %in% colnames(annotated))
  # Should use pathway_id as pathway_name
  expect_equal(annotated$pathway_name, annotated$pathway_id)
})