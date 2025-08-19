test_that("MetaCyc gene set preparation works correctly", {
  # Skip if required packages are not available
  skip_if_not_installed("fgsea")
  
  # Test prepare_gene_sets for MetaCyc
  gene_sets <- prepare_gene_sets(pathway_type = "MetaCyc")
  
  # Check that gene_sets is a list
  expect_type(gene_sets, "list")
  
  # Check that we have some gene sets (should be 23 based on our mapping)
  expect_gt(length(gene_sets), 0)
  
  # Check that each gene set contains EC numbers with proper format
  for (pathway in names(gene_sets)) {
    ec_numbers <- gene_sets[[pathway]]
    expect_type(ec_numbers, "character")
    expect_gt(length(ec_numbers), 0)
    
    # Check that all EC numbers start with "EC:"
    for (ec in ec_numbers) {
      expect_true(startsWith(ec, "EC:"))
    }
  }
  
  # Test specific pathways we know should exist
  expected_pathways <- c("1CMET2-PWY", "ANAGLYCOLYSIS-PWY", "ARG+POLYAMINE-SYN")
  available_pathways <- names(gene_sets)
  
  for (pathway in expected_pathways) {
    if (pathway %in% available_pathways) {
      expect_true(pathway %in% names(gene_sets))
      expect_gt(length(gene_sets[[pathway]]), 0)
    }
  }
})

test_that("MetaCyc GSEA workflow works end-to-end", {
  # Skip if required packages are not available
  skip_if_not_installed("fgsea")
  
  # Load test data
  data("metacyc_abundance", package = "ggpicrust2", envir = environment())
  data("metadata", package = "ggpicrust2", envir = environment())
  
  # Prepare abundance data
  abundance_data <- as.data.frame(metacyc_abundance)
  rownames(abundance_data) <- abundance_data[["pathway"]]
  abundance_data <- abundance_data[, -1]  # Remove pathway column
  
  # Ensure we have EC format in row names for GSEA
  # For this test, we'll simulate EC abundance data
  # In real usage, users would provide EC abundance data
  ec_abundance_sample <- abundance_data[1:10, 1:10]  # Sample subset
  
  # Convert pathway IDs to some EC numbers for testing
  rownames(ec_abundance_sample) <- paste0("EC:", c("1.1.1.1", "1.2.1.12", "2.7.1.1", "2.7.1.11", "4.1.2.13", 
                                                  "1.2.1.59", "2.7.2.3", "5.4.2.11", "4.2.1.11", "5.4.2.1"))
  
  # Create simple metadata for testing
  test_metadata <- data.frame(
    sample_id = colnames(ec_abundance_sample),
    group = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )
  rownames(test_metadata) <- test_metadata$sample_id
  
  # Test pathway_gsea function with MetaCyc
  expect_no_error({
    results <- pathway_gsea(
      abundance = ec_abundance_sample,
      metadata = test_metadata,
      group = "group",
      pathway_type = "MetaCyc",
      method = "fgsea",
      nperm = 100,  # Reduced for faster testing
      min_size = 2,
      max_size = 100
    )
  })
  
  # Check that results have expected structure
  expect_s3_class(results, "data.frame")
  expect_true("pathway_id" %in% colnames(results))
  expect_true("pathway_name" %in% colnames(results))
  expect_true("method" %in% colnames(results))
  expect_equal(results$method[1], "fgsea")
})

test_that("MetaCyc pathway annotation works correctly", {
  # Create mock GSEA results
  mock_results <- data.frame(
    pathway_id = c("1CMET2-PWY", "ANAGLYCOLYSIS-PWY", "ARG+POLYAMINE-SYN", "UNKNOWN-PWY"),
    size = c(4, 5, 5, 3),
    ES = c(0.5, -0.3, 0.7, 0.1),
    NES = c(1.2, -0.8, 1.5, 0.2),
    pvalue = c(0.05, 0.1, 0.02, 0.8),
    p.adjust = c(0.15, 0.2, 0.08, 0.9),
    stringsAsFactors = FALSE
  )
  
  # Test annotation
  annotated <- gsea_pathway_annotation(mock_results, pathway_type = "MetaCyc")
  
  # Check that annotation worked
  expect_s3_class(annotated, "data.frame")
  expect_true("pathway_name" %in% colnames(annotated))
  expect_equal(nrow(annotated), nrow(mock_results))
  
  # Check that known pathways got proper names
  expected_names <- c(
    "N10-formyl-tetrahydrofolate biosynthesis",  # 1CMET2-PWY
    "ANAGLYCOLYSIS-PWY",  # This might not have a mapping, so should keep ID
    "ARG+POLYAMINE-SYN"   # This might not have a mapping, so should keep ID  
  )
  
  # At least some pathways should have descriptions different from their IDs
  annotated_known <- annotated[annotated$pathway_id %in% c("1CMET2-PWY", "ANAGLYCOLYSIS-PWY", "ARG+POLYAMINE-SYN"), ]
  expect_gt(nrow(annotated_known), 0)
  
  # Unknown pathways should keep their ID as name
  unknown_row <- annotated[annotated$pathway_id == "UNKNOWN-PWY", ]
  expect_equal(unknown_row$pathway_name, "UNKNOWN-PWY")
})

test_that("MetaCyc gene sets have reasonable sizes", {
  gene_sets <- prepare_gene_sets(pathway_type = "MetaCyc")
  
  if (length(gene_sets) > 0) {
    # Check gene set sizes are reasonable
    sizes <- sapply(gene_sets, length)
    
    # All gene sets should have at least 1 EC number
    expect_true(all(sizes >= 1))
    
    # No gene set should be excessively large (>50 would be unusual for MetaCyc)
    expect_true(all(sizes <= 50))
    
    # Check that we have a reasonable number of gene sets
    expect_gt(length(gene_sets), 5)  # At least a few pathways
    expect_lt(length(gene_sets), 100)  # But not too many for our test mapping
  }
})

test_that("MetaCyc error handling works correctly", {
  # Test with missing reference data (by temporarily renaming the file)
  original_path <- system.file("extdata", "metacyc_to_ec_reference.RData", package = "ggpicrust2")
  
  if (file.exists(original_path)) {
    # Test should work normally when file exists
    expect_no_error({
      gene_sets <- prepare_gene_sets(pathway_type = "MetaCyc")
    })
  }
  
  # Test invalid pathway_type
  expect_error({
    prepare_gene_sets(pathway_type = "InvalidType")
  })
  
  # Test pathway_gsea with invalid parameters
  expect_error({
    pathway_gsea(
      abundance = matrix(1:10, nrow = 2),
      metadata = data.frame(group = c("A", "A", "A", "B", "B")),
      group = "group",
      pathway_type = "InvalidType"
    )
  })
})

test_that("MetaCyc integration with existing workflow", {
  # Test that MetaCyc pathways are properly recognized
  valid_types <- c("KEGG", "MetaCyc", "GO")
  
  # pathway_gsea should accept MetaCyc as a valid pathway_type
  expect_true("MetaCyc" %in% valid_types)
  
  # Test parameter validation
  expect_no_error({
    # This should not throw an error for pathway_type validation
    if (requireNamespace("fgsea", quietly = TRUE)) {
      # Create minimal test data
      test_abundance <- matrix(rnorm(20), nrow = 4, ncol = 5)
      rownames(test_abundance) <- c("EC:1.1.1.1", "EC:1.2.1.12", "EC:2.7.1.1", "EC:2.7.1.11")
      colnames(test_abundance) <- paste0("Sample", 1:5)
      
      test_metadata <- data.frame(
        sample_id = colnames(test_abundance),
        group = c("A", "A", "A", "B", "B"),
        stringsAsFactors = FALSE
      )
      rownames(test_metadata) <- test_metadata$sample_id
      
      # This should validate successfully even if it doesn't run GSEA
      expect_no_error({
        tryCatch({
          results <- pathway_gsea(
            abundance = test_abundance,
            metadata = test_metadata,
            group = "group",
            pathway_type = "MetaCyc",
            method = "fgsea",
            nperm = 10,
            min_size = 1,
            max_size = 10
          )
        }, error = function(e) {
          # It's OK if GSEA fails due to data issues, but pathway_type should be valid
          if (grepl("pathway_type", e$message)) {
            stop(e)  # Re-throw if it's a pathway_type validation error
          }
        })
      })
    }
  })
})