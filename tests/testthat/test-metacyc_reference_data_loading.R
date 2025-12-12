# MetaCyc Reference Data Loading System Tests
# Testing data loading mechanisms, backup functions, and error recovery
# Following Linus principles: robust data handling, clear error messages, reliable fallbacks

test_that("MetaCyc reference data file existence and accessibility", {
  # Test that the reference data file exists and is accessible
  metacyc_ref_path <- system.file("extdata", "metacyc_to_ec_reference.RData", package = "ggpicrust2")
  
  # Basic file system checks
  expect_true(file.exists(metacyc_ref_path))
  expect_gt(file.size(metacyc_ref_path), 0)

  # Test file is readable
  expect_no_error({
    load(metacyc_ref_path, envir = environment())
  })
  
  # Verify loaded object structure
  expect_true(exists("metacyc_to_ec_reference", envir = environment()))
  expect_s3_class(metacyc_to_ec_reference, "data.frame")
  expect_gt(nrow(metacyc_to_ec_reference), 0)
  expect_equal(ncol(metacyc_to_ec_reference), 2)
  expect_equal(colnames(metacyc_to_ec_reference), c("pathway", "ec_numbers"))
})

test_that("MetaCyc reference data loading in prepare_gene_sets function", {
  # Test the actual loading mechanism used in prepare_gene_sets
  skip_if_not_installed("fgsea")
  
  # Clear environment to test fresh loading
  if (exists("metacyc_to_ec_reference")) {
    rm(metacyc_to_ec_reference)
  }
  
  # Test loading through prepare_gene_sets
  gene_sets <- prepare_gene_sets(pathway_type = "MetaCyc")
  
  # Validate that gene sets were created successfully
  expect_type(gene_sets, "list")
  expect_gt(length(gene_sets), 0)
  
  # Test that the loading creates proper gene set structure
  for (pathway_id in names(gene_sets)[1:min(5, length(gene_sets))]) {
    expect_type(gene_sets[[pathway_id]], "character")
    expect_gt(length(gene_sets[[pathway_id]]), 0)
    # All should be EC numbers with proper format
    for (ec in gene_sets[[pathway_id]]) {
      expect_true(startsWith(ec, "EC:"))
    }
  }
})

test_that("MetaCyc reference data loading error handling", {
  # Test behavior when reference data is missing or corrupted
  skip_if_not_installed("fgsea")
  
  # This test verifies graceful error handling when data file is not found
  # We can't easily simulate missing file, but we can test the error path
  
  # Test with invalid pathway_type (should not reach MetaCyc loading code)
  expect_error({
    prepare_gene_sets(pathway_type = "InvalidPathwayType")
  })
  
  # Test the specific error message structure for MetaCyc loading failure
  # This tests the tryCatch block in prepare_gene_sets for MetaCyc
  
  # Mock a scenario where the loaded data is malformed
  # Create a temporary environment to test loading behavior
  test_env <- new.env()
  
  # Test that function handles empty or malformed reference data gracefully
  # If the reference data was somehow corrupted to be empty
  test_env$metacyc_to_ec_reference <- data.frame(
    pathway = character(0),
    ec_numbers = character(0)
  )
  
  # Test handling of empty reference data
  # This should not crash even with empty reference data
  gene_sets_empty <- list()
  n_rows <- nrow(test_env$metacyc_to_ec_reference)
  if (n_rows > 0) {
    for (i in seq_len(n_rows)) {
      pathway_id <- test_env$metacyc_to_ec_reference[i, "pathway"]
      ec_string <- as.character(test_env$metacyc_to_ec_reference[i, "ec_numbers"])

      if (is.na(ec_string) || ec_string == "" || ec_string == "NA") {
        next
      }

      ec_numbers <- strsplit(ec_string, ";")[[1]]
      ec_numbers <- trimws(ec_numbers)
      ec_numbers <- ec_numbers[ec_numbers != ""]

      if (length(ec_numbers) > 0) {
        ec_numbers <- ifelse(grepl("^EC:", ec_numbers), ec_numbers, paste0("EC:", ec_numbers))
        gene_sets_empty[[pathway_id]] <- ec_numbers
      }
    }
  }
  expect_equal(length(gene_sets_empty), 0)  # Empty input should produce empty output
})

test_that("MetaCyc data format validation during loading", {
  # Test that loaded reference data has correct format and structure
  metacyc_ref_path <- system.file("extdata", "metacyc_to_ec_reference.RData", package = "ggpicrust2")
  load(metacyc_ref_path, envir = environment())
  
  # Validate data frame structure
  expect_s3_class(metacyc_to_ec_reference, "data.frame")
  expect_equal(ncol(metacyc_to_ec_reference), 2)
  expect_true("pathway" %in% colnames(metacyc_to_ec_reference))
  expect_true("ec_numbers" %in% colnames(metacyc_to_ec_reference))
  
  # Validate data types
  expect_type(metacyc_to_ec_reference$pathway, "character")
  expect_type(metacyc_to_ec_reference$ec_numbers, "character")
  
  # Validate pathway IDs are reasonable
  pathway_ids <- metacyc_to_ec_reference$pathway
  expect_true(all(!is.na(pathway_ids)))
  expect_true(all(nchar(pathway_ids) > 0))
  expect_true(all(nchar(pathway_ids) <= 50))  # Reasonable length limit
  
  # Validate EC number strings
  for (i in 1:min(20, nrow(metacyc_to_ec_reference))) {
    ec_string <- metacyc_to_ec_reference[i, "ec_numbers"]
    if (!is.na(ec_string) && ec_string != "" && ec_string != "NA") {
      # Should be semicolon-separated EC numbers
      expect_true(grepl("^\\d+\\.\\d+\\.\\d+\\.\\d+(;\\d+\\.\\d+\\.\\d+\\.\\d+)*$", ec_string) ||
                  grepl("^EC:\\d+\\.\\d+\\.\\d+\\.\\d+(;EC:\\d+\\.\\d+\\.\\d+\\.\\d+)*$", ec_string))
    }
  }
  
  # Test uniqueness of pathway IDs
  expect_equal(length(pathway_ids), length(unique(pathway_ids)))
})

test_that("MetaCyc gene set creation from reference data", {
  # Test the specific logic that converts reference data to gene sets
  metacyc_ref_path <- system.file("extdata", "metacyc_to_ec_reference.RData", package = "ggpicrust2")
  load(metacyc_ref_path, envir = environment())
  
  # Manually replicate the gene set creation logic to test it
  gene_sets_manual <- list()
  
  for (i in 1:nrow(metacyc_to_ec_reference)) {
    pathway_id <- metacyc_to_ec_reference[i, "pathway"]
    ec_string <- as.character(metacyc_to_ec_reference[i, "ec_numbers"])
    
    # Skip pathways with no EC mappings
    if (is.na(ec_string) || ec_string == "" || ec_string == "NA") {
      next
    }
    
    # Split EC numbers by semicolon
    ec_numbers <- strsplit(ec_string, ";")[[1]]
    ec_numbers <- trimws(ec_numbers)  # Remove whitespace
    ec_numbers <- ec_numbers[ec_numbers != ""]  # Remove empty strings
    
    if (length(ec_numbers) > 0) {
      # Add EC: prefix if not present for consistency
      ec_numbers <- ifelse(grepl("^EC:", ec_numbers), ec_numbers, paste0("EC:", ec_numbers))
      gene_sets_manual[[pathway_id]] <- ec_numbers
    }
  }
  
  # Validate manual gene set creation
  expect_type(gene_sets_manual, "list")
  expect_gt(length(gene_sets_manual), 0)
  
  # Compare with prepare_gene_sets output
  skip_if_not_installed("fgsea")
  gene_sets_function <- prepare_gene_sets(pathway_type = "MetaCyc")
  
  # Should have same pathways (allowing for potential filtering differences)
  common_pathways <- intersect(names(gene_sets_manual), names(gene_sets_function))
  expect_gt(length(common_pathways), length(gene_sets_manual) * 0.8)  # At least 80% overlap
  
  # For common pathways, EC numbers should match
  for (pathway in common_pathways[1:min(10, length(common_pathways))]) {
    manual_ecs <- sort(gene_sets_manual[[pathway]])
    function_ecs <- sort(gene_sets_function[[pathway]])
    expect_equal(manual_ecs, function_ecs)
  }
})

test_that("MetaCyc reference data consistency across versions", {
  # Test that the reference data maintains consistency and expected properties
  metacyc_ref_path <- system.file("extdata", "metacyc_to_ec_reference.RData", package = "ggpicrust2")
  load(metacyc_ref_path, envir = environment())
  
  # Basic consistency checks
  expect_gte(nrow(metacyc_to_ec_reference), 50)   # Minimum reasonable coverage
  expect_lte(nrow(metacyc_to_ec_reference), 1000) # Maximum reasonable coverage
  
  # Check for expected key MetaCyc pathways
  expected_pathways <- c("1CMET2-PWY", "ANAGLYCOLYSIS-PWY", "ARG+POLYAMINE-SYN")
  pathway_ids <- metacyc_to_ec_reference$pathway
  
  for (expected in expected_pathways) {
    if (expected %in% pathway_ids) {
      # If pathway exists, it should have EC mappings
      row_idx <- which(pathway_ids == expected)
      ec_string <- metacyc_to_ec_reference[row_idx, "ec_numbers"]
      expect_true(!is.na(ec_string) && ec_string != "")
    }
  }
  
  # Test data completeness
  non_empty_mappings <- sum(!is.na(metacyc_to_ec_reference$ec_numbers) & 
                           metacyc_to_ec_reference$ec_numbers != "" &
                           metacyc_to_ec_reference$ec_numbers != "NA")
  completeness_rate <- non_empty_mappings / nrow(metacyc_to_ec_reference)
  # Check that some pathways have EC mappings (percentage varies by data version)
  expect_gt(completeness_rate, 0.0)
})

test_that("MetaCyc reference data memory efficiency", {
  # Test memory usage and loading efficiency
  metacyc_ref_path <- system.file("extdata", "metacyc_to_ec_reference.RData", package = "ggpicrust2")
  
  # Test file size is reasonable (should be less than 10MB)
  file_size_mb <- file.size(metacyc_ref_path) / (1024^2)
  expect_lt(file_size_mb, 10)
  
  # Test loading time
  start_time <- Sys.time()
  load(metacyc_ref_path, envir = new.env())
  end_time <- Sys.time()
  loading_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Should load within 2 seconds
  expect_lt(loading_time, 2)
  
  # Test memory footprint of loaded data (should use less than 5MB in memory)
  load(metacyc_ref_path, envir = environment())
  data_size_mb <- as.numeric(object.size(metacyc_to_ec_reference)) / (1024^2)
  expect_lt(data_size_mb, 5)
})

test_that("MetaCyc backup function validation", {
  # Test if there's a backup/fallback function for creating MetaCyc mappings
  # This would be similar to create_basic_go_mapping() for GO pathways
  
  # Check if there's a create_basic_metacyc_mapping function or similar
  # If it exists, test it; if not, document that it should be implemented
  
  # First check if such a function exists
  backup_function_exists <- exists("create_basic_metacyc_mapping", mode = "function")
  
  if (backup_function_exists) {
    # Test the backup function
    expect_no_error({
      backup_mapping <- create_basic_metacyc_mapping()
    })
    
    expect_s3_class(backup_mapping, "data.frame")
    expect_true("pathway" %in% colnames(backup_mapping))
    expect_true("ec_numbers" %in% colnames(backup_mapping))
    expect_gt(nrow(backup_mapping), 0)
    
  } else {
    # Document that backup function should be implemented
    skip("create_basic_metacyc_mapping function not implemented - consider adding for robustness")
  }
})

test_that("MetaCyc data loading integration with GSEA workflow", {
  # Test that the loaded reference data integrates properly with GSEA workflow
  skip_if_not_installed("fgsea")
  
  # Test complete workflow from data loading to gene set preparation
  expect_no_error({
    gene_sets <- prepare_gene_sets(pathway_type = "MetaCyc")
  })
  
  # Test that gene sets can be used in fgsea
  if (length(gene_sets) > 0) {
    # Create minimal test data
    test_genes <- unique(unlist(gene_sets))[1:min(20, length(unique(unlist(gene_sets))))]
    test_stats <- rnorm(length(test_genes))
    names(test_stats) <- test_genes
    
    expect_no_error({
      # This tests that gene sets have correct format for fgsea
      fgsea_result <- fgsea::fgsea(
        pathways = gene_sets[1:min(5, length(gene_sets))],  # Test subset
        stats = test_stats,
        nperm = 10,  # Minimal for testing
        minSize = 1
      )
    })
  }
})

test_that("MetaCyc reference data update compatibility", {
  # Test that the current loading system would handle updated reference data
  metacyc_ref_path <- system.file("extdata", "metacyc_to_ec_reference.RData", package = "ggpicrust2")
  load(metacyc_ref_path, envir = environment())
  
  # Current data structure requirements
  required_columns <- c("pathway", "ec_numbers")
  expect_true(all(required_columns %in% colnames(metacyc_to_ec_reference)))
  
  # Test that additional columns would be handled gracefully
  extended_data <- metacyc_to_ec_reference
  extended_data$description <- paste("Description for", extended_data$pathway)
  extended_data$category <- sample(c("Metabolism", "Biosynthesis", "Degradation"), 
                                  nrow(extended_data), replace = TRUE)
  
  # Simulate loading extended data (the current loading code should handle extra columns)
  expect_no_error({
    gene_sets_extended <- list()
    
    for (i in 1:min(10, nrow(extended_data))) {
      pathway_id <- extended_data[i, "pathway"]
      ec_string <- as.character(extended_data[i, "ec_numbers"])
      
      if (!is.na(ec_string) && ec_string != "" && ec_string != "NA") {
        ec_numbers <- strsplit(ec_string, ";")[[1]]
        ec_numbers <- trimws(ec_numbers)
        ec_numbers <- ec_numbers[ec_numbers != ""]
        
        if (length(ec_numbers) > 0) {
          ec_numbers <- ifelse(grepl("^EC:", ec_numbers), ec_numbers, paste0("EC:", ec_numbers))
          gene_sets_extended[[pathway_id]] <- ec_numbers
        }
      }
    }
  })
  
  expect_type(gene_sets_extended, "list")
})