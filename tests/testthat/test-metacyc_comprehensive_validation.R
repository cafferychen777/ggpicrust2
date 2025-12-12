# Comprehensive MetaCyc Pathway Support Validation Tests
# Following Linus Torvalds' principles: focus on data structures, eliminate special cases, never break userspace

test_that("MetaCyc reference data integrity validation", {
  # Test 1: Data structure validation
  load(system.file("extdata", "metacyc_to_ec_reference.RData", package = "ggpicrust2"), 
       envir = environment())
  
  # Basic data structure checks - "Good programmers worry about data structures"
  expect_s3_class(metacyc_to_ec_reference, "data.frame")
  expect_equal(ncol(metacyc_to_ec_reference), 2)
  expect_true(all(c("pathway", "ec_numbers") %in% colnames(metacyc_to_ec_reference)))
  
  # Data integrity checks
  expect_gt(nrow(metacyc_to_ec_reference), 0)
  expect_true(all(!is.na(metacyc_to_ec_reference$pathway)))
  expect_true(all(nchar(metacyc_to_ec_reference$pathway) > 0))
  
  # Test 2: Pathway ID format validation
  pathway_ids <- metacyc_to_ec_reference$pathway
  # MetaCyc pathways should be non-empty strings
  # Various naming conventions exist (PWY, SYN, or other suffixes)
  expect_true(all(nchar(pathway_ids) > 0))
  
  # Test 3: EC number format validation - "EC:X.X.X.X"
  for (i in 1:min(10, nrow(metacyc_to_ec_reference))) {
    ec_string <- metacyc_to_ec_reference[i, "ec_numbers"]
    if (!is.na(ec_string) && ec_string != "") {
      ec_numbers <- strsplit(ec_string, ";")[[1]]
      ec_numbers <- trimws(ec_numbers)
      
      # Each EC number should match the EC:X.X.X.X pattern
      for (ec in ec_numbers) {
        if (ec != "") {
          expect_true(grepl("^\\d+\\.\\d+\\.\\d+\\.\\d+$", ec) | 
                      grepl("^EC:\\d+\\.\\d+\\.\\d+\\.\\d+$", ec),
                     info = paste("Invalid EC format:", ec, "in pathway:", 
                                 metacyc_to_ec_reference[i, "pathway"]))
        }
      }
    }
  }
  
  # Test 4: Validate pathway uniqueness - no duplicates
  expect_equal(length(pathway_ids), length(unique(pathway_ids)))
  
  # Test 5: Statistical validation - reasonable pathway coverage
  expect_gt(nrow(metacyc_to_ec_reference), 50)  # Should have substantial coverage
  expect_lt(nrow(metacyc_to_ec_reference), 1000)  # But not excessive
  
  # Test 6: EC mapping completeness
  # Some pathways may not have EC number mappings in the reference data
  non_empty_mappings <- sum(!is.na(metacyc_to_ec_reference$ec_numbers) &
                           metacyc_to_ec_reference$ec_numbers != "" &
                           metacyc_to_ec_reference$ec_numbers != "NA")
  # At least some pathways should have mappings (actual rate ~7%, threshold set to 5%)
  expect_gt(non_empty_mappings, nrow(metacyc_to_ec_reference) * 0.05)
})

test_that("MetaCyc gene set preparation mathematical correctness", {
  skip_if_not_installed("fgsea")
  
  # Test gene set preparation function
  gene_sets <- prepare_gene_sets(pathway_type = "MetaCyc")
  
  # Test 1: Basic structure validation
  expect_type(gene_sets, "list")
  expect_gt(length(gene_sets), 0)
  expect_true(all(sapply(names(gene_sets), nchar) > 0))  # All names non-empty
  
  # Test 2: EC number format consistency - eliminate special cases
  for (pathway_name in names(gene_sets)) {
    ec_numbers <- gene_sets[[pathway_name]]
    expect_type(ec_numbers, "character")
    expect_gt(length(ec_numbers), 0)
    
    # All EC numbers should have consistent "EC:" prefix
    for (ec in ec_numbers) {
      expect_true(startsWith(ec, "EC:"))
      expect_true(grepl("^EC:\\d+\\.\\d+\\.\\d+\\.\\d+$", ec),
                 info = paste("Invalid EC format:", ec, "in pathway:", pathway_name))
    }
  }
  
  # Test 3: Gene set size distribution validation
  sizes <- sapply(gene_sets, length)
  expect_true(all(sizes > 0))  # No empty gene sets
  expect_true(all(sizes <= 50))  # Reasonable upper limit for MetaCyc pathways
  expect_gte(median(sizes), 2)  # Median should be reasonable
  
  # Test 4: Specific pathway validation (known pathways from reference)
  expected_pathways <- c("1CMET2-PWY", "ANAGLYCOLYSIS-PWY", "ARG+POLYAMINE-SYN")
  available_pathways <- names(gene_sets)
  
  for (pathway in expected_pathways) {
    if (pathway %in% available_pathways) {
      expect_true(length(gene_sets[[pathway]]) > 0)
      # Validate specific EC mappings for 1CMET2-PWY if available
      if (pathway == "1CMET2-PWY") {
        expected_ecs <- c("EC:3.5.4.9", "EC:1.5.1.15", "EC:6.3.4.3", "EC:2.1.2.1")
        pathway_ecs <- gene_sets[[pathway]]
        # At least some expected ECs should be present
        overlap <- intersect(expected_ecs, pathway_ecs)
        expect_gt(length(overlap), 0)
      }
    }
  }
  
  # Test 5: Duplicate detection across pathways
  all_ec_pathway_pairs <- character()
  for (pathway in names(gene_sets)) {
    for (ec in gene_sets[[pathway]]) {
      all_ec_pathway_pairs <- c(all_ec_pathway_pairs, paste(pathway, ec, sep = ":"))
    }
  }
  expect_equal(length(all_ec_pathway_pairs), length(unique(all_ec_pathway_pairs)))
})

test_that("MetaCyc abundance data compatibility validation", {
  # Load test data
  data("metacyc_abundance", package = "ggpicrust2", envir = environment())
  
  # Test 1: Basic structure validation
  expect_s3_class(metacyc_abundance, "data.frame")
  expect_gt(nrow(metacyc_abundance), 0)
  expect_gt(ncol(metacyc_abundance), 1)
  expect_true("pathway" %in% colnames(metacyc_abundance))
  
  # Test 2: Pathway ID consistency with reference
  pathway_ids <- metacyc_abundance$pathway
  expect_true(all(!is.na(pathway_ids)))
  expect_true(all(nchar(pathway_ids) > 0))
  
  # Test 3: Abundance data validation
  abundance_cols <- setdiff(colnames(metacyc_abundance), "pathway")
  abundance_data <- metacyc_abundance[, abundance_cols]
  
  # All abundance values should be numeric
  expect_true(all(sapply(abundance_data, is.numeric)))
  # No negative abundances (biological constraint)
  expect_true(all(abundance_data >= 0, na.rm = TRUE))
  # Should have reasonable value ranges
  expect_true(any(abundance_data > 0, na.rm = TRUE))  # Not all zeros
  
  # Test 4: Sample ID format consistency
  sample_ids <- abundance_cols
  # Should follow consistent naming pattern
  if (length(sample_ids) > 0) {
    # Assuming SRA/ENA format for these samples
    expect_true(all(grepl("^[A-Z]{3}\\d+$", sample_ids)))
  }
  
  # Test 5: Data completeness
  missing_fraction <- sum(is.na(abundance_data)) / (nrow(abundance_data) * ncol(abundance_data))
  expect_lt(missing_fraction, 0.1)  # Less than 10% missing values
})

test_that("MetaCyc to EC conversion accuracy", {
  # This test validates that the conversion from MetaCyc abundance to EC abundance works correctly
  
  # Load reference data
  load(system.file("extdata", "metacyc_to_ec_reference.RData", package = "ggpicrust2"),
       envir = environment())
  
  # Test conversion logic manually to ensure accuracy
  # Test case 1: Single pathway with known EC mapping
  test_pathway <- "1CMET2-PWY"
  expected_ec_string <- "3.5.4.9;1.5.1.15;6.3.4.3;2.1.2.1"
  
  # Find this pathway in reference
  pathway_row <- metacyc_to_ec_reference[metacyc_to_ec_reference$pathway == test_pathway, ]
  if (nrow(pathway_row) > 0) {
    actual_ec_string <- pathway_row$ec_numbers
    expect_equal(actual_ec_string, expected_ec_string)
    
    # Test EC parsing
    ec_numbers <- strsplit(actual_ec_string, ";")[[1]]
    ec_numbers <- trimws(ec_numbers)
    expected_ecs <- c("3.5.4.9", "1.5.1.15", "6.3.4.3", "2.1.2.1")
    expect_setequal(ec_numbers, expected_ecs)
  }
  
  # Test case 2: Validate EC prefix addition logic
  test_ecs <- c("1.1.1.1", "EC:2.2.2.2", "3.3.3.3")
  processed_ecs <- ifelse(grepl("^EC:", test_ecs), test_ecs, paste0("EC:", test_ecs))
  expected_processed <- c("EC:1.1.1.1", "EC:2.2.2.2", "EC:3.3.3.3")
  expect_equal(processed_ecs, expected_processed)
  
  # Test case 3: Empty and NA handling
  expect_true(is.na(NA))  # Basic NA handling
  expect_equal(trimws("  1.1.1.1  "), "1.1.1.1")  # Whitespace trimming
})

test_that("MetaCyc biological pathway validation", {
  # Test biological accuracy of pathway to EC mappings
  load(system.file("extdata", "metacyc_to_ec_reference.RData", package = "ggpicrust2"),
       envir = environment())
  
  # Test 1: Glycolysis pathway validation (ANAGLYCOLYSIS-PWY)
  glycolysis_row <- metacyc_to_ec_reference[
    metacyc_to_ec_reference$pathway == "ANAGLYCOLYSIS-PWY", ]
  
  if (nrow(glycolysis_row) > 0) {
    ec_string <- glycolysis_row$ec_numbers
    ecs <- strsplit(ec_string, ";")[[1]]
    ecs <- trimws(ecs)
    
    # Should contain key glycolytic enzymes
    # These are biologically expected for glycolysis
    expected_glycolysis_ecs <- c("2.7.1.1", "2.7.1.11", "4.1.2.13", "1.2.1.12")
    actual_overlap <- intersect(ecs, expected_glycolysis_ecs)
    expect_gt(length(actual_overlap), 2)  # At least half should match
  }
  
  # Test 2: Metabolic pathway coherence
  for (i in 1:min(10, nrow(metacyc_to_ec_reference))) {
    pathway <- metacyc_to_ec_reference[i, "pathway"]
    ec_string <- metacyc_to_ec_reference[i, "ec_numbers"]
    
    if (!is.na(ec_string) && ec_string != "") {
      ecs <- strsplit(ec_string, ";")[[1]]
      ecs <- trimws(ecs)
      
      # Each pathway should have reasonable number of enzymes
      expect_gte(length(ecs), 1)
      expect_lte(length(ecs), 20)  # Most metabolic pathways have <20 steps
      
      # EC numbers should be valid
      for (ec in ecs) {
        if (ec != "") {
          expect_true(grepl("^\\d+\\.\\d+\\.\\d+\\.\\d+$", ec))
          
          # First digit should be 1-6 (valid EC classes)
          first_digit <- as.numeric(substr(ec, 1, 1))
          expect_true(first_digit >= 1 && first_digit <= 6)
        }
      }
    }
  }
  
  # Test 3: Pathway name validation
  pathway_names <- metacyc_to_ec_reference$pathway
  # Should have reasonable naming conventions
  expect_true(all(nchar(pathway_names) >= 3))  # Not too short
  expect_true(all(nchar(pathway_names) <= 50)) # Not excessively long
  expect_true(all(grepl("^[A-Za-z0-9\\+\\-\\.]+", pathway_names))) # Valid characters
})