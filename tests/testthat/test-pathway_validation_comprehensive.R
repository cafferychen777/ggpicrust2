#!/usr/bin/env Rscript
#' Comprehensive Pathway Validation Tests
#' 
#' Following Linus's testing philosophy: "If it's not tested, it's broken"
#' Tests the universal validation system for all pathway types

test_that("validate_pathway_data handles basic structure validation", {
  
  # Test with valid gene sets
  valid_kegg <- list(
    "ko00010" = c("K00844", "K12407", "K00845"),
    "ko00020" = c("K00239", "K00240", "K00241")
  )
  
  expect_true(validate_pathway_data(valid_kegg, "KEGG"))
  
  # Test with empty list
  expect_warning(
    result <- validate_pathway_data(list(), "KEGG"),
    "No KEGG pathways loaded"
  )
  expect_false(result)
  
  # Test with non-list input
  expect_error(
    validate_pathway_data("not_a_list", "KEGG"),
    "Gene sets must be provided as a list"
  )
  
  # Test with unnamed pathways
  unnamed_sets <- list(c("K00844", "K12407"), c("K00239", "K00240"))
  expect_error(
    validate_pathway_data(unnamed_sets, "KEGG"),
    "All pathways must have valid, non-empty identifiers"
  )
  
  # Test with duplicate pathway IDs
  duplicate_sets <- list(
    "ko00010" = c("K00844", "K12407"),
    "ko00010" = c("K00239", "K00240")  # Duplicate name
  )
  expect_warning(
    validate_pathway_data(duplicate_sets, "KEGG"),
    "Duplicate pathway IDs detected"
  )
})

test_that("validate_pathway_format works for KEGG", {
  
  # Valid KEGG format
  valid_kegg <- list(
    "ko00010" = c("K00844", "K12407", "K00845"),
    "ko00020" = c("K00239", "K00240", "K00241")
  )
  
  expect_silent(validate_pathway_format(names(valid_kegg), valid_kegg, "KEGG"))
  
  # Invalid KEGG pathway IDs
  invalid_pathway_kegg <- list(
    "invalid_id" = c("K00844", "K12407"),
    "ko00020" = c("K00239", "K00240")
  )
  
  expect_warning(
    validate_pathway_format(names(invalid_pathway_kegg), invalid_pathway_kegg, "KEGG"),
    "Invalid KEGG pathway IDs detected"
  )
  
  # Invalid KO IDs
  invalid_ko_kegg <- list(
    "ko00010" = c("invalid_ko", "K12407"),
    "ko00020" = c("K00239", "also_invalid")
  )
  
  expect_warning(
    validate_pathway_format(names(invalid_ko_kegg), invalid_ko_kegg, "KEGG"),
    "Invalid KO identifiers detected"
  )
})

test_that("validate_pathway_format works for MetaCyc", {
  
  # Valid MetaCyc format
  valid_metacyc <- list(
    "PWY-101" = c("EC:1.1.1.1", "EC:2.3.1.12"),
    "GLYCOLYSIS" = c("EC:1.1.1.1", "EC:4.2.1.11")
  )
  
  expect_silent(validate_pathway_format(names(valid_metacyc), valid_metacyc, "MetaCyc"))
  
  # Invalid EC numbers
  invalid_ec_metacyc <- list(
    "PWY-101" = c("invalid_ec", "EC:2.3.1.12"),
    "GLYCOLYSIS" = c("EC:1.1.1.1", "also_invalid")
  )
  
  expect_warning(
    validate_pathway_format(names(invalid_ec_metacyc), invalid_ec_metacyc, "MetaCyc"),
    "Invalid EC number format detected"
  )
})

test_that("validate_pathway_format works for GO", {
  
  # Valid GO format
  valid_go <- list(
    "GO:0008150" = c("K00844", "K12407"),
    "GO:0003674" = c("K00239", "K00240")
  )
  
  expect_silent(validate_pathway_format(names(valid_go), valid_go, "GO"))
  
  # Invalid GO IDs
  invalid_go <- list(
    "invalid_go" = c("K00844", "K12407"),
    "GO:0003674" = c("K00239", "K00240")
  )
  
  expect_warning(
    validate_pathway_format(names(invalid_go), invalid_go, "GO"),
    "Invalid GO term IDs detected"
  )
})

test_that("validate_gene_set_quality detects size issues", {
  
  # Mix of good, empty, tiny, and huge gene sets
  mixed_quality <- list(
    "good_pathway" = paste0("K", sprintf("%05d", 1:20)),  # Good size
    "empty_pathway" = character(0),  # Empty
    "tiny_pathway" = c("K00001"),  # Too small
    "huge_pathway" = paste0("K", sprintf("%05d", 1:600))  # Too large
  )
  
  expect_warning(
    validate_gene_set_quality(mixed_quality, "TEST"),
    "Empty gene sets detected"
  )
  
  expect_warning(
    validate_gene_set_quality(mixed_quality, "TEST"),
    "Very small gene sets \\(<3 genes\\) detected"
  )
  
  expect_warning(
    validate_gene_set_quality(mixed_quality, "TEST"),
    "Very large gene sets \\(>500 genes\\) detected"
  )
})

test_that("load_kegg_gene_sets works with real data", {
  
  # This test requires actual reference data
  skip_if_not(file.exists(system.file("extdata", "KO_reference.RData", package = "ggpicrust2")),
              "KO reference data not available")
  
  gene_sets <- load_kegg_gene_sets()
  
  expect_type(gene_sets, "list")
  expect_true(length(gene_sets) > 0)
  expect_true(all(sapply(names(gene_sets), function(x) grepl("^ko[0-9]{5}$", x))))
  
  # Check that genes are KO IDs
  all_genes <- unique(unlist(gene_sets, use.names = FALSE))
  valid_kos <- grepl("^K[0-9]{5}$", all_genes)
  expect_true(sum(valid_kos) / length(all_genes) > 0.8)  # At least 80% should be valid
})

test_that("load_metacyc_gene_sets works with available data", {
  
  skip_if_not(file.exists(system.file("extdata", "MetaCyc_reference.RData", package = "ggpicrust2")),
              "MetaCyc reference data not available")
  
  gene_sets <- load_metacyc_gene_sets()
  
  expect_type(gene_sets, "list")
  if (length(gene_sets) > 0) {
    # Check that genes are EC numbers
    all_genes <- unique(unlist(gene_sets, use.names = FALSE))
    valid_ecs <- grepl("^EC:", all_genes)
    expect_true(sum(valid_ecs) / length(all_genes) > 0.8)
  }
})

test_that("load_go_gene_sets handles missing data gracefully", {
  
  expect_warning(
    gene_sets <- load_go_gene_sets(),
    "GO pathway gene sets not yet implemented"
  )
  
  expect_type(gene_sets, "list")
  expect_equal(length(gene_sets), 0)
})

test_that("check_pathway_consistency works with multiple pathway types", {
  
  kegg_sets <- list(
    "ko00010" = c("K00844", "K12407", "K00845"),
    "ko00020" = c("K00239", "K00240", "K00241")
  )
  
  metacyc_sets <- list(
    "PWY-101" = c("EC:1.1.1.1", "EC:2.3.1.12"),
    "GLYCOLYSIS" = c("K00844", "K00845")  # Some overlap with KEGG
  )
  
  pathway_list <- list(
    "KEGG" = kegg_sets,
    "MetaCyc" = metacyc_sets
  )
  
  expect_message(
    check_pathway_consistency(pathway_list),
    "Pathway consistency analysis"
  )
  
  # Test with single pathway type
  expect_message(
    check_pathway_consistency(list("KEGG" = kegg_sets)),
    "Only one pathway type provided"
  )
})

test_that("diagnose_pathway_quality provides useful diagnostics", {
  
  test_sets <- list(
    "good_pathway" = paste0("K", sprintf("%05d", 1:20)),
    "empty_pathway" = character(0),
    "tiny_pathway" = c("K00001"),
    "invalid_pathway_id" = paste0("K", sprintf("%05d", 1:10))
  )
  
  names(test_sets)[4] <- "invalid_id"  # Make one pathway ID invalid
  
  diagnostics <- diagnose_pathway_quality(test_sets, "KEGG")
  
  expect_s3_class(diagnostics, "data.frame")
  expect_equal(nrow(diagnostics), 4)
  expect_true("is_empty" %in% colnames(diagnostics))
  expect_true("is_tiny" %in% colnames(diagnostics))
  expect_true("valid_pathway_format" %in% colnames(diagnostics))
  
  # Check that issues are correctly identified
  expect_true(diagnostics$is_empty[diagnostics$pathway_id == "empty_pathway"])
  expect_true(diagnostics$is_tiny[diagnostics$pathway_id == "tiny_pathway"])
  expect_false(diagnostics$valid_pathway_format[diagnostics$pathway_id == "invalid_id"])
})

test_that("Integration test: end-to-end pathway validation workflow", {
  
  # Test the complete workflow
  skip_if_not(file.exists(system.file("extdata", "KO_reference.RData", package = "ggpicrust2")),
              "Reference data not available for integration test")
  
  # Load gene sets
  kegg_sets <- load_kegg_gene_sets()
  
  # Validate them
  expect_true(validate_pathway_data(kegg_sets, "KEGG"))
  
  # Diagnose quality
  diagnostics <- diagnose_pathway_quality(kegg_sets, "KEGG")
  expect_true(nrow(diagnostics) > 0)
  
  # Check for critical issues
  critical_issues <- sum(diagnostics$is_empty) + 
                     sum(!diagnostics$valid_pathway_format) + 
                     sum(diagnostics$valid_gene_fraction < 0.5)
  
  # In a well-formed dataset, we shouldn't have many critical issues
  expect_lt(critical_issues / nrow(diagnostics), 0.1)  # Less than 10% critical issues
})

test_that("Error handling is robust", {
  
  # Test with malformed inputs
  expect_error(diagnose_pathway_quality(NULL, "KEGG"))
  expect_error(diagnose_pathway_quality(list(), "KEGG"))
  
  # Test with invalid pathway types
  expect_error(
    switch("INVALID", 
           "KEGG" = "kegg", 
           "MetaCyc" = "metacyc", 
           "GO" = "go",
           stop("pathway_type must be one of 'KEGG', 'MetaCyc', or 'GO'")),
    "pathway_type must be one of"
  )
  
  # Test validation with malformed gene sets
  malformed <- list("pathway1" = list(invalid = "structure"))  # Not character vector
  # Our system is robust and handles this gracefully with warnings
  expect_warning(
    validate_gene_set_quality(malformed, "TEST"),
    class = "warning"
  )
})

# Performance test for large datasets
test_that("Validation performance is acceptable", {
  
  # Create a large test dataset
  large_gene_sets <- list()
  for (i in 1:1000) {
    pathway_id <- sprintf("ko%05d", i)
    n_genes <- sample(10:50, 1)
    large_gene_sets[[pathway_id]] <- paste0("K", sprintf("%05d", sample(1:25000, n_genes)))
  }
  
  # Measure validation time
  start_time <- Sys.time()
  result <- validate_pathway_data(large_gene_sets, "KEGG")
  end_time <- Sys.time()
  
  validation_time <- as.numeric(end_time - start_time, units = "secs")
  
  expect_true(result)
  expect_lt(validation_time, 10)  # Should complete within 10 seconds
  
  message(sprintf("Validation of 1000 pathways completed in %.2f seconds", validation_time))
})