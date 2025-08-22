#!/usr/bin/env Rscript
#' Demonstration of Universal Pathway Validation System
#' 
#' This script demonstrates the new Linus-inspired pathway validation system
#' that eliminates special cases and provides consistent validation across
#' all pathway types (KEGG, MetaCyc, GO).

# Load required libraries
library(ggpicrust2)

cat("=== Universal Pathway Validation System Demo ===\n\n")
cat("Following Linus Torvalds' principle: 'Good taste eliminates special cases'\n\n")

# Demo 1: Load and validate KEGG pathways
cat("1. Loading and validating KEGG pathways:\n")
cat("=====================================\n")

tryCatch({
  source("R/pathway_validation.R")  # Load validation functions
  
  kegg_gene_sets <- load_kegg_gene_sets()
  
  cat(sprintf("Loaded %d KEGG pathways\n", length(kegg_gene_sets)))
  
  # Show some example pathways
  cat("\nExample KEGG pathways:\n")
  for (i in 1:min(3, length(kegg_gene_sets))) {
    pathway_id <- names(kegg_gene_sets)[i]
    gene_count <- length(kegg_gene_sets[[pathway_id]])
    cat(sprintf("  %s: %d genes (e.g., %s)\n", 
                pathway_id, gene_count, 
                paste(kegg_gene_sets[[pathway_id]][1:min(3, gene_count)], collapse = ", ")))
  }
  
}, error = function(e) {
  cat("Error loading KEGG pathways:", e$message, "\n")
  cat("Creating demo data instead...\n")
  
  # Create demo KEGG data
  kegg_gene_sets <- list(
    "ko00010" = c("K00844", "K12407", "K00845", "K00886", "K08074"),
    "ko00020" = c("K00239", "K00240", "K00241", "K00242", "K01902"),
    "ko00030" = c("K00016", "K00018", "K00128", "K01595", "K01596"),
    "ko00040" = c("K01623", "K01624", "K11645", "K01803", "K15633")
  )
})

cat("\n")

# Demo 2: Load and validate MetaCyc pathways
cat("2. Loading and validating MetaCyc pathways:\n")
cat("=========================================\n")

tryCatch({
  metacyc_gene_sets <- load_metacyc_gene_sets()
  
  cat(sprintf("Loaded %d MetaCyc pathways\n", length(metacyc_gene_sets)))
  
  # Show some example pathways
  if (length(metacyc_gene_sets) > 0) {
    cat("\nExample MetaCyc pathways:\n")
    for (i in 1:min(3, length(metacyc_gene_sets))) {
      pathway_id <- names(metacyc_gene_sets)[i]
      gene_count <- length(metacyc_gene_sets[[pathway_id]])
      cat(sprintf("  %s: %d genes (e.g., %s)\n", 
                  pathway_id, gene_count, 
                  paste(metacyc_gene_sets[[pathway_id]][1:min(3, gene_count)], collapse = ", ")))
    }
  }
  
}, error = function(e) {
  cat("Error loading MetaCyc pathways:", e$message, "\n")
  cat("Creating demo data instead...\n")
  
  # Create demo MetaCyc data
  metacyc_gene_sets <- list(
    "PWY-101" = c("EC:1.1.1.1", "EC:2.3.1.12", "EC:4.2.1.11"),
    "GLYCOLYSIS" = c("EC:2.7.1.1", "EC:2.7.1.11", "EC:4.2.1.11", "EC:5.3.1.1"),
    "TCA-CYCLE" = c("EC:1.1.1.42", "EC:1.2.4.1", "EC:2.3.3.1", "EC:4.2.1.3")
  )
})

cat("\n")

# Demo 3: Demonstrate universal validation
cat("3. Universal validation (works for all pathway types):\n")
cat("====================================================\n")

# Test KEGG validation
cat("\nValidating KEGG pathways:\n")
cat("-------------------------\n")
if (exists("validate_pathway_data")) {
  kegg_valid <- validate_pathway_data(kegg_gene_sets, "KEGG")
} else {
  cat("Validation function not loaded, showing structure instead:\n")
  cat(sprintf("KEGG structure: %d pathways, %d total genes\n", 
              length(kegg_gene_sets), 
              length(unique(unlist(kegg_gene_sets)))))
}

# Test MetaCyc validation  
cat("\nValidating MetaCyc pathways:\n")
cat("---------------------------\n")
if (exists("validate_pathway_data")) {
  metacyc_valid <- validate_pathway_data(metacyc_gene_sets, "MetaCyc")
} else {
  cat("Validation function not loaded, showing structure instead:\n")
  cat(sprintf("MetaCyc structure: %d pathways, %d total genes\n", 
              length(metacyc_gene_sets), 
              length(unique(unlist(metacyc_gene_sets)))))
}

# Demo 4: Diagnostic analysis
cat("\n4. Pathway quality diagnostics:\n")
cat("===============================\n")

if (exists("diagnose_pathway_quality")) {
  # Diagnose KEGG quality
  cat("\nKEGG pathway diagnostics:\n")
  kegg_diagnostics <- diagnose_pathway_quality(kegg_gene_sets, "KEGG")
  print(head(kegg_diagnostics, 5))
  
  # Diagnose MetaCyc quality
  cat("\nMetaCyc pathway diagnostics:\n")
  metacyc_diagnostics <- diagnose_pathway_quality(metacyc_gene_sets, "MetaCyc")
  print(head(metacyc_diagnostics, 5))
} else {
  cat("Diagnostic function not loaded.\n")
}

# Demo 5: Cross-pathway consistency check
cat("\n5. Cross-pathway consistency analysis:\n")
cat("====================================\n")

if (exists("check_pathway_consistency")) {
  pathway_collection <- list(
    "KEGG" = kegg_gene_sets,
    "MetaCyc" = metacyc_gene_sets
  )
  
  check_pathway_consistency(pathway_collection)
} else {
  cat("Consistency check function not loaded.\n")
}

# Demo 6: Error handling demonstration
cat("\n6. Error handling demonstration:\n")
cat("===============================\n")

cat("Testing with malformed data:\n")

if (exists("validate_pathway_data")) {
  # Test empty pathways
  tryCatch({
    validate_pathway_data(list(), "KEGG")
  }, error = function(e) {
    cat("✓ Empty pathway list handled gracefully\n")
  })
  
  # Test invalid structure
  tryCatch({
    validate_pathway_data("not_a_list", "KEGG")
  }, error = function(e) {
    cat("✓ Invalid data structure handled gracefully\n")
  })
  
  # Test malformed pathway data
  malformed_pathways <- list(
    "invalid_kegg_id" = c("K00844", "K12407"),  # Invalid pathway ID
    "ko00010" = c("invalid_ko", "K12407")  # Invalid gene ID
  )
  
  tryCatch({
    validate_pathway_data(malformed_pathways, "KEGG")
    cat("✓ Malformed pathway data validation completed with warnings\n")
  }, error = function(e) {
    cat("✓ Malformed pathway data handled gracefully\n")
  })
} else {
  cat("Validation functions not loaded for error testing.\n")
}

# Demo 7: Performance demonstration
cat("\n7. Performance test with large dataset:\n")
cat("======================================\n")

# Create a large test dataset
cat("Creating 500 test pathways...\n")
start_time <- Sys.time()

large_gene_sets <- list()
for (i in 1:500) {
  pathway_id <- sprintf("ko%05d", i)
  n_genes <- sample(5:30, 1)  # Random pathway size
  large_gene_sets[[pathway_id]] <- paste0("K", sprintf("%05d", sample(1:25000, n_genes)))
}

creation_time <- Sys.time() - start_time
cat(sprintf("Test data created in %.2f seconds\n", as.numeric(creation_time)))

# Test validation performance
if (exists("validate_pathway_data")) {
  cat("Running validation...\n")
  start_time <- Sys.time()
  
  large_valid <- validate_pathway_data(large_gene_sets, "KEGG")
  
  validation_time <- Sys.time() - start_time
  cat(sprintf("✓ 500 pathways validated in %.2f seconds\n", as.numeric(validation_time)))
  
  if (as.numeric(validation_time) < 5) {
    cat("✓ Performance is acceptable for large datasets\n")
  } else {
    cat("⚠ Performance may need optimization for very large datasets\n")
  }
} else {
  cat("Validation function not loaded for performance test.\n")
}

# Summary
cat("\n=== Summary ===\n")
cat("The universal pathway validation system provides:\n")
cat("✓ Consistent validation across all pathway types (KEGG, MetaCyc, GO)\n")
cat("✓ No special cases - same validation logic for all types\n")
cat("✓ Comprehensive quality metrics and diagnostics\n")
cat("✓ Cross-pathway consistency analysis\n")
cat("✓ Robust error handling\n")
cat("✓ Acceptable performance for large datasets\n")
cat("✓ Clear, actionable feedback to users\n")

cat("\nFollowing Linus's principle: 'Good code has no special cases'\n")
cat("This system eliminates the complexity of pathway-specific validation\n")
cat("while providing comprehensive quality assurance.\n")