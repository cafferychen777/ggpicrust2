#!/usr/bin/env Rscript
# Comprehensive GO Pathway Support Validation for ggpicrust2
# Testing all aspects of GO functionality as per Linus principles:
# 1. Fail fast with clear reasons
# 2. No special cases - consistent behavior
# 3. Test real problems, not imaginary ones

# Load required packages
library(ggpicrust2)
library(testthat)
library(fgsea)

# Source the pathway_gsea.R to get functions
source("R/pathway_gsea.R")

cat("=== Comprehensive GO Pathway Support Validation ===\n")
cat("Starting validation at:", as.character(Sys.time()), "\n\n")

# Track validation results
validation_results <- list()

# 1. VALIDATE GO DATA STRUCTURE
cat("1. GO DATA STRUCTURE VALIDATION\n")
cat("================================\n")

test_that("GO mapping structure is correct", {
  # Load basic GO mapping
  go_mapping <- create_basic_go_mapping()
  
  # Basic structure tests
  expect_is(go_mapping, "data.frame")
  expect_true(all(c("go_id", "go_name", "category", "ko_members") %in% colnames(go_mapping)))
  
  # GO ID format validation - strict compliance
  invalid_go_ids <- go_mapping$go_id[!grepl("^GO:\\d{7}$", go_mapping$go_id)]
  if (length(invalid_go_ids) > 0) {
    cat("FAIL: Invalid GO ID formats found:", paste(invalid_go_ids, collapse = ", "), "\n")
    validation_results$go_id_format <- "FAIL"
  } else {
    cat("PASS: All GO IDs follow GO:XXXXXXX format\n")
    validation_results$go_id_format <- "PASS"
  }
  
  # Category validation
  valid_categories <- c("BP", "MF", "CC")
  invalid_categories <- go_mapping$category[!go_mapping$category %in% valid_categories]
  if (length(invalid_categories) > 0) {
    cat("FAIL: Invalid categories found:", paste(unique(invalid_categories), collapse = ", "), "\n")
    validation_results$go_categories <- "FAIL"
  } else {
    cat("PASS: All categories are BP, MF, or CC\n")
    validation_results$go_categories <- "PASS"
  }
  
  # KO format validation
  all_kos <- unique(unlist(strsplit(go_mapping$ko_members, ";")))
  invalid_kos <- all_kos[!grepl("^K\\d{5}$", all_kos)]
  if (length(invalid_kos) > 0) {
    cat("FAIL: Invalid KO formats found:", paste(head(invalid_kos, 5), collapse = ", "), "\n")
    validation_results$ko_format <- "FAIL"
  } else {
    cat("PASS: All KO IDs follow K##### format\n")
    validation_results$ko_format <- "PASS"
  }
})

# Count terms by category
go_mapping <- create_basic_go_mapping()
category_counts <- table(go_mapping$category)
cat("\nGO terms by category:\n")
for (cat_name in names(category_counts)) {
  cat(sprintf("%s: %d terms\n", cat_name, category_counts[cat_name]))
}
validation_results$category_counts <- category_counts

# 2. CATEGORY-SPECIFIC TESTING
cat("\n2. CATEGORY-SPECIFIC TESTING\n")
cat("==============================\n")

# Test each category
for (category in c("BP", "MF", "CC")) {
  cat(sprintf("\nTesting %s category:\n", category))
  
  # Test gene set preparation
  tryCatch({
    gene_sets <- prepare_gene_sets("GO", go_category = category)
    
    if (length(gene_sets) == 0) {
      cat(sprintf("FAIL: No gene sets returned for %s\n", category))
      validation_results[[paste0(category, "_gene_sets")]] <- "FAIL"
    } else {
      cat(sprintf("PASS: %d gene sets loaded for %s\n", length(gene_sets), category))
      validation_results[[paste0(category, "_gene_sets")]] <- "PASS"
      
      # Validate gene set sizes
      sizes <- sapply(gene_sets, length)
      cat(sprintf("  Gene set sizes: min=%d, max=%d, mean=%.1f\n", 
                  min(sizes), max(sizes), mean(sizes)))
      
      # Check for reasonable sizes (3-50 KOs per term)
      oversized <- names(gene_sets)[sizes > 50]
      undersized <- names(gene_sets)[sizes < 3]
      
      if (length(oversized) > 0) {
        cat(sprintf("  WARNING: %d terms with >50 KOs: %s\n", 
                    length(oversized), paste(head(oversized, 3), collapse = ", ")))
      }
      
      if (length(undersized) > 0) {
        cat(sprintf("  WARNING: %d terms with <3 KOs: %s\n", 
                    length(undersized), paste(head(undersized, 3), collapse = ", ")))
      }
    }
  }, error = function(e) {
    cat(sprintf("FAIL: Error loading %s gene sets: %s\n", category, e$message))
    validation_results[[paste0(category, "_gene_sets")]] <- "FAIL"
  })
}

# 3. GSEA WORKFLOW INTEGRATION
cat("\n3. GSEA WORKFLOW INTEGRATION\n")
cat("=============================\n")

# Load test data
data("ko_abundance", package = "ggpicrust2")
data("metadata", package = "ggpicrust2")

# Prepare abundance data
abundance_data <- as.data.frame(ko_abundance)
rownames(abundance_data) <- abundance_data[, "#NAME"]
abundance_data <- abundance_data[, -1]

# Use subset for testing (faster)
abundance_subset <- abundance_data[1:100, ]

cat("Test data prepared:\n")
cat(sprintf("  Features: %d\n", nrow(abundance_subset)))
cat(sprintf("  Samples: %d\n", ncol(abundance_subset)))
cat(sprintf("  Groups in metadata: %s\n", paste(unique(metadata$Environment), collapse = ", ")))

# Test GSEA with each category and ranking method
ranking_methods <- c("signal2noise", "t_test", "log2_ratio", "diff_abundance")
go_categories <- c("BP", "MF", "CC")

for (category in go_categories) {
  cat(sprintf("\nTesting GSEA with GO %s category:\n", category))
  
  for (rank_method in ranking_methods) {
    cat(sprintf("  Ranking method: %s... ", rank_method))
    
    tryCatch({
      gsea_results <- pathway_gsea(
        abundance = abundance_subset,
        metadata = metadata,
        group = "Environment",
        pathway_type = "GO",
        go_category = category,
        method = "fgsea",
        rank_method = rank_method,
        nperm = 100,  # Reduced for testing
        min_size = 2,
        max_size = 500,
        seed = 42
      )
      
      if (nrow(gsea_results) > 0) {
        cat(sprintf("PASS (%d results)\n", nrow(gsea_results)))
        validation_results[[paste0("gsea_", category, "_", rank_method)]] <- "PASS"
        
        # Validate result structure
        required_cols <- c("pathway_id", "NES", "pvalue", "p.adjust")
        missing_cols <- required_cols[!required_cols %in% colnames(gsea_results)]
        if (length(missing_cols) > 0) {
          cat(sprintf("    WARNING: Missing columns: %s\n", paste(missing_cols, collapse = ", ")))
        }
        
        # Check statistical validity
        if (any(gsea_results$pvalue < 0 | gsea_results$pvalue > 1)) {
          cat("    WARNING: Invalid p-values detected\n")
        }
        
        if (any(is.na(gsea_results$NES))) {
          cat("    WARNING: NA values in NES column\n")
        }
        
      } else {
        cat("PASS (no enriched pathways)\n")
        validation_results[[paste0("gsea_", category, "_", rank_method)]] <- "PASS"
      }
      
    }, error = function(e) {
      cat(sprintf("FAIL (%s)\n", e$message))
      validation_results[[paste0("gsea_", category, "_", rank_method)]] <- "FAIL"
    })
  }
}

# 4. GO ANNOTATION SYSTEM TEST
cat("\n4. GO ANNOTATION SYSTEM TEST\n")
cat("=============================\n")

# Test annotation for each category
for (category in go_categories) {
  cat(sprintf("\nTesting annotation for %s category:\n", category))
  
  tryCatch({
    # Run GSEA first
    gsea_results <- pathway_gsea(
      abundance = abundance_subset,
      metadata = metadata,
      group = "Environment",
      pathway_type = "GO",
      go_category = category,
      method = "fgsea",
      nperm = 50,
      min_size = 2,
      max_size = 500,
      seed = 42
    )
    
    if (nrow(gsea_results) > 0) {
      # Test annotation
      annotated_results <- gsea_pathway_annotation(
        gsea_results = gsea_results,
        pathway_type = "GO"
      )
      
      cat(sprintf("  PASS: Annotated %d results\n", nrow(annotated_results)))
      validation_results[[paste0("annotation_", category)]] <- "PASS"
      
      # Check annotation quality
      if ("pathway_name" %in% colnames(annotated_results)) {
        # Count how many have proper names vs just GO IDs
        proper_names <- sum(!is.na(annotated_results$pathway_name) & 
                           annotated_results$pathway_name != annotated_results$pathway_id)
        cat(sprintf("  Proper pathway names: %d/%d\n", proper_names, nrow(annotated_results)))
        
        if ("category" %in% colnames(annotated_results)) {
          cat(sprintf("  Category information included: YES\n"))
        } else {
          cat(sprintf("  Category information included: NO\n"))
        }
      }
      
    } else {
      cat("  No results to annotate (not an error)\n")
      validation_results[[paste0("annotation_", category)]] <- "PASS"
    }
    
  }, error = function(e) {
    cat(sprintf("  FAIL: %s\n", e$message))
    validation_results[[paste0("annotation_", category)]] <- "FAIL"
  })
}

# 5. KO MAPPING ACCURACY
cat("\n5. KO MAPPING ACCURACY\n")
cat("=======================\n")

go_mapping <- create_basic_go_mapping()

# Test biological relevance of mappings
cat("Testing KO mapping accuracy:\n")

# Check for realistic gene set sizes
gene_set_sizes <- sapply(strsplit(go_mapping$ko_members, ";"), length)
cat(sprintf("  Gene set size distribution:\n"))
cat(sprintf("    Min: %d KOs\n", min(gene_set_sizes)))
cat(sprintf("    Max: %d KOs\n", max(gene_set_sizes)))
cat(sprintf("    Mean: %.1f KOs\n", mean(gene_set_sizes)))
cat(sprintf("    Median: %.1f KOs\n", median(gene_set_sizes)))

# Check for overlapping KOs between terms (expected for biological data)
all_ko_assignments <- data.frame(
  go_id = rep(go_mapping$go_id, gene_set_sizes),
  ko_id = unlist(strsplit(go_mapping$ko_members, ";"))
)

ko_counts <- table(all_ko_assignments$ko_id)
multi_assigned_kos <- sum(ko_counts > 1)
total_unique_kos <- length(ko_counts)

cat(sprintf("  KO assignment overlap:\n"))
cat(sprintf("    Total unique KOs: %d\n", total_unique_kos))
cat(sprintf("    Multi-assigned KOs: %d (%.1f%%)\n", 
            multi_assigned_kos, 100 * multi_assigned_kos / total_unique_kos))

if (multi_assigned_kos / total_unique_kos > 0.7) {
  cat("    WARNING: Very high overlap (>70%) - may indicate mapping issues\n")
  validation_results$ko_mapping_accuracy <- "WARNING"
} else if (multi_assigned_kos / total_unique_kos < 0.1) {
  cat("    WARNING: Very low overlap (<10%) - may be too restrictive\n")
  validation_results$ko_mapping_accuracy <- "WARNING"
} else {
  cat("    PASS: Reasonable overlap for biological data\n")
  validation_results$ko_mapping_accuracy <- "PASS"
}

# 6. CROSS-CATEGORY ANALYSIS
cat("\n6. CROSS-CATEGORY ANALYSIS\n")
cat("===========================\n")

# Test running all categories together
cat("Testing all categories simultaneously:\n")

tryCatch({
  gsea_all <- pathway_gsea(
    abundance = abundance_subset,
    metadata = metadata,
    group = "Environment",
    pathway_type = "GO",
    go_category = "all",
    method = "fgsea",
    nperm = 100,
    min_size = 2,
    max_size = 500,
    seed = 42
  )
  
  cat(sprintf("  PASS: %d total results from all categories\n", nrow(gsea_all)))
  validation_results$cross_category <- "PASS"
  
  # Check category distribution if annotation is available
  if (nrow(gsea_all) > 0) {
    annotated_all <- gsea_pathway_annotation(gsea_all, pathway_type = "GO")
    
    if ("category" %in% colnames(annotated_all)) {
      category_dist <- table(annotated_all$category)
      cat("  Results by category:\n")
      for (cat_name in names(category_dist)) {
        cat(sprintf("    %s: %d\n", cat_name, category_dist[cat_name]))
      }
    }
  }
  
}, error = function(e) {
  cat(sprintf("  FAIL: %s\n", e$message))
  validation_results$cross_category <- "FAIL"
})

# 7. ERROR HANDLING
cat("\n7. ERROR HANDLING\n")
cat("==================\n")

cat("Testing error handling:\n")

# Test invalid GO category
cat("  Invalid go_category parameter... ")
error_caught <- tryCatch({
  pathway_gsea(
    abundance = abundance_subset[1:10, ],
    metadata = metadata,
    group = "Environment",
    pathway_type = "GO",
    go_category = "INVALID",
    method = "fgsea",
    nperm = 10
  )
  FALSE  # Should not reach here
}, error = function(e) {
  TRUE   # Error was caught
})

if (error_caught) {
  cat("PASS (error caught)\n")
  validation_results$error_invalid_category <- "PASS"
} else {
  cat("FAIL (error not caught)\n")
  validation_results$error_invalid_category <- "FAIL"
}

# Test with insufficient data
cat("  Insufficient sample size... ")
error_caught <- tryCatch({
  pathway_gsea(
    abundance = abundance_subset[1:10, 1:2],  # Only 2 samples
    metadata = metadata[1:2, ],
    group = "Environment",
    pathway_type = "GO",
    go_category = "BP",
    method = "fgsea",
    nperm = 10
  )
  FALSE
}, error = function(e) {
  TRUE
})

if (error_caught) {
  cat("PASS (error caught)\n")
  validation_results$error_insufficient_samples <- "PASS"
} else {
  cat("FAIL (error not caught)\n")
  validation_results$error_insufficient_samples <- "FAIL"
}

# 8. PERFORMANCE BENCHMARKING
cat("\n8. PERFORMANCE BENCHMARKING\n")
cat("============================\n")

cat("Testing performance with realistic data sizes:\n")

# Test with larger dataset
large_subset <- abundance_data[1:500, ]  # 500 features

performance_times <- list()

for (category in c("BP", "MF", "CC")) {
  cat(sprintf("  %s category... ", category))
  
  start_time <- Sys.time()
  
  tryCatch({
    gsea_results <- pathway_gsea(
      abundance = large_subset,
      metadata = metadata,
      group = "Environment",
      pathway_type = "GO",
      go_category = category,
      method = "fgsea",
      nperm = 100,
      min_size = 3,
      max_size = 50,
      seed = 42
    )
    
    end_time <- Sys.time()
    elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
    performance_times[[category]] <- elapsed
    
    cat(sprintf("%.2f seconds (%d results)\n", elapsed, nrow(gsea_results)))
    
    # Performance check
    if (elapsed > 60) {
      cat(sprintf("    WARNING: Slow performance (>60s)\n"))
    }
    
  }, error = function(e) {
    cat(sprintf("FAIL: %s\n", e$message))
    performance_times[[category]] <- NA
  })
}

validation_results$performance_times <- performance_times

# 9. REPRODUCIBILITY TEST
cat("\n9. REPRODUCIBILITY TEST\n")
cat("========================\n")

cat("Testing result reproducibility:\n")

# Run same analysis twice with same seed
results1 <- pathway_gsea(
  abundance = abundance_subset,
  metadata = metadata,
  group = "Environment",
  pathway_type = "GO",
  go_category = "BP",
  method = "fgsea",
  nperm = 100,
  seed = 42
)

results2 <- pathway_gsea(
  abundance = abundance_subset,
  metadata = metadata,
  group = "Environment",
  pathway_type = "GO",
  go_category = "BP",
  method = "fgsea",
  nperm = 100,
  seed = 42
)

if (nrow(results1) > 0 && nrow(results2) > 0) {
  # Check if results are identical
  if (nrow(results1) == nrow(results2)) {
    # Sort by pathway_id for comparison
    results1_sorted <- results1[order(results1$pathway_id), ]
    results2_sorted <- results2[order(results2$pathway_id), ]
    
    # Compare key columns
    nes_identical <- all.equal(results1_sorted$NES, results2_sorted$NES, tolerance = 1e-10)
    pval_identical <- all.equal(results1_sorted$pvalue, results2_sorted$pvalue, tolerance = 1e-10)
    
    if (isTRUE(nes_identical) && isTRUE(pval_identical)) {
      cat("  PASS: Results are reproducible\n")
      validation_results$reproducibility <- "PASS"
    } else {
      cat("  FAIL: Results differ between runs\n")
      validation_results$reproducibility <- "FAIL"
    }
  } else {
    cat("  FAIL: Different number of results\n")
    validation_results$reproducibility <- "FAIL"
  }
} else {
  cat("  PASS: No results to compare (consistent)\n")
  validation_results$reproducibility <- "PASS"
}

# FINAL VALIDATION SUMMARY
cat("\n" + rep("=", 50) + "\n")
cat("FINAL VALIDATION SUMMARY\n")
cat(rep("=", 50) + "\n")

pass_count <- sum(sapply(validation_results, function(x) {
  if (is.character(x)) x == "PASS" else FALSE
}))

fail_count <- sum(sapply(validation_results, function(x) {
  if (is.character(x)) x == "FAIL" else FALSE
}))

warning_count <- sum(sapply(validation_results, function(x) {
  if (is.character(x)) x == "WARNING" else FALSE
}))

total_tests <- pass_count + fail_count + warning_count

cat(sprintf("Total tests: %d\n", total_tests))
cat(sprintf("PASSED: %d (%.1f%%)\n", pass_count, 100 * pass_count / total_tests))
cat(sprintf("FAILED: %d (%.1f%%)\n", fail_count, 100 * fail_count / total_tests))
cat(sprintf("WARNINGS: %d (%.1f%%)\n", warning_count, 100 * warning_count / total_tests))

# Detailed results
cat("\nDetailed test results:\n")
for (test_name in names(validation_results)) {
  result <- validation_results[[test_name]]
  if (is.character(result)) {
    status_symbol <- switch(result,
                           "PASS" = "✓",
                           "FAIL" = "✗",
                           "WARNING" = "⚠",
                           "?")
    cat(sprintf("  %s %s: %s\n", status_symbol, test_name, result))
  }
}

# Critical issues check
critical_failures <- c(
  "go_id_format", "go_categories", "ko_format",
  "cross_category", "error_invalid_category"
)

critical_failed <- sapply(critical_failures, function(test) {
  if (test %in% names(validation_results)) {
    validation_results[[test]] == "FAIL"
  } else {
    FALSE
  }
})

if (any(critical_failed)) {
  cat("\n⚠ CRITICAL ISSUES DETECTED ⚠\n")
  cat("The following critical tests failed:\n")
  for (test in critical_failures[critical_failed]) {
    cat(sprintf("  - %s\n", test))
  }
  cat("\nGO functionality is NOT ready for release.\n")
} else {
  cat("\n✓ VALIDATION COMPLETE ✓\n")
  if (fail_count == 0) {
    cat("All tests passed. GO functionality is ready for production use.\n")
  } else {
    cat(sprintf("Minor issues detected (%d failures), but core functionality is working.\n", fail_count))
    cat("Review failed tests before release.\n")
  }
}

cat("\nValidation completed at:", as.character(Sys.time()), "\n")

# Save validation results
saveRDS(validation_results, "go_validation_results.rds")
cat("Validation results saved to: go_validation_results.rds\n")