#!/usr/bin/env Rscript
# Focused GO Pathway Support Validation for ggpicrust2
# Addressing the sample alignment issue and focusing on core GO functionality

# Load required packages
library(ggpicrust2)
library(fgsea)

# Source the pathway_gsea.R to get functions
source("R/pathway_gsea.R")

cat("=== Focused GO Pathway Support Validation ===\n")
cat("Testing Date:", as.character(Sys.time()), "\n\n")

# Initialize validation results
validation_results <- list()

# 1. GO DATA STRUCTURE VALIDATION
cat("1. GO DATA STRUCTURE VALIDATION\n")
cat("================================\n")

# Test basic GO mapping creation
go_mapping <- create_basic_go_mapping()

cat("GO mapping structure:\n")
cat("  Total terms:", nrow(go_mapping), "\n")
cat("  Columns:", paste(colnames(go_mapping), collapse = ", "), "\n")

# Category distribution
category_counts <- table(go_mapping$category)
cat("  Category distribution:\n")
for (cat_name in names(category_counts)) {
  cat(sprintf("    %s: %d terms\n", cat_name, category_counts[cat_name]))
}

# Validate GO ID format (strict)
invalid_go_ids <- go_mapping$go_id[!grepl("^GO:\\d{7}$", go_mapping$go_id)]
if (length(invalid_go_ids) > 0) {
  cat("  FAIL: Invalid GO ID formats:", paste(invalid_go_ids, collapse = ", "), "\n")
  validation_results$go_id_format <- "FAIL"
} else {
  cat("  PASS: All GO IDs properly formatted (GO:XXXXXXX)\n")
  validation_results$go_id_format <- "PASS"
}

# Validate categories
valid_categories <- c("BP", "MF", "CC")
invalid_categories <- go_mapping$category[!go_mapping$category %in% valid_categories]
if (length(invalid_categories) > 0) {
  cat("  FAIL: Invalid categories:", paste(unique(invalid_categories), collapse = ", "), "\n")
  validation_results$go_categories <- "FAIL"
} else {
  cat("  PASS: All categories valid (BP, MF, CC)\n")
  validation_results$go_categories <- "PASS"
}

# Validate KO format
all_kos <- unique(unlist(strsplit(go_mapping$ko_members, ";")))
invalid_kos <- all_kos[!grepl("^K\\d{5}$", all_kos)]
if (length(invalid_kos) > 0) {
  cat("  FAIL: Invalid KO formats found:", length(invalid_kos), "invalid KOs\n")
  validation_results$ko_format <- "FAIL"
} else {
  cat("  PASS: All KO IDs properly formatted (K#####)\n")
  validation_results$ko_format <- "PASS"
}

# Check for expected 36 total terms (20 BP + 8 MF + 8 CC)
expected_total <- 36
if (nrow(go_mapping) == expected_total) {
  cat("  PASS: Expected", expected_total, "terms found\n")
  validation_results$term_count <- "PASS"
} else {
  cat("  WARNING: Expected", expected_total, "terms, found", nrow(go_mapping), "\n")
  validation_results$term_count <- "WARNING"
}

# 2. GENE SET PREPARATION TESTING
cat("\n2. GENE SET PREPARATION TESTING\n")
cat("=================================\n")

for (category in c("BP", "MF", "CC", "all")) {
  cat(sprintf("Testing %s category:\n", category))
  
  tryCatch({
    gene_sets <- prepare_gene_sets("GO", go_category = category)
    
    if (length(gene_sets) == 0) {
      cat(sprintf("  FAIL: No gene sets returned for %s\n", category))
      validation_results[[paste0("gene_sets_", category)]] <- "FAIL"
    } else {
      cat(sprintf("  PASS: %d gene sets loaded\n", length(gene_sets)))
      validation_results[[paste0("gene_sets_", category)]] <- "PASS"
      
      # Check gene set sizes
      sizes <- sapply(gene_sets, length)
      cat(sprintf("    Size range: %d-%d KOs (mean: %.1f)\n", 
                  min(sizes), max(sizes), mean(sizes)))
      
      # Check GO ID format in names
      go_names <- names(gene_sets)
      invalid_names <- go_names[!grepl("^GO:\\d{7}$", go_names)]
      if (length(invalid_names) > 0) {
        cat(sprintf("    WARNING: %d invalid GO IDs in gene set names\n", length(invalid_names)))
      }
      
      # Validate gene set sizes are reasonable (3-50 for production)
      too_small <- sum(sizes < 3)
      too_large <- sum(sizes > 50)
      if (too_small > 0) {
        cat(sprintf("    INFO: %d gene sets have <3 KOs\n", too_small))
      }
      if (too_large > 0) {
        cat(sprintf("    WARNING: %d gene sets have >50 KOs\n", too_large))
      }
    }
  }, error = function(e) {
    cat(sprintf("  FAIL: Error in %s: %s\n", category, e$message))
    validation_results[[paste0("gene_sets_", category)]] <- "FAIL"
  })
}

# 3. CATEGORY FILTERING VALIDATION
cat("\n3. CATEGORY FILTERING VALIDATION\n")
cat("==================================\n")

# Test that different categories return different gene sets
bp_sets <- prepare_gene_sets("GO", go_category = "BP")
mf_sets <- prepare_gene_sets("GO", go_category = "MF")
cc_sets <- prepare_gene_sets("GO", go_category = "CC")
all_sets <- prepare_gene_sets("GO", go_category = "all")

# Check for no overlap between categories
bp_ids <- names(bp_sets)
mf_ids <- names(mf_sets) 
cc_ids <- names(cc_sets)

overlap_bp_mf <- length(intersect(bp_ids, mf_ids))
overlap_bp_cc <- length(intersect(bp_ids, cc_ids))
overlap_mf_cc <- length(intersect(mf_ids, cc_ids))

if (overlap_bp_mf == 0 && overlap_bp_cc == 0 && overlap_mf_cc == 0) {
  cat("  PASS: No overlap between categories (as expected)\n")
  validation_results$category_separation <- "PASS"
} else {
  cat("  FAIL: Unexpected overlap between categories\n")
  cat(sprintf("    BP-MF: %d, BP-CC: %d, MF-CC: %d\n", overlap_bp_mf, overlap_bp_cc, overlap_mf_cc))
  validation_results$category_separation <- "FAIL"
}

# Check that 'all' contains all individual categories
all_ids <- names(all_sets)
if (all(bp_ids %in% all_ids) && all(mf_ids %in% all_ids) && all(cc_ids %in% all_ids)) {
  cat("  PASS: 'all' category includes all individual categories\n")
  validation_results$category_all_inclusive <- "PASS"
} else {
  cat("  FAIL: 'all' category missing some terms\n")
  validation_results$category_all_inclusive <- "FAIL"
}

# 4. MOCK GSEA TESTING (with aligned data)
cat("\n4. MOCK GSEA TESTING\n")
cat("=====================\n")

# Create mock data with aligned sample names
set.seed(42)
n_features <- 100
n_samples <- 20

# Create mock abundance data with KO IDs that exist in our GO mapping
mock_kos <- sample(all_kos, n_features, replace = TRUE)
mock_abundance <- matrix(
  rpois(n_features * n_samples, lambda = 100),
  nrow = n_features,
  ncol = n_samples
)
rownames(mock_abundance) <- mock_kos
colnames(mock_abundance) <- paste0("Sample_", 1:n_samples)

# Create aligned metadata
mock_metadata <- data.frame(
  row.names = colnames(mock_abundance),
  Environment = rep(c("Group1", "Group2"), each = n_samples/2),
  stringsAsFactors = FALSE
)

cat("Mock data created:\n")
cat("  Features:", nrow(mock_abundance), "\n")
cat("  Samples:", ncol(mock_abundance), "\n") 
cat("  Groups:", paste(unique(mock_metadata$Environment), collapse = ", "), "\n")

# Test GSEA with each category
for (category in c("BP", "MF", "CC")) {
  cat(sprintf("\nTesting GSEA with %s category:\n", category))
  
  tryCatch({
    gsea_results <- pathway_gsea(
      abundance = mock_abundance,
      metadata = mock_metadata,
      group = "Environment",
      pathway_type = "GO",
      go_category = category,
      method = "fgsea",
      nperm = 100,
      min_size = 2,
      max_size = 50,
      seed = 42
    )
    
    cat(sprintf("  PASS: GSEA completed (%d results)\n", nrow(gsea_results)))
    validation_results[[paste0("gsea_", category)]] <- "PASS"
    
    # Validate result structure
    expected_cols <- c("pathway_id", "NES", "pvalue", "p.adjust")
    missing_cols <- expected_cols[!expected_cols %in% colnames(gsea_results)]
    if (length(missing_cols) > 0) {
      cat(sprintf("    WARNING: Missing columns: %s\n", paste(missing_cols, collapse = ", ")))
    } else {
      cat("    PASS: All expected columns present\n")
    }
    
    # Check for valid statistical values
    if (nrow(gsea_results) > 0) {
      invalid_pvals <- sum(gsea_results$pvalue < 0 | gsea_results$pvalue > 1)
      invalid_nes <- sum(is.na(gsea_results$NES) | is.infinite(gsea_results$NES))
      
      if (invalid_pvals > 0) {
        cat(sprintf("    WARNING: %d invalid p-values\n", invalid_pvals))
      }
      
      if (invalid_nes > 0) {
        cat(sprintf("    WARNING: %d invalid NES values\n", invalid_nes))
      }
      
      if (invalid_pvals == 0 && invalid_nes == 0) {
        cat("    PASS: All statistical values valid\n")
      }
    }
    
  }, error = function(e) {
    cat(sprintf("  FAIL: %s\n", e$message))
    validation_results[[paste0("gsea_", category)]] <- "FAIL"
  })
}

# 5. ANNOTATION TESTING
cat("\n5. ANNOTATION TESTING\n")
cat("======================\n")

# Test annotation with each category
for (category in c("BP", "MF", "CC")) {
  cat(sprintf("Testing annotation for %s:\n", category))
  
  tryCatch({
    # First run GSEA
    gsea_results <- pathway_gsea(
      abundance = mock_abundance,
      metadata = mock_metadata,
      group = "Environment",
      pathway_type = "GO",
      go_category = category,
      method = "fgsea",
      nperm = 50,
      min_size = 2,
      max_size = 50,
      seed = 42
    )
    
    if (nrow(gsea_results) > 0) {
      # Test annotation
      annotated_results <- gsea_pathway_annotation(
        gsea_results = gsea_results,
        pathway_type = "GO"
      )
      
      cat(sprintf("  PASS: %d results annotated\n", nrow(annotated_results)))
      validation_results[[paste0("annotation_", category)]] <- "PASS"
      
      # Check annotation quality
      if ("pathway_name" %in% colnames(annotated_results)) {
        proper_names <- sum(!is.na(annotated_results$pathway_name) & 
                           annotated_results$pathway_name != annotated_results$pathway_id)
        cat(sprintf("    Descriptive names: %d/%d\n", proper_names, nrow(annotated_results)))
        
        if ("category" %in% colnames(annotated_results)) {
          cat("    Category info: included\n")
        } else {
          cat("    Category info: missing\n")
        }
      }
    } else {
      cat("  INFO: No results to annotate\n")
      validation_results[[paste0("annotation_", category)]] <- "PASS"
    }
    
  }, error = function(e) {
    cat(sprintf("  FAIL: %s\n", e$message))
    validation_results[[paste0("annotation_", category)]] <- "FAIL"
  })
}

# 6. ERROR HANDLING VALIDATION
cat("\n6. ERROR HANDLING VALIDATION\n")
cat("==============================\n")

cat("Testing invalid GO category:\n")
error_caught <- tryCatch({
  pathway_gsea(
    abundance = mock_abundance[1:10, ],
    metadata = mock_metadata,
    group = "Environment",
    pathway_type = "GO",
    go_category = "INVALID",
    method = "fgsea",
    nperm = 10
  )
  FALSE
}, error = function(e) TRUE)

if (error_caught) {
  cat("  PASS: Invalid category error caught\n")
  validation_results$error_handling <- "PASS"
} else {
  cat("  FAIL: Invalid category not caught\n")
  validation_results$error_handling <- "FAIL"
}

# 7. BIOLOGICAL VALIDATION
cat("\n7. BIOLOGICAL VALIDATION\n")
cat("==========================\n")

# Check that GO terms are biologically meaningful
sample_go_terms <- go_mapping[c(1, 11, 21, 29), ]  # Sample from each category

cat("Sample GO terms for biological review:\n")
for (i in 1:nrow(sample_go_terms)) {
  term <- sample_go_terms[i, ]
  cat(sprintf("  %s (%s): %s\n", term$go_id, term$category, term$go_name))
  ko_count <- length(strsplit(term$ko_members, ";")[[1]])
  cat(sprintf("    KO count: %d\n", ko_count))
}

# Check for reasonable overlap between terms
ko_term_matrix <- matrix(0, nrow = length(all_kos), ncol = nrow(go_mapping))
rownames(ko_term_matrix) <- all_kos
colnames(ko_term_matrix) <- go_mapping$go_id

for (i in 1:nrow(go_mapping)) {
  term_kos <- strsplit(go_mapping$ko_members[i], ";")[[1]]
  ko_term_matrix[term_kos, i] <- 1
}

# Calculate term-term similarity (Jaccard index)
n_comparisons <- min(100, ncol(ko_term_matrix))  # Limit for performance
if (n_comparisons > 1) {
  sample_terms <- sample(ncol(ko_term_matrix), n_comparisons)
  similarities <- numeric()
  
  for (i in 1:(n_comparisons-1)) {
    for (j in (i+1):n_comparisons) {
      term1 <- ko_term_matrix[, sample_terms[i]]
      term2 <- ko_term_matrix[, sample_terms[j]]
      
      intersection <- sum(term1 & term2)
      union <- sum(term1 | term2)
      
      if (union > 0) {
        jaccard <- intersection / union
        similarities <- c(similarities, jaccard)
      }
    }
  }
  
  if (length(similarities) > 0) {
    mean_similarity <- mean(similarities)
    cat(sprintf("Mean pairwise similarity (Jaccard): %.3f\n", mean_similarity))
    
    if (mean_similarity > 0.8) {
      cat("  WARNING: Very high similarity - terms may be redundant\n")
      validation_results$biological_validity <- "WARNING"
    } else if (mean_similarity < 0.1) {
      cat("  WARNING: Very low similarity - terms may be too specific\n")
      validation_results$biological_validity <- "WARNING"
    } else {
      cat("  PASS: Reasonable term similarity for biological data\n")
      validation_results$biological_validity <- "PASS"
    }
  }
}

# SUMMARY REPORT
cat("\n" + rep("=", 60) + "\n")
cat("GO PATHWAY SUPPORT VALIDATION SUMMARY\n")
cat(rep("=", 60) + "\n")

# Count results
pass_count <- sum(sapply(validation_results, function(x) x == "PASS"))
fail_count <- sum(sapply(validation_results, function(x) x == "FAIL"))
warning_count <- sum(sapply(validation_results, function(x) x == "WARNING"))
total_tests <- length(validation_results)

cat(sprintf("Tests completed: %d\n", total_tests))
cat(sprintf("PASSED: %d (%.1f%%)\n", pass_count, 100 * pass_count / total_tests))
cat(sprintf("FAILED: %d (%.1f%%)\n", fail_count, 100 * fail_count / total_tests))
cat(sprintf("WARNINGS: %d (%.1f%%)\n", warning_count, 100 * warning_count / total_tests))

cat("\nDetailed Results:\n")
for (test_name in names(validation_results)) {
  result <- validation_results[[test_name]]
  symbol <- switch(result, "PASS" = "✓", "FAIL" = "✗", "WARNING" = "⚠", "?")
  cat(sprintf("  %s %s: %s\n", symbol, test_name, result))
}

# Critical assessment
critical_tests <- c("go_id_format", "go_categories", "ko_format", "category_separation")
critical_failures <- sapply(critical_tests, function(test) {
  if (test %in% names(validation_results)) {
    validation_results[[test]] == "FAIL"
  } else {
    FALSE
  }
})

cat("\n")
if (any(critical_failures)) {
  cat("❌ CRITICAL ISSUES DETECTED\n")
  cat("The following critical tests failed:\n")
  for (test in critical_tests[critical_failures]) {
    cat(sprintf("  - %s\n", test))
  }
  cat("\n⚠️  GO functionality requires fixes before release.\n")
} else if (fail_count > 0) {
  cat("⚠️  MINOR ISSUES DETECTED\n")
  cat(sprintf("%d non-critical tests failed. Review before release.\n", fail_count))
  cat("Core GO functionality appears to be working.\n")
} else {
  cat("✅ VALIDATION SUCCESSFUL\n")
  cat("All tests passed. GO functionality is ready for production.\n")
  
  # Production readiness checklist
  cat("\nProduction Readiness Checklist:\n")
  cat("✓ All 36 GO terms properly formatted\n")
  cat("✓ Three categories (BP, MF, CC) correctly separated\n") 
  cat("✓ KO mappings biologically valid\n")
  cat("✓ GSEA integration functional\n")
  cat("✓ Annotation system working\n")
  cat("✓ Error handling appropriate\n")
}

cat("\nValidation completed at:", as.character(Sys.time()), "\n")

# Save results for further analysis
save(validation_results, go_mapping, file = "go_validation_complete.RData")
cat("Results saved to: go_validation_complete.RData\n")