#!/usr/bin/env Rscript
# Test the fixes for MetaCyc pathway_errorbar issue

cat("=== Testing fixes for MetaCyc pathway_errorbar issue ===\n")

# Load data
load('data/metacyc_abundance.RData')
load('data/metadata.RData')

# Test 1: Test the improved load_reference_data function
cat("\n=== Test 1: Fixed load_reference_data function ===\n")

source('R/pathway_annotation.R')

tryCatch({
  ref_data <- load_reference_data("MetaCyc")
  cat("SUCCESS: Reference data loaded successfully!\n")
  cat(sprintf("Loaded %d pathways\n", nrow(ref_data)))
  cat(sprintf("Column names: %s\n", paste(colnames(ref_data), collapse=", ")))
}, error = function(e) {
  cat("ERROR: Failed to load reference data:\n")
  cat(e$message, "\n")
})

# Test 2: Test the complete annotation workflow
cat("\n=== Test 2: Complete annotation workflow ===\n")

# Create realistic DAA results
set.seed(42)
pathway_features <- metacyc_abundance$pathway[1:30]

daa_results_df <- data.frame(
  feature = pathway_features,
  method = rep("LinDA", length(pathway_features)),
  group1 = rep("Pro-survival", length(pathway_features)),
  group2 = rep("Non-survival", length(pathway_features)),
  p_values = runif(length(pathway_features), 0.001, 0.9),
  stringsAsFactors = FALSE
)

# Create p_adjust values with some significant results
daa_results_df$p_adjust <- pmax(daa_results_df$p_values * runif(length(pathway_features), 0.8, 1.2), 0.001)

cat(sprintf("Created DAA results with %d features\n", nrow(daa_results_df)))
cat(sprintf("Significant features (p < 0.05): %d\n", sum(daa_results_df$p_adjust < 0.05)))

# Test pathway_annotation with the fixed function
tryCatch({
  daa_annotated_results_df <- pathway_annotation(
    pathway = "MetaCyc",
    daa_results_df = daa_results_df,
    ko_to_kegg = FALSE
  )
  
  cat("SUCCESS: pathway_annotation completed!\n")
  
  # Check annotation results
  missing_descriptions <- sum(is.na(daa_annotated_results_df$description))
  cat(sprintf("Features with missing descriptions: %d out of %d\n", 
              missing_descriptions, nrow(daa_annotated_results_df)))
  
  if (missing_descriptions == 0) {
    cat("EXCELLENT: All features were successfully annotated!\n")
  }
  
}, error = function(e) {
  cat("ERROR in pathway_annotation:\n")
  cat(e$message, "\n")
  daa_annotated_results_df <- daa_results_df
  daa_annotated_results_df$description <- NA_character_
})

# Test 3: Test the improved pathway_errorbar error messages
cat("\n=== Test 3: Testing pathway_errorbar with improved error handling ===\n")

# Source the improved pathway_errorbar function
source('R/pathway_errorbar.R')

# Create Group vector
Group <- factor(metadata$Environment)
names(Group) <- metadata$sample_name

# Test scenario 1: With missing descriptions (should give helpful error)
cat("\nScenario 1: Testing with missing descriptions...\n")
daa_broken <- daa_results_df
daa_broken$description <- NA_character_  # Simulate broken annotation

# Prepare abundance matrix
abundance_matrix <- as.matrix(metacyc_abundance[, -1])
rownames(abundance_matrix) <- metacyc_abundance$pathway

tryCatch({
  p <- pathway_errorbar(
    abundance = abundance_matrix,
    daa_results_df = daa_broken,
    Group = Group,
    p_values_threshold = 0.05,
    ko_to_kegg = FALSE,
    x_lab = "description"
  )
  cat("Unexpected SUCCESS - this should have failed\n")
}, error = function(e) {
  cat("Expected ERROR with improved message:\n")
  cat(e$message, "\n")
})

# Test scenario 2: With working annotations (should succeed)
if (exists("daa_annotated_results_df") && "description" %in% colnames(daa_annotated_results_df)) {
  cat("\nScenario 2: Testing with working annotations...\n")
  
  # Check if we have significant features
  sig_count <- sum(daa_annotated_results_df$p_adjust < 0.05, na.rm = TRUE)
  cat(sprintf("Significant features available: %d\n", sig_count))
  
  if (sig_count > 0) {
    tryCatch({
      # Test with a more permissive threshold to ensure we have features
      p <- pathway_errorbar(
        abundance = abundance_matrix,
        daa_results_df = daa_annotated_results_df,
        Group = Group,
        p_values_threshold = 0.1,  # More permissive threshold
        ko_to_kegg = FALSE,
        x_lab = "description"
      )
      cat("SUCCESS: pathway_errorbar completed successfully!\n")
      cat("Plot class:", class(p), "\n")
    }, error = function(e) {
      cat("ERROR in pathway_errorbar:\n")
      cat(e$message, "\n")
    })
  } else {
    cat("No significant features to test with - trying relaxed threshold...\n")
    
    tryCatch({
      p <- pathway_errorbar(
        abundance = abundance_matrix,
        daa_results_df = daa_annotated_results_df,
        Group = Group,
        p_values_threshold = 0.5,  # Very permissive threshold
        ko_to_kegg = FALSE,
        x_lab = "description"
      )
      cat("SUCCESS with relaxed threshold!\n")
    }, error = function(e) {
      cat("Still failed with relaxed threshold:\n")
      cat(e$message, "\n")
    })
  }
}

cat("\n=== Testing completed ===\n")
cat("Summary of fixes:\n")
cat("1. Improved load_reference_data() with better fallback strategies\n")
cat("2. Enhanced error messages in pathway_errorbar()\n")
cat("3. Better diagnostic information when annotation fails\n")
cat("4. Clearer guidance for MetaCyc-specific issues\n")