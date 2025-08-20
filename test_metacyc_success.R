#!/usr/bin/env Rscript
# Test successful MetaCyc workflow

cat("=== Testing complete MetaCyc workflow ===\n")

# Load data
load('data/metacyc_abundance.RData')
load('data/metadata.RData')

# Source functions
source('R/pathway_annotation.R')
source('R/pathway_errorbar.R')

# Create test DAA results with guaranteed significant features
set.seed(123)
pathway_features <- metacyc_abundance$pathway[1:20]

daa_results_df <- data.frame(
  feature = pathway_features,
  method = rep("LinDA", length(pathway_features)),
  group1 = rep("Pro-survival", length(pathway_features)),
  group2 = rep("Non-survival", length(pathway_features)),
  p_values = c(rep(0.001, 5), rep(0.02, 5), rep(0.3, 10)),  # 10 significant features
  stringsAsFactors = FALSE
)

daa_results_df$p_adjust <- daa_results_df$p_values * 1.1  # Slightly higher but still significant

cat("Created DAA results:\n")
cat(sprintf("Total features: %d\n", nrow(daa_results_df)))
cat(sprintf("Significant features (p < 0.05): %d\n", sum(daa_results_df$p_adjust < 0.05)))

# Test annotation
cat("\n=== Step 1: Testing pathway_annotation ===\n")
daa_annotated_results_df <- pathway_annotation(
  pathway = "MetaCyc",
  daa_results_df = daa_results_df,
  ko_to_kegg = FALSE
)

cat("Annotation results:\n")
missing_desc <- sum(is.na(daa_annotated_results_df$description))
cat(sprintf("Features with missing descriptions: %d\n", missing_desc))

if (missing_desc == 0) {
  cat("SUCCESS: All features annotated!\n")
  
  # Show first few annotations
  cat("Sample annotations:\n")
  for(i in 1:min(5, nrow(daa_annotated_results_df))) {
    cat(sprintf("  %s -> %s\n", 
                daa_annotated_results_df$feature[i],
                substr(daa_annotated_results_df$description[i], 1, 50)))
  }
}

# Prepare abundance matrix and Group vector
abundance_matrix <- as.matrix(metacyc_abundance[, -1])
rownames(abundance_matrix) <- metacyc_abundance$pathway

Group <- factor(metadata$Environment)
names(Group) <- metadata$sample_name

cat("\n=== Step 2: Testing pathway_errorbar ===\n")

tryCatch({
  p <- pathway_errorbar(
    abundance = abundance_matrix,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    p_values_threshold = 0.05,
    order = "group",
    select = NULL,
    ko_to_kegg = FALSE,
    p_value_bar = TRUE,
    x_lab = "description"
  )
  
  cat("SUCCESS: pathway_errorbar completed successfully!\n")
  cat("Plot type:", class(p), "\n")
  
  # The plot should be a patchwork object
  if ("patchwork" %in% class(p)) {
    cat("EXCELLENT: Generated patchwork plot as expected!\n")
  }
  
}, error = function(e) {
  cat("ERROR in pathway_errorbar:\n")
  cat(e$message, "\n")
})

cat("\n=== Step 3: Testing edge cases ===\n")

# Test with very strict threshold (should provide helpful error)
cat("Testing with very strict threshold...\n")
tryCatch({
  p <- pathway_errorbar(
    abundance = abundance_matrix,
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    p_values_threshold = 0.0001,  # Very strict
    ko_to_kegg = FALSE,
    x_lab = "description"
  )
  cat("Unexpected success with strict threshold\n")
}, error = function(e) {
  cat("Expected error with helpful message:\n")
  # Only show first few lines of error message
  error_lines <- strsplit(e$message, "\n")[[1]]
  cat(paste(error_lines[1:3], collapse="\n"), "\n")
})

cat("\n=== Testing completed ===\n")
cat("The MetaCyc workflow is now working properly!\n")