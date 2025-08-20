#!/usr/bin/env Rscript
# Manual test to reproduce the exact MetaCyc pathway_errorbar issue

# Source required functions first (basic R way)
source('R/pathway_annotation.R')

# Load data
load('data/metacyc_abundance.RData')
load('data/metadata.RData')

cat("=== MetaCyc pathway_errorbar reproduction test ===\n")

# Step 1: Create simulated DAA results (since we can't run pathway_daa without dependencies)
cat("Creating simulated DAA results...\n")

# Use actual pathways from the MetaCyc data
pathway_features <- metacyc_abundance$pathway[1:20]  # First 20 pathways

# Create realistic DAA results structure
daa_results_df <- data.frame(
  feature = pathway_features,
  method = rep("LinDA", length(pathway_features)),
  group1 = rep("Pro-survival", length(pathway_features)),
  group2 = rep("Non-survival", length(pathway_features)),
  p_values = runif(length(pathway_features), 0.001, 0.8),
  p_adjust = runif(length(pathway_features), 0.001, 0.8),
  stringsAsFactors = FALSE
)

cat("DAA results structure:\n")
print(str(daa_results_df))

# Step 2: Test pathway annotation
cat("\n=== Testing pathway_annotation function ===\n")

tryCatch({
  daa_annotated_results_df <- pathway_annotation(
    pathway = "MetaCyc",
    daa_results_df = daa_results_df,
    ko_to_kegg = FALSE
  )
  
  cat("Annotation completed successfully!\n")
  cat("Annotated structure:\n")
  print(str(daa_annotated_results_df))
  
  # Check for missing descriptions
  missing_descriptions <- sum(is.na(daa_annotated_results_df$description))
  cat(sprintf("Features with missing descriptions: %d out of %d\n", 
              missing_descriptions, nrow(daa_annotated_results_df)))
  
  # Show which features have descriptions
  cat("\nAnnotation results:\n")
  for(i in 1:min(10, nrow(daa_annotated_results_df))) {
    cat(sprintf("%s -> %s\n", 
                daa_annotated_results_df$feature[i],
                ifelse(is.na(daa_annotated_results_df$description[i]), 
                       "MISSING", daa_annotated_results_df$description[i])))
  }
  
}, error = function(e) {
  cat("ERROR in pathway_annotation:\n")
  print(e$message)
  stop("Cannot proceed without annotation")
})

# Step 3: Simulate the pathway_errorbar filtering logic
cat("\n=== Testing pathway_errorbar filtering logic ===\n")

# Key parameters from pathway_errorbar
p_values_threshold <- 0.05
x_lab <- "description"  # Default for ko_to_kegg = FALSE

cat(sprintf("Using p_values_threshold: %g\n", p_values_threshold))
cat(sprintf("Using x_lab: %s\n", x_lab))

# Step 3a: Filter by p-values
daa_results_filtered_df <- daa_annotated_results_df[daa_annotated_results_df$p_adjust < p_values_threshold,]
cat(sprintf("After p-value filtering: %d features\n", nrow(daa_results_filtered_df)))

# Step 3b: Exclude rows with missing pathway annotation (key step!)
if(nrow(daa_results_filtered_df) > 0) {
  # This is the critical line from pathway_errorbar.R:263
  daa_results_filtered_df <- daa_results_filtered_df[!is.na(daa_results_filtered_df[,x_lab]),]
  cat(sprintf("After excluding missing %s: %d features\n", x_lab, nrow(daa_results_filtered_df)))
  
  if(nrow(daa_results_filtered_df) == 0) {
    cat("*** THIS IS THE ISSUE! ***\n")
    cat("All features were filtered out due to missing descriptions!\n")
    
    # Show what happened
    sig_features <- daa_annotated_results_df[daa_annotated_results_df$p_adjust < p_values_threshold,]
    cat(sprintf("Significant features before description filter: %d\n", nrow(sig_features)))
    
    if(nrow(sig_features) > 0) {
      missing_desc_count <- sum(is.na(sig_features$description))
      cat(sprintf("Significant features with missing descriptions: %d\n", missing_desc_count))
      
      cat("\nDetailed breakdown:\n")
      for(i in 1:nrow(sig_features)) {
        cat(sprintf("  %s: p_adjust=%.4f, description=%s\n",
                    sig_features$feature[i],
                    sig_features$p_adjust[i],
                    ifelse(is.na(sig_features$description[i]), "MISSING", "present")))
      }
    }
  } else {
    cat("SUCCESS: Some features remain after all filtering!\n")
  }
} else {
  cat("No significant features found at p < 0.05\n")
}

cat("\n=== Test completed ===\n")