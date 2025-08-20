#!/usr/bin/env Rscript
# Test script to reproduce the MetaCyc pathway_errorbar issue

# Load necessary libraries
library(dplyr)

# Load the data
load('data/metacyc_abundance.RData')
load('data/metadata.RData')

# Source the functions
source('R/pathway_daa.R')
source('R/pathway_annotation.R') 
source('R/pathway_errorbar.R')

# Print data structure for debugging
cat("MetaCyc abundance data structure:\n")
print(str(metacyc_abundance))
print(head(metacyc_abundance[,1:5]))

cat("\nMetadata structure:\n")
print(str(metadata))
print(head(metadata))

# Step 1: Perform differential abundance analysis
cat("\n=== Step 1: DAA Analysis ===\n")
print("Running pathway_daa...")

daa_results_df <- pathway_daa(
  abundance = metacyc_abundance %>% column_to_rownames("pathway"),
  metadata = metadata,
  group = "Environment",
  daa_method = "LinDA",
  select = NULL,
  reference = NULL
)

cat("DAA results structure:\n")
print(str(daa_results_df))
print(head(daa_results_df))

# Step 2: Annotation
cat("\n=== Step 2: Annotation ===\n")
print("Running pathway_annotation...")

daa_annotated_results_df <- pathway_annotation(
  pathway = "MetaCyc",
  daa_results_df = daa_results_df,
  ko_to_kegg = FALSE
)

cat("Annotated results structure:\n")
print(str(daa_annotated_results_df))
print(head(daa_annotated_results_df))

# Check for missing descriptions
missing_descriptions <- sum(is.na(daa_annotated_results_df$description))
cat(sprintf("Features with missing descriptions: %d\n", missing_descriptions))

# Check for significant features
sig_features <- sum(daa_annotated_results_df$p_adjust < 0.05, na.rm = TRUE)
cat(sprintf("Significant features (p < 0.05): %d\n", sig_features))

# Step 3: Create Group vector
Group <- metadata$Environment
names(Group) <- metadata$sample_name

# Step 4: Try pathway_errorbar
cat("\n=== Step 3: pathway_errorbar ===\n")
print("Running pathway_errorbar...")

tryCatch({
  p <- pathway_errorbar(
    abundance = metacyc_abundance %>% column_to_rownames("pathway"),
    daa_results_df = daa_annotated_results_df,
    Group = Group,
    p_values_threshold = 0.05,
    order = "group",
    select = NULL,
    ko_to_kegg = FALSE,
    p_value_bar = TRUE,
    colors = NULL,
    x_lab = "description"
  )
  cat("SUCCESS: pathway_errorbar completed successfully!\n")
  print(class(p))
}, error = function(e) {
  cat("ERROR in pathway_errorbar:\n")
  print(e$message)
  
  # Additional debugging
  cat("\nDebugging info:\n")
  
  # Check filtered results
  filtered_df <- daa_annotated_results_df[daa_annotated_results_df$p_adjust < 0.05,]
  cat(sprintf("Rows after p-value filtering: %d\n", nrow(filtered_df)))
  
  # Check for missing descriptions after filtering
  if(nrow(filtered_df) > 0) {
    missing_after_filter <- sum(is.na(filtered_df$description))
    cat(sprintf("Missing descriptions after p-value filtering: %d\n", missing_after_filter))
    
    # Check after excluding missing descriptions
    final_df <- filtered_df[!is.na(filtered_df$description),]
    cat(sprintf("Final rows after excluding missing descriptions: %d\n", nrow(final_df)))
  }
})

cat("\n=== Test completed ===\n")