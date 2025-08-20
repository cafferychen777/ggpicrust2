#!/usr/bin/env Rscript
# Direct test of the annotation issue

# Load data
load('data/metacyc_abundance.RData')
load('data/metadata.RData')
load('inst/extdata/MetaCyc_reference.RData')

cat("=== Direct MetaCyc annotation test ===\n")

# Use actual pathways from the MetaCyc data
pathway_features <- metacyc_abundance$pathway[1:20]

# Create simulated DAA results
daa_results_df <- data.frame(
  feature = pathway_features,
  method = rep("LinDA", length(pathway_features)),
  group1 = rep("Pro-survival", length(pathway_features)),
  group2 = rep("Non-survival", length(pathway_features)),
  p_values = runif(length(pathway_features), 0.001, 0.8),
  p_adjust = runif(length(pathway_features), 0.001, 0.8),
  stringsAsFactors = FALSE
)

cat("Testing pathways:\n")
print(daa_results_df$feature)

# Direct annotation using the reference data
cat("\n=== Manual annotation process ===\n")

# Standardize column names as done in pathway_annotation.R
ref_data <- MetaCyc_reference
if (all(c("X1", "X2") %in% colnames(ref_data))) {
  colnames(ref_data) <- c("id", "description")
}

# Match features with reference data
features <- daa_results_df$feature
matches <- match(features, ref_data$id)

# Create description column
descriptions <- rep(NA_character_, length(features))
valid_matches <- !is.na(matches)
if (any(valid_matches)) {
  descriptions[valid_matches] <- ref_data$description[matches[valid_matches]]
}

# Add descriptions to DAA results
daa_annotated_results_df <- daa_results_df
daa_annotated_results_df$description <- descriptions

cat("Annotation results:\n")
for(i in 1:nrow(daa_annotated_results_df)) {
  cat(sprintf("%s -> %s\n", 
              daa_annotated_results_df$feature[i],
              ifelse(is.na(daa_annotated_results_df$description[i]), 
                     "NOT FOUND", daa_annotated_results_df$description[i])))
}

# Check missing descriptions
missing_descriptions <- sum(is.na(daa_annotated_results_df$description))
cat(sprintf("\nFeatures with missing descriptions: %d out of %d\n", 
            missing_descriptions, nrow(daa_annotated_results_df)))

# Now test the pathway_errorbar filtering logic
cat("\n=== Testing pathway_errorbar filtering ===\n")

p_values_threshold <- 0.05
x_lab <- "description"

# Filter by p-values
daa_results_filtered_df <- daa_annotated_results_df[daa_annotated_results_df$p_adjust < p_values_threshold,]
cat(sprintf("After p-value filtering (p < %g): %d features\n", p_values_threshold, nrow(daa_results_filtered_df)))

if(nrow(daa_results_filtered_df) > 0) {
  cat("Significant features:\n")
  for(i in 1:nrow(daa_results_filtered_df)) {
    cat(sprintf("  %s: p_adjust=%.4f, description=%s\n",
                daa_results_filtered_df$feature[i],
                daa_results_filtered_df$p_adjust[i],
                ifelse(is.na(daa_results_filtered_df$description[i]), "MISSING", "present")))
  }
  
  # Apply the critical filtering step from pathway_errorbar.R line 263
  # daa_results_df <- daa_results_df[!is.na(daa_results_df[,x_lab]),]
  before_count <- nrow(daa_results_filtered_df)
  daa_results_filtered_df <- daa_results_filtered_df[!is.na(daa_results_filtered_df[,x_lab]),]
  after_count <- nrow(daa_results_filtered_df)
  
  cat(sprintf("\nAfter excluding missing %s: %d features (was %d)\n", 
              x_lab, after_count, before_count))
  
  if(after_count == 0) {
    cat("*** ISSUE REPRODUCED! ***\n")
    cat("All significant features were excluded due to missing descriptions!\n")
    cat("This would cause pathway_errorbar to fail with:\n")
    cat("'Visualization with 'pathway_errorbar' cannot be performed because there are no features with statistical significance.'\n")
  } else {
    cat("Remaining features after all filtering:\n")
    for(i in 1:nrow(daa_results_filtered_df)) {
      cat(sprintf("  %s: %s\n",
                  daa_results_filtered_df$feature[i],
                  daa_results_filtered_df$description[i]))
    }
  }
} else {
  cat("No significant features found - would also cause pathway_errorbar to fail\n")
}

cat("\n=== Test completed ===\n")