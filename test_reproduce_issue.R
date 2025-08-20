#!/usr/bin/env Rscript
# Test to reproduce the exact issue users face

cat("=== Reproducing real-world MetaCyc issue ===\n")

# Load data
load('data/metacyc_abundance.RData')
load('data/metadata.RData')

# Create realistic DAA results (some significant, some not)
set.seed(42)  # For reproducible results
pathway_features <- metacyc_abundance$pathway[1:50]  # Use more pathways

daa_results_df <- data.frame(
  feature = pathway_features,
  method = rep("LinDA", length(pathway_features)),
  group1 = rep("Pro-survival", length(pathway_features)),
  group2 = rep("Non-survival", length(pathway_features)),
  p_values = runif(length(pathway_features), 0.001, 0.9),
  stringsAsFactors = FALSE
)

# Make p_adjust values - some significant, some not
daa_results_df$p_adjust <- pmax(daa_results_df$p_values * runif(length(pathway_features), 0.8, 1.2), 0.001)

cat("DAA results summary:\n")
cat(sprintf("Total features: %d\n", nrow(daa_results_df)))
cat(sprintf("Significant features (p < 0.05): %d\n", sum(daa_results_df$p_adjust < 0.05)))
cat(sprintf("Significant features (p < 0.1): %d\n", sum(daa_results_df$p_adjust < 0.1)))

# Scenario 1: Test with broken annotation (simulating package installation issue)
cat("\n=== Scenario 1: Broken annotation (simulates real user issue) ===\n")

# Simulate what happens when pathway_annotation fails
daa_annotated_results_df_broken <- daa_results_df
daa_annotated_results_df_broken$description <- NA_character_  # All descriptions missing!

cat("Simulating the case where pathway_annotation fails to add descriptions...\n")

# Test pathway_errorbar filtering logic with missing descriptions
p_values_threshold <- 0.05
x_lab <- "description"

# Filter by p-values first
sig_features_broken <- daa_annotated_results_df_broken[daa_annotated_results_df_broken$p_adjust < p_values_threshold,]
cat(sprintf("Significant features before description filter: %d\n", nrow(sig_features_broken)))

# Apply the critical filter that removes missing descriptions
if(nrow(sig_features_broken) > 0) {
  final_features_broken <- sig_features_broken[!is.na(sig_features_broken[,x_lab]),]
  cat(sprintf("Significant features after description filter: %d\n", nrow(final_features_broken)))
  
  if(nrow(final_features_broken) == 0) {
    cat("*** ISSUE REPRODUCED! ***\n")
    cat("This is the exact error users see: 'no features with statistical significance'\n")
  }
}

# Scenario 2: Test with working annotation
cat("\n=== Scenario 2: Working annotation ===\n")

# Load reference data manually (simulating fixed annotation)
load('inst/extdata/MetaCyc_reference.RData')
ref_data <- MetaCyc_reference
if (all(c("X1", "X2") %in% colnames(ref_data))) {
  colnames(ref_data) <- c("id", "description")
}

# Properly annotate the results
features <- daa_results_df$feature
matches <- match(features, ref_data$id)
descriptions <- rep(NA_character_, length(features))
valid_matches <- !is.na(matches)
if (any(valid_matches)) {
  descriptions[valid_matches] <- ref_data$description[matches[valid_matches]]
}

daa_annotated_results_df_working <- daa_results_df
daa_annotated_results_df_working$description <- descriptions

missing_count <- sum(is.na(descriptions))
cat(sprintf("Features with missing descriptions: %d out of %d\n", missing_count, length(descriptions)))

# Test pathway_errorbar filtering with proper annotations
sig_features_working <- daa_annotated_results_df_working[daa_annotated_results_df_working$p_adjust < p_values_threshold,]
cat(sprintf("Significant features before description filter: %d\n", nrow(sig_features_working)))

if(nrow(sig_features_working) > 0) {
  # Show which significant features have descriptions
  cat("Significant features and their annotation status:\n")
  for(i in 1:nrow(sig_features_working)) {
    desc_status <- ifelse(is.na(sig_features_working$description[i]), "MISSING", "present")
    cat(sprintf("  %s: p=%.4f, description=%s\n", 
                sig_features_working$feature[i], 
                sig_features_working$p_adjust[i], 
                desc_status))
  }
  
  final_features_working <- sig_features_working[!is.na(sig_features_working[,x_lab]),]
  cat(sprintf("Significant features after description filter: %d\n", nrow(final_features_working)))
  
  if(nrow(final_features_working) > 0) {
    cat("SUCCESS: pathway_errorbar would work with these features!\n")
  } else {
    cat("STILL FAILS: Even with annotation, no features remain after filtering\n")
  }
}

# Scenario 3: Test edge case with relaxed threshold
cat("\n=== Scenario 3: Relaxed threshold ===\n")
p_values_threshold_relaxed <- 0.1

sig_features_relaxed <- daa_annotated_results_df_working[daa_annotated_results_df_working$p_adjust < p_values_threshold_relaxed,]
cat(sprintf("Significant features with p < 0.1: %d\n", nrow(sig_features_relaxed)))

if(nrow(sig_features_relaxed) > 0) {
  final_features_relaxed <- sig_features_relaxed[!is.na(sig_features_relaxed[,x_lab]),]
  cat(sprintf("Features remaining after description filter: %d\n", nrow(final_features_relaxed)))
}

cat("\n=== Analysis Complete ===\n")
cat("ROOT CAUSE: The issue occurs when pathway_annotation fails to load reference data\n")
cat("due to package installation problems, resulting in missing descriptions.\n")
cat("SOLUTION: Fix the load_reference_data function to handle installation issues gracefully.\n")