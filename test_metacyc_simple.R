#!/usr/bin/env Rscript
# Simplified test to reproduce MetaCyc issue without external dependencies

# Load the data
load('data/metacyc_abundance.RData')
load('data/metadata.RData')

# Check if we can load pathway_annotation reference
cat("=== Testing MetaCyc annotation reference ===\n")

# Try to load reference data
tryCatch({
  load('inst/extdata/MetaCyc_reference.RData')
  cat("MetaCyc reference loaded successfully.\n")
  cat("Reference structure:\n")
  print(str(MetaCyc_reference))
  cat("First few entries:\n")
  print(head(MetaCyc_reference))
  
  # Check column names
  cat(sprintf("Column names: %s\n", paste(colnames(MetaCyc_reference), collapse=", ")))
  
}, error = function(e) {
  cat("ERROR loading MetaCyc reference:\n")
  print(e$message)
})

# Test annotation function directly
cat("\n=== Testing annotation logic ===\n")

# Create minimal test data similar to daa_results
test_features <- c("1CMET2-PWY", "ALL-CHORISMATE-PWY", "ANAEROFRUCAT-PWY", "NONEXISTENT-PWY")

# Simulate annotation matching
matches <- match(test_features, MetaCyc_reference$X1)
descriptions <- rep(NA_character_, length(test_features))
valid_matches <- !is.na(matches)

if (any(valid_matches)) {
  descriptions[valid_matches] <- MetaCyc_reference$X2[matches[valid_matches]]
}

cat("Test results:\n")
for(i in 1:length(test_features)) {
  cat(sprintf("Feature: %s -> Description: %s\n", test_features[i], 
              ifelse(is.na(descriptions[i]), "NOT FOUND", descriptions[i])))
}

cat("\n=== Test column name standardization ===\n")
ref_data <- MetaCyc_reference

# Test the standardization logic from pathway_annotation.R
if (all(c("X1", "X2") %in% colnames(ref_data))) {
  cat("Found X1, X2 columns - standardizing to id, description\n")
  colnames(ref_data) <- c("id", "description")
  cat("New column names:", paste(colnames(ref_data), collapse=", "), "\n")
} else {
  cat("Columns already standardized or different format\n")
}

cat("\n=== Test completed ===\n")