#!/usr/bin/env Rscript
# Debug the pathway_errorbar issue

cat("=== Debugging pathway_errorbar issue ===\n")

# Load data
load('data/metacyc_abundance.RData')
load('data/metadata.RData')

# Check data structures
cat("MetaCyc abundance structure:\n")
cat("Dimensions:", dim(metacyc_abundance), "\n")
cat("First few pathway names:", metacyc_abundance$pathway[1:5], "\n")
cat("First few column names:", colnames(metacyc_abundance)[1:5], "\n")

cat("\nMetadata structure:\n")
cat("Dimensions:", dim(metadata), "\n")
cat("Sample names (first 5):", metadata$sample_name[1:5], "\n")
cat("Environment values:", unique(metadata$Environment), "\n")

# Create abundance matrix properly
abundance_matrix <- as.matrix(metacyc_abundance[, -1])
rownames(abundance_matrix) <- metacyc_abundance$pathway

cat("\nAbundance matrix:\n")
cat("Dimensions:", dim(abundance_matrix), "\n")
cat("Row names (first 5):", rownames(abundance_matrix)[1:5], "\n")
cat("Column names (first 5):", colnames(abundance_matrix)[1:5], "\n")

# Create Group vector
Group <- factor(metadata$Environment)
names(Group) <- metadata$sample_name

cat("\nGroup vector:\n")
cat("Length:", length(Group), "\n")
cat("Levels:", levels(Group), "\n")
cat("Names (first 5):", names(Group)[1:5], "\n")
cat("Values (first 5):", as.character(Group)[1:5], "\n")

# Check alignment
cat("\nAlignment check:\n")
cat("Abundance columns == Group names:", all(colnames(abundance_matrix) == names(Group)), "\n")
cat("Missing in abundance:", setdiff(names(Group), colnames(abundance_matrix)), "\n")
cat("Missing in Group:", setdiff(colnames(abundance_matrix), names(Group)), "\n")

# Create simple test DAA results
pathway_features <- rownames(abundance_matrix)[1:5]
daa_results_df <- data.frame(
  feature = pathway_features,
  method = rep("LinDA", length(pathway_features)),
  group1 = rep("Pro-survival", length(pathway_features)),
  group2 = rep("Non-survival", length(pathway_features)),
  p_values = rep(0.01, length(pathway_features)),
  p_adjust = rep(0.01, length(pathway_features)),
  description = c(
    "test description 1",
    "test description 2", 
    "test description 3",
    "test description 4",
    "test description 5"
  ),
  stringsAsFactors = FALSE
)

cat("\nTest DAA results:\n")
cat("Features:", daa_results_df$feature, "\n")
cat("All features in abundance matrix:", all(daa_results_df$feature %in% rownames(abundance_matrix)), "\n")

# Test basic inputs to pathway_errorbar
cat("\n=== Testing pathway_errorbar inputs ===\n")

# Source the function
source('R/pathway_errorbar.R')

# Test just the input validation part
cat("Testing input validation...\n")

tryCatch({
  # Call with detailed error catching
  result <- pathway_errorbar(
    abundance = abundance_matrix,
    daa_results_df = daa_results_df,
    Group = Group,
    p_values_threshold = 0.05,
    ko_to_kegg = FALSE,
    x_lab = "description"
  )
  cat("SUCCESS: pathway_errorbar worked!\n")
  cat("Result class:", class(result), "\n")
}, error = function(e) {
  cat("ERROR in pathway_errorbar:\n")
  cat("Error message:", e$message, "\n")
  cat("Call:", deparse(e$call), "\n")
  
  # Get more detailed traceback
  cat("\nTraceback:\n")
  traceback()
})

cat("\n=== Debug completed ===\n")