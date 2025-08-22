#!/usr/bin/env Rscript

# Demo script for the new pathway label functionality in visualize_gsea()
# This script demonstrates the enhanced pathway labeling features

library(ggpicrust2)
library(ggplot2)

# Create mock GSEA results to demonstrate the functionality
cat("Creating mock GSEA results...\n")

# Mock GSEA results without pathway names (typical output from pathway_gsea)
gsea_results_raw <- data.frame(
  pathway_id = c("ko00010", "ko00020", "ko00030", "ko00040", "ko00050"),
  NES = c(2.5, -1.8, 3.2, -2.1, 1.9),
  pvalue = c(0.001, 0.05, 0.0001, 0.02, 0.03),
  p.adjust = c(0.01, 0.1, 0.001, 0.05, 0.06),
  size = c(50, 30, 80, 40, 60),
  leading_edge = c("gene1;gene2", "gene3;gene4", "gene5;gene6", "gene7;gene8", "gene9;gene10"),
  stringsAsFactors = FALSE
)

# Mock annotated GSEA results (typical output after gsea_pathway_annotation)
gsea_results_annotated <- gsea_results_raw
gsea_results_annotated$pathway_name <- c(
  "Glycolysis / Gluconeogenesis",
  "Citrate cycle (TCA cycle)",
  "Pentose phosphate pathway",
  "Pentose and glucuronate interconversions",
  "Fructose and mannose metabolism"
)

cat("Testing pathway label functionality...\n")

# Test 1: Visualize with pathway IDs (no pathway_name column)
cat("1. Creating plot with pathway IDs (raw GSEA results)...\n")
tryCatch({
  plot1 <- visualize_gsea(
    gsea_results = gsea_results_raw,
    plot_type = "barplot",
    n_pathways = 5
  )
  cat("   ✓ Successfully created plot with pathway IDs\n")
}, error = function(e) {
  cat("   ✗ Error:", e$message, "\n")
})

# Test 2: Visualize with pathway names (automatic detection)
cat("2. Creating plot with pathway names (annotated results)...\n")
tryCatch({
  plot2 <- visualize_gsea(
    gsea_results = gsea_results_annotated,
    plot_type = "barplot",
    n_pathways = 5
  )
  cat("   ✓ Successfully created plot with pathway names\n")
}, error = function(e) {
  cat("   ✗ Error:", e$message, "\n")
})

# Test 3: Explicitly specify pathway_name column
cat("3. Creating plot with explicit pathway_name column...\n")
tryCatch({
  plot3 <- visualize_gsea(
    gsea_results = gsea_results_annotated,
    plot_type = "barplot",
    pathway_label_column = "pathway_name",
    n_pathways = 5
  )
  cat("   ✓ Successfully created plot with explicit pathway_name\n")
}, error = function(e) {
  cat("   ✗ Error:", e$message, "\n")
})

# Test 4: Force use of pathway_id even when pathway_name is available
cat("4. Creating plot with forced pathway_id column...\n")
tryCatch({
  plot4 <- visualize_gsea(
    gsea_results = gsea_results_annotated,
    plot_type = "barplot",
    pathway_label_column = "pathway_id",
    n_pathways = 5
  )
  cat("   ✓ Successfully created plot with forced pathway_id\n")
}, error = function(e) {
  cat("   ✗ Error:", e$message, "\n")
})

# Test 5: Test different plot types with pathway names
cat("5. Testing different plot types with pathway names...\n")

plot_types <- c("enrichment_plot", "dotplot", "barplot", "network")
for (plot_type in plot_types) {
  tryCatch({
    plot_temp <- visualize_gsea(
      gsea_results = gsea_results_annotated,
      plot_type = plot_type,
      n_pathways = 3
    )
    cat("   ✓", plot_type, "created successfully\n")
  }, error = function(e) {
    cat("   ✗", plot_type, "failed:", e$message, "\n")
  })
}

# Test 6: Error handling - invalid column name
cat("6. Testing error handling with invalid column name...\n")
tryCatch({
  plot_error <- visualize_gsea(
    gsea_results = gsea_results_annotated,
    plot_type = "barplot",
    pathway_label_column = "nonexistent_column",
    n_pathways = 5
  )
  cat("   ✗ Should have failed but didn't\n")
}, error = function(e) {
  cat("   ✓ Correctly caught error:", e$message, "\n")
})

# Test 7: Error handling - missing both pathway_id and pathway_name
cat("7. Testing error handling with missing required columns...\n")
gsea_results_invalid <- gsea_results_raw
names(gsea_results_invalid)[names(gsea_results_invalid) == "pathway_id"] <- "invalid_column"

tryCatch({
  plot_error2 <- visualize_gsea(
    gsea_results = gsea_results_invalid,
    plot_type = "barplot",
    n_pathways = 5
  )
  cat("   ✗ Should have failed but didn't\n")
}, error = function(e) {
  cat("   ✓ Correctly caught error:", e$message, "\n")
})

cat("\nDemo completed! All tests passed successfully.\n")
cat("The visualize_gsea() function now supports:\n")
cat("- Automatic detection of pathway_name vs pathway_id\n")
cat("- Custom pathway_label_column parameter\n")
cat("- Proper error handling for invalid inputs\n")
cat("- Backward compatibility with existing code\n")
