#!/usr/bin/env Rscript

# Demo script to show PICRUSt 2.6.2 compatibility fixes
# This script demonstrates how the improved ggpicrust2 handles problematic data

library(ggpicrust2)
library(dplyr)

cat("=== PICRUSt 2.6.2 Compatibility Demo ===\n\n")

# Simulate PICRUSt 2.6.2 problematic data
cat("1. Creating simulated PICRUSt 2.6.2 data with many zero-abundance features...\n")

# Create abundance data that mimics the PICRUSt 2.6.2 issue
abundance_problematic <- data.frame(
  sample1 = c(0, 0, 0, 0, 0, 10, 5, 0, 0, 0),
  sample2 = c(0, 0, 0, 0, 0, 20, 8, 0, 0, 0),
  sample3 = c(0, 0, 0, 0, 0, 15, 6, 0, 0, 0),
  sample4 = c(0, 0, 0, 0, 0, 25, 10, 0, 0, 0),
  sample5 = c(0, 0, 0, 0, 0, 12, 7, 0, 0, 0),
  sample6 = c(0, 0, 0, 0, 0, 18, 9, 0, 0, 0),
  row.names = paste0("K", sprintf("%05d", 1:10))
)

metadata <- data.frame(
  sample = paste0("sample", 1:6),
  group = c("control", "control", "control", "treatment", "treatment", "treatment"),
  stringsAsFactors = FALSE
)

cat("Data summary:\n")
cat("- Total features:", nrow(abundance_problematic), "\n")
cat("- Features with zero abundance:", sum(rowSums(abundance_problematic) == 0), "\n")
cat("- Proportion of zero features:", 
    round(sum(rowSums(abundance_problematic) == 0) / nrow(abundance_problematic) * 100, 1), "%\n\n")

# Test the improved pathway_daa function
cat("2. Testing improved pathway_daa function with ALDEx2...\n")

tryCatch({
  result_aldex2 <- pathway_daa(
    abundance = abundance_problematic,
    metadata = metadata,
    group = "group",
    daa_method = "ALDEx2"
  )
  
  cat("✓ ALDEx2 analysis completed successfully!\n")
  cat("- Results shape:", nrow(result_aldex2), "rows x", ncol(result_aldex2), "columns\n")
  cat("- Methods included:", paste(unique(result_aldex2$method), collapse = ", "), "\n\n")
  
}, error = function(e) {
  cat("✗ ALDEx2 analysis failed:", e$message, "\n\n")
})

# Test with LinDA
cat("3. Testing improved pathway_daa function with LinDA...\n")

tryCatch({
  result_linda <- pathway_daa(
    abundance = abundance_problematic,
    metadata = metadata,
    group = "group",
    daa_method = "LinDA"
  )
  
  cat("✓ LinDA analysis completed successfully!\n")
  cat("- Results shape:", nrow(result_linda), "rows x", ncol(result_linda), "columns\n\n")
  
}, error = function(e) {
  cat("✗ LinDA analysis failed:", e$message, "\n")
  cat("This is expected for very sparse data. Try ALDEx2 or DESeq2 instead.\n\n")
})

# Test ko2kegg_abundance with problematic data
cat("4. Testing improved ko2kegg_abundance function...\n")

ko_data_problematic <- data.frame(
  `#NAME` = paste0("K", sprintf("%05d", 1:5)),
  sample1 = c(0, 0, 0, 0, 0),
  sample2 = c(0, 0, 0, 0, 0),
  sample3 = c(0, 0, 0, 0, 0),
  check.names = FALSE
)

tryCatch({
  kegg_result <- ko2kegg_abundance(data = ko_data_problematic)
  cat("✓ ko2kegg_abundance completed with compatibility handling!\n")
  cat("- Output pathways:", nrow(kegg_result), "\n\n")
  
}, error = function(e) {
  cat("✗ ko2kegg_abundance failed:", e$message, "\n\n")
})

cat("=== Summary ===\n")
cat("The improved ggpicrust2 now includes:\n")
cat("1. Better data validation and filtering\n")
cat("2. Improved error messages with troubleshooting hints\n")
cat("3. Graceful handling of PICRUSt 2.6.2 format issues\n")
cat("4. Automatic detection and warning of compatibility problems\n")
cat("5. Fallback strategies for problematic data\n\n")

cat("For users experiencing PICRUSt 2.6.2 issues:\n")
cat("- Try ALDEx2 or DESeq2 instead of LinDA\n")
cat("- Check the compatibility guide: inst/PICRUST_COMPATIBILITY.md\n")
cat("- Report persistent issues with sample data\n")
