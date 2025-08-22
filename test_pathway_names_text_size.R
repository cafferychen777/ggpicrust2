#!/usr/bin/env Rscript

# Test script for the new pathway_names_text_size parameter
# This script demonstrates the new functionality for controlling pathway names text size

library(ggpicrust2)
library(dplyr)

# Load example data
data("kegg_abundance")
data("metadata")
data("daa_annotated_results_df")

# Prepare the data
Group <- metadata$Environment
names(Group) <- metadata$sample_name

# Filter to a subset for clearer demonstration
selected_pathways <- daa_annotated_results_df %>%
  arrange(p_adjust) %>%
  slice(1:10) %>%
  pull(feature)

# Test 1: Default behavior (auto text size)
cat("Creating plot with default pathway names text size (auto)...\n")
p1 <- pathway_errorbar(
  abundance = kegg_abundance,
  daa_results_df = daa_annotated_results_df,
  Group = Group,
  ko_to_kegg = TRUE,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = selected_pathways,
  x_lab = "pathway_name",
  pathway_names_text_size = "auto"  # Default
)

# Test 2: Small text size
cat("Creating plot with small pathway names text size (8)...\n")
p2 <- pathway_errorbar(
  abundance = kegg_abundance,
  daa_results_df = daa_annotated_results_df,
  Group = Group,
  ko_to_kegg = TRUE,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = selected_pathways,
  x_lab = "pathway_name",
  pathway_names_text_size = 8
)

# Test 3: Large text size
cat("Creating plot with large pathway names text size (14)...\n")
p3 <- pathway_errorbar(
  abundance = kegg_abundance,
  daa_results_df = daa_annotated_results_df,
  Group = Group,
  ko_to_kegg = TRUE,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = selected_pathways,
  x_lab = "pathway_name",
  pathway_names_text_size = 14
)

# Test 4: Combined with pathway_class_text_size
cat("Creating plot with both pathway names and pathway class text sizes customized...\n")
p4 <- pathway_errorbar(
  abundance = kegg_abundance,
  daa_results_df = daa_annotated_results_df,
  Group = Group,
  ko_to_kegg = TRUE,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = selected_pathways,
  x_lab = "pathway_name",
  pathway_names_text_size = 12,
  pathway_class_text_size = 6
)

# Save plots for comparison
cat("Saving plots for comparison...\n")

# Create output directory if it doesn't exist
if (!dir.exists("pathway_names_text_size_test_results")) {
  dir.create("pathway_names_text_size_test_results")
}

# Save plots
ggsave("pathway_names_text_size_test_results/01_default_auto.pdf", p1, width = 12, height = 8)
ggsave("pathway_names_text_size_test_results/02_small_size_8.pdf", p2, width = 12, height = 8)
ggsave("pathway_names_text_size_test_results/03_large_size_14.pdf", p3, width = 12, height = 8)
ggsave("pathway_names_text_size_test_results/04_combined_custom.pdf", p4, width = 12, height = 8)

cat("Test completed! Check the 'pathway_names_text_size_test_results' directory for output plots.\n")
cat("The plots demonstrate:\n")
cat("1. Default auto text size\n")
cat("2. Small text size (8)\n")
cat("3. Large text size (14)\n")
cat("4. Combined customization of both pathway names and pathway class text sizes\n")
