#!/usr/bin/env Rscript

# Example demonstrating the new pathway_names_text_size parameter
# This addresses GitHub issue #173: modify pathway_class and pathway_map text size in errorbar plot

library(ggpicrust2)
library(dplyr)

# Load example data
data("kegg_abundance")
data("metadata")
data("daa_annotated_results_df")

# Prepare the data
Group <- metadata$Environment
names(Group) <- metadata$sample_name

# Select a subset of pathways for demonstration
selected_pathways <- daa_annotated_results_df %>%
  arrange(p_adjust) %>%
  slice(1:8) %>%
  pull(feature)

# Example 1: Default behavior (auto text size)
# This will automatically calculate appropriate text sizes based on the number of pathways
p1 <- pathway_errorbar(
  abundance = kegg_abundance,
  daa_results_df = daa_annotated_results_df,
  Group = Group,
  ko_to_kegg = TRUE,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = selected_pathways,
  x_lab = "pathway_name"
  # pathway_names_text_size = "auto" is the default
)

# Example 2: Increase pathway names text size for better readability
# This addresses the issue where pathway names are too small to read
p2 <- pathway_errorbar(
  abundance = kegg_abundance,
  daa_results_df = daa_annotated_results_df,
  Group = Group,
  ko_to_kegg = TRUE,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = selected_pathways,
  x_lab = "pathway_name",
  pathway_names_text_size = 12  # Increase from default (~10) to 12
)

# Example 3: Customize both pathway names and pathway class text sizes
# This gives full control over text readability
p3 <- pathway_errorbar(
  abundance = kegg_abundance,
  daa_results_df = daa_annotated_results_df,
  Group = Group,
  ko_to_kegg = TRUE,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = selected_pathways,
  x_lab = "pathway_name",
  pathway_names_text_size = 14,      # Large pathway names
  pathway_class_text_size = 6        # Large pathway class annotations
)

# Example 4: For plots with many pathways, use smaller text
many_pathways <- daa_annotated_results_df %>%
  arrange(p_adjust) %>%
  slice(1:15) %>%
  pull(feature)

p4 <- pathway_errorbar(
  abundance = kegg_abundance,
  daa_results_df = daa_annotated_results_df,
  Group = Group,
  ko_to_kegg = TRUE,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = many_pathways,
  x_lab = "pathway_name",
  pathway_names_text_size = 8,       # Smaller text for many pathways
  pathway_class_text_size = 4
)

cat("Examples created successfully!\n")
cat("New parameter 'pathway_names_text_size' allows control of y-axis pathway labels text size.\n")
cat("This addresses GitHub issue #173 about small, hard-to-read pathway text.\n")
cat("\nUsage:\n")
cat("- pathway_names_text_size = 'auto' (default): Automatically calculated size\n")
cat("- pathway_names_text_size = 12: Custom size (numeric value)\n")
cat("- Works together with existing pathway_class_text_size parameter\n")
