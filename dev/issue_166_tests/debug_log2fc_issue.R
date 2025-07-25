# Debug script to identify the log2_fold_change inconsistency issue

# Load the package
devtools::load_all()

# Load example data
data("ko_abundance")
data("metadata")

print("=== Debugging log2_fold_change Inconsistency ===")

# Convert KO abundance to KEGG pathways
kegg_abundance <- ko2kegg_abundance(data = ko_abundance)

# Get DAA results
daa_results <- pathway_daa(
  abundance = kegg_abundance,
  metadata = metadata,
  group = "Environment",
  daa_method = "ALDEx2"
)

# Filter for single method
daa_single_method <- daa_results[daa_results$method == "ALDEx2_Welch's t test", ]

print("\n1. Checking metadata and Group vector consistency")
print("Original metadata:")
print(head(metadata))
print("\nGroup vector used in pathway_errorbar_table:")
Group <- metadata$Environment
print(Group)
print("\nSample names in abundance data:")
print(colnames(kegg_abundance)[1:10])

print("\n2. Testing calculate_abundance_stats directly with both approaches")

# Test feature ko00010
test_feature <- "ko00010"
group1_name <- "Pro-inflammatory"
group2_name <- "Pro-survival"

print(paste("\nTesting feature:", test_feature))
print(paste("Group1:", group1_name, "Group2:", group2_name))

# Method 1 approach (as used in pathway_daa)
print("\n--- Method 1 approach (pathway_daa) ---")
stats1 <- calculate_abundance_stats(
  abundance = kegg_abundance,
  metadata = metadata,
  group = "Environment",
  features = test_feature,
  group1 = group1_name,
  group2 = group2_name
)
print("Method 1 result:")
print(stats1)

# Method 2 approach (as used in pathway_errorbar_table)
print("\n--- Method 2 approach (pathway_errorbar_table) ---")
metadata2 <- data.frame(
  sample = colnames(kegg_abundance),
  group_col = Group
)
print("Constructed metadata for Method 2:")
print(head(metadata2))

stats2 <- calculate_abundance_stats(
  abundance = kegg_abundance,
  metadata = metadata2,
  group = "group_col",
  features = test_feature,
  group1 = group1_name,
  group2 = group2_name
)
print("Method 2 result:")
print(stats2)

print("\n3. Detailed comparison")
print(paste("Method 1 log2FC:", stats1$log2_fold_change))
print(paste("Method 2 log2FC:", stats2$log2_fold_change))
print(paste("Difference:", abs(stats1$log2_fold_change - stats2$log2_fold_change)))

print("\n4. Checking sample-group mapping")
print("Method 1 metadata sample-group mapping:")
method1_mapping <- metadata[, c("sample_name", "Environment")]
print(head(method1_mapping))

print("Method 2 metadata sample-group mapping:")
method2_mapping <- metadata2[, c("sample", "group_col")]
print(head(method2_mapping))

# Check if the mappings are identical
print("\n5. Verifying mapping consistency")
# Reorder method2 to match method1 order
method1_ordered <- method1_mapping[order(method1_mapping$sample_name), ]
method2_ordered <- method2_mapping[order(method2_mapping$sample), ]

print("Are sample names identical?")
print(identical(method1_ordered$sample_name, method2_ordered$sample))

print("Are group assignments identical?")
print(identical(as.character(method1_ordered$Environment), as.character(method2_ordered$group_col)))

print("\n=== Conclusion ===")
if (abs(stats1$log2_fold_change - stats2$log2_fold_change) < 1e-6) {
  print("✅ Methods are consistent - issue might be elsewhere")
} else {
  print("❌ Methods are inconsistent - found the root cause!")
  print("Need to investigate further...")
}
