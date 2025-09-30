#!/usr/bin/env Rscript
# Fix metadata.RData to have correct rownames
# This fixes Issue #XXX where pathway_gsea fails with "0/50 overlapping samples"

cat("========================================\n")
cat("Fixing metadata.RData rownames\n")
cat("========================================\n\n")

# Load current metadata
load("data/metadata.RData")

cat("BEFORE fix:\n")
cat("  Dimensions:", dim(metadata), "\n")
cat("  First 5 rownames:", paste(rownames(metadata)[1:5], collapse=", "), "\n")
cat("  First 5 sample_name:", paste(metadata$sample_name[1:5], collapse=", "), "\n")
cat("  Rownames match sample_name?", all(rownames(metadata) == metadata$sample_name), "\n\n")

# Check if already fixed
if (all(rownames(metadata) == metadata$sample_name)) {
  cat("✅ Rownames already correct! No fix needed.\n")
  quit(status = 0)
}

# Apply fix: Set rownames to sample IDs
rownames(metadata) <- metadata$sample_name

cat("AFTER fix:\n")
cat("  First 5 rownames:", paste(rownames(metadata)[1:5], collapse=", "), "\n")
cat("  Rownames match sample_name?", all(rownames(metadata) == metadata$sample_name), "\n\n")

# Save fixed metadata
save(metadata, file = "data/metadata.RData", compress = "xz")

cat("========================================\n")
cat("✅ SUCCESS: metadata.RData fixed!\n")
cat("========================================\n\n")

cat("Now all pathway_gsea examples will work correctly.\n")
cat("Users can use: pathway_gsea(abundance = ..., metadata = metadata, ...)\n")
cat("without needing to manually set rownames.\n")