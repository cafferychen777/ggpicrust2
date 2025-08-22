#!/usr/bin/env Rscript
# Debug data preprocessing issue

library(ggpicrust2)

# Load data
data("ko_abundance")
data("metadata")

cat("Original KO abundance structure:\n")
cat("Dimensions:", dim(ko_abundance), "\n")
cat("First column name:", colnames(ko_abundance)[1], "\n")
cat("First few KO IDs from first column:\n")
print(head(ko_abundance[,1]))
cat("Current rownames:", head(rownames(ko_abundance)), "\n")

# Correct way to handle #NAME column
if (colnames(ko_abundance)[1] == "#NAME") {
  # Convert to data.frame and set rownames to the KO IDs before removing the column
  ko_abundance <- as.data.frame(ko_abundance)
  rownames(ko_abundance) <- ko_abundance[, 1]
  ko_abundance <- ko_abundance[, -1]
}

cat("\nAfter processing:\n")
cat("Dimensions:", dim(ko_abundance), "\n")
cat("Rownames (KO IDs):", head(rownames(ko_abundance)), "\n")
cat("Column names (samples):", head(colnames(ko_abundance)), "\n")

# Create proper matched data
n_samples <- 10
ko_abundance_subset <- ko_abundance[, 1:n_samples, drop = FALSE]
sample_ids <- colnames(ko_abundance_subset)

metadata_subset <- data.frame(
  Environment = rep(c("Group1", "Group2"), length.out = n_samples),
  row.names = sample_ids
)

cat("\nTest data created:\n")
cat("Abundance subset rownames:", head(rownames(ko_abundance_subset)), "\n")
cat("Metadata subset rownames:", head(rownames(metadata_subset)), "\n")
cat("Sample overlap:", length(intersect(colnames(ko_abundance_subset), rownames(metadata_subset))), "\n")

# Test ranking calculation
rank_result <- ggpicrust2:::calculate_rank_metric(
  abundance = as.matrix(ko_abundance_subset),
  metadata = metadata_subset,
  group = "Environment",
  method = "signal2noise"
)

cat("\nRanking results:\n")
cat("Length:", length(rank_result), "\n")
cat("Has names:", !is.null(names(rank_result)), "\n")
cat("First few names:", head(names(rank_result)), "\n")
cat("First few values:", head(rank_result), "\n")
cat("Non-finite values:", sum(!is.finite(rank_result)), "\n")

# Test if this works with fgsea
if (!is.null(names(rank_result)) && sum(is.finite(rank_result)) > 0) {
  cat("\nTesting with fgsea...\n")
  
  gene_sets <- ggpicrust2:::prepare_gene_sets("KEGG")
  
  if (requireNamespace("fgsea", quietly = TRUE)) {
    fgsea_result <- tryCatch({
      fgsea::fgsea(
        pathways = gene_sets,
        stats = rank_result,
        minSize = 5,
        maxSize = 500
      )
    }, error = function(e) {
      cat("FGSEA failed:", e$message, "\n")
      NULL
    })
    
    if (!is.null(fgsea_result)) {
      cat("FGSEA SUCCESS! Results:", nrow(fgsea_result), "pathways\n")
    }
  }
}