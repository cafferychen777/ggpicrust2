#!/usr/bin/env Rscript
# Debug FGSEA naming issue

library(ggpicrust2)

# Create test data
data("ko_abundance")
data("metadata")

# Remove #NAME column if present
if (colnames(ko_abundance)[1] == "#NAME") {
  ko_abundance <- ko_abundance[, -1]
}

# Create matched datasets
n_samples <- min(ncol(ko_abundance), nrow(metadata), 10)  # Small test
sample_ids <- colnames(ko_abundance)[1:n_samples]
metadata_subset <- metadata[1:n_samples, , drop = FALSE]
rownames(metadata_subset) <- sample_ids
ko_abundance_subset <- ko_abundance[, 1:n_samples, drop = FALSE]

# Ensure groups exist
metadata_subset$Environment <- rep(c("Group1", "Group2"), length.out = n_samples)

cat("Data preparation complete\n")
cat("Ko abundance shape:", dim(ko_abundance_subset), "\n")
cat("Metadata shape:", dim(metadata_subset), "\n")
cat("Sample overlap:", length(intersect(colnames(ko_abundance_subset), rownames(metadata_subset))), "\n")

# Test ranking calculation directly
cat("Testing ranking calculation...\n")

abundance_mat <- as.matrix(ko_abundance_subset)
group <- "Environment"
metadata <- metadata_subset

# Test calculate_rank_metric function directly
rank_result <- tryCatch({
  ggpicrust2:::calculate_rank_metric(abundance_mat, metadata, group, method = "signal2noise")
}, error = function(e) {
  cat("Ranking calculation failed:", e$message, "\n")
  NULL
})

if (!is.null(rank_result)) {
  cat("Ranking calculation successful\n")
  cat("Rank result length:", length(rank_result), "\n")
  cat("Has names:", !is.null(names(rank_result)), "\n")
  cat("First few names:", head(names(rank_result)), "\n")
  cat("First few values:", head(rank_result), "\n")
  
  # Test gene sets preparation
  cat("Testing gene sets preparation...\n")
  gene_sets <- tryCatch({
    ggpicrust2:::prepare_gene_sets("KEGG")
  }, error = function(e) {
    cat("Gene sets preparation failed:", e$message, "\n")
    NULL
  })
  
  if (!is.null(gene_sets)) {
    cat("Gene sets preparation successful\n")
    cat("Number of gene sets:", length(gene_sets), "\n")
    cat("First gene set names:", head(names(gene_sets)), "\n")
    
    # Test fgsea directly
    cat("Testing fgsea directly...\n")
    
    # Load fgsea
    if (requireNamespace("fgsea", quietly = TRUE)) {
      fgsea_result <- tryCatch({
        fgsea::fgsea(
          pathways = gene_sets,
          stats = rank_result,
          minSize = 5,
          maxSize = 500,
          nperm = 100
        )
      }, error = function(e) {
        cat("FGSEA failed:", e$message, "\n")
        e
      })
      
      if (!inherits(fgsea_result, "error")) {
        cat("FGSEA successful!\n")
        cat("Number of results:", nrow(fgsea_result), "\n")
      }
    } else {
      cat("fgsea package not available\n")
    }
  }
} else {
  cat("Cannot proceed without ranking results\n")
}