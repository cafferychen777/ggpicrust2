# Simple test to validate GO functionality

library(ggpicrust2)

# Test 1: Basic GO mapping creation
cat("Testing basic GO mapping creation...\n")
source("R/pathway_gsea.R")
go_mapping <- create_basic_go_mapping()
cat(sprintf("Created GO mapping with %d terms\n", nrow(go_mapping)))
cat(sprintf("Categories: %s\n", paste(unique(go_mapping$category), collapse = ", ")))

# Test 2: Gene sets preparation
cat("\nTesting gene sets preparation...\n")
gene_sets_bp <- prepare_gene_sets("GO", go_category = "BP")
cat(sprintf("BP gene sets: %d\n", length(gene_sets_bp)))

gene_sets_mf <- prepare_gene_sets("GO", go_category = "MF")
cat(sprintf("MF gene sets: %d\n", length(gene_sets_mf)))

gene_sets_cc <- prepare_gene_sets("GO", go_category = "CC")
cat(sprintf("CC gene sets: %d\n", length(gene_sets_cc)))

# Test 3: Sample GO terms
cat("\nSample GO terms:\n")
sample_terms <- head(names(gene_sets_bp), 3)
for (term in sample_terms) {
  kos <- gene_sets_bp[[term]]
  cat(sprintf("%s: %s\n", term, paste(head(kos, 5), collapse = ", ")))
}

# Test 4: Test with actual data (if available)
cat("\nTesting with sample data...\n")
tryCatch({
  data("ko_abundance", package = "ggpicrust2")
  data("metadata", package = "ggpicrust2")
  
  # Prepare small test dataset
  abundance_data <- as.data.frame(ko_abundance)
  
  # Check column names and sample names  
  cat(sprintf("KO abundance columns: %s\n", paste(head(colnames(abundance_data)), collapse = ", ")))
  cat(sprintf("Metadata rownames: %s\n", paste(head(rownames(metadata)), collapse = ", ")))
  
  # Prepare abundance properly
  rownames(abundance_data) <- abundance_data[, "#NAME"]
  abundance_data <- abundance_data[, -1]
  
  cat(sprintf("After processing - KO abundance columns: %s\n", paste(head(colnames(abundance_data)), collapse = ", ")))
  
  # Find common samples
  common_samples <- intersect(colnames(abundance_data), rownames(metadata))
  cat(sprintf("Common samples: %d\n", length(common_samples)))
  
  if (length(common_samples) >= 6) {
    # Subset to common samples
    abundance_subset <- abundance_data[1:50, common_samples]
    metadata_subset <- metadata[common_samples, , drop = FALSE]
    
    cat("Running GO GSEA test...\n")
    gsea_results <- pathway_gsea(
      abundance = abundance_subset,
      metadata = metadata_subset,
      group = "Environment",
      pathway_type = "GO",
      go_category = "BP",
      method = "fgsea",
      nperm = 50,
      min_size = 2,
      max_size = 500
    )
    
    cat(sprintf("GSEA results: %d pathways\n", nrow(gsea_results)))
    if (nrow(gsea_results) > 0) {
      cat("Sample results:\n")
      print(head(gsea_results[, c("pathway_id", "NES", "pvalue")], 3))
      
      # Test annotation
      annotated_results <- gsea_pathway_annotation(gsea_results, "GO")
      cat(sprintf("Annotation added pathway names: %d\n", 
                  sum(!is.na(annotated_results$pathway_name))))
    }
  } else {
    cat("Insufficient common samples for GSEA test\n")
  }
  
}, error = function(e) {
  cat(sprintf("Data test failed: %s\n", e$message))
})

cat("\nGO functionality test completed!\n")