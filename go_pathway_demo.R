# Demo: GO Pathway Analysis with ggpicrust2

library(ggpicrust2)
source("R/pathway_gsea.R")

cat("=== ggpicrust2 GO Pathway Analysis Demo ===\n\n")

# 1. Show available GO categories
cat("1. Available GO Categories:\n")
cat("   - BP: Biological Process\n")
cat("   - MF: Molecular Function\n")
cat("   - CC: Cellular Component\n\n")

# 2. Create sample data for demo
cat("2. Creating sample data...\n")
set.seed(42)

# Create mock KO abundance data
sample_kos <- c("K00134", "K01810", "K00927", "K01623", "K01803",  # Glycolysis
                "K01902", "K01903", "K00031", "K00164", "K00382",  # TCA cycle
                "K00059", "K00625", "K01895", "K07512", "K00626",  # Fatty acid
                "K02519", "K02543", "K02992", "K02946", "K02874")  # Translation

sample_names <- paste0("Sample_", 1:10)
group_data <- data.frame(
  Environment = rep(c("Condition_A", "Condition_B"), each = 5),
  row.names = sample_names
)

# Create abundance matrix
abundance_matrix <- matrix(
  rpois(length(sample_kos) * length(sample_names), lambda = 100),
  nrow = length(sample_kos),
  ncol = length(sample_names),
  dimnames = list(sample_kos, sample_names)
)

# Add some differential abundance
condition_a_samples <- sample_names[1:5]
condition_b_samples <- sample_names[6:10]

# Make first 5 KOs higher in condition A
abundance_matrix[1:5, condition_a_samples] <- abundance_matrix[1:5, condition_a_samples] * 2
# Make last 5 KOs higher in condition B  
abundance_matrix[16:20, condition_b_samples] <- abundance_matrix[16:20, condition_b_samples] * 2

cat(sprintf("Created abundance matrix: %d KOs × %d samples\n", nrow(abundance_matrix), ncol(abundance_matrix)))

# 3. Show GO gene sets
cat("\n3. GO Gene Sets by Category:\n")

bp_sets <- prepare_gene_sets("GO", go_category = "BP")
cat(sprintf("   Biological Process (BP): %d terms\n", length(bp_sets)))

mf_sets <- prepare_gene_sets("GO", go_category = "MF") 
cat(sprintf("   Molecular Function (MF): %d terms\n", length(mf_sets)))

cc_sets <- prepare_gene_sets("GO", go_category = "CC")
cat(sprintf("   Cellular Component (CC): %d terms\n", length(cc_sets)))

# Show sample terms
cat("\n   Sample BP terms:\n")
sample_bp_terms <- head(names(bp_sets), 5)
for (i in seq_along(sample_bp_terms)) {
  term <- sample_bp_terms[i]
  ko_count <- length(bp_sets[[term]])
  cat(sprintf("   %d. %s (%d KOs)\n", i, term, ko_count))
}

# 4. Run GO GSEA Analysis
cat("\n4. Running GO GSEA Analysis...\n")

tryCatch({
  # Run GSEA for each GO category
  for (category in c("BP", "MF", "CC")) {
    cat(sprintf("\n   Analyzing %s category...\n", category))
    
    gsea_results <- pathway_gsea(
      abundance = abundance_matrix,
      metadata = group_data,
      group = "Environment", 
      pathway_type = "GO",
      go_category = category,
      method = "fgsea",
      nperm = 100,
      min_size = 2,
      max_size = 50
    )
    
    cat(sprintf("   Found %d enriched pathways in %s\n", nrow(gsea_results), category))
    
    if (nrow(gsea_results) > 0) {
      # Show top results
      top_results <- head(gsea_results[order(gsea_results$pvalue), ], 3)
      cat(sprintf("   Top significant pathways:\n"))
      
      for (i in 1:nrow(top_results)) {
        row <- top_results[i, ]
        cat(sprintf("     %s (NES=%.3f, p=%.3f)\n", 
                   row$pathway_id, row$NES, row$pvalue))
      }
      
      # Test annotation
      annotated_results <- gsea_pathway_annotation(gsea_results, "GO")
      annotated_count <- sum(!is.na(annotated_results$pathway_name) & 
                           annotated_results$pathway_name != annotated_results$pathway_id)
      cat(sprintf("   Annotated pathway names: %d/%d\n", annotated_count, nrow(gsea_results)))
    }
  }
  
}, error = function(e) {
  cat(sprintf("GSEA analysis failed: %s\n", e$message))
})

# 5. Usage examples
cat("\n5. Usage Examples:\n")
cat("\n   # Basic GO analysis (BP category only)\n")
cat("   gsea_results <- pathway_gsea(\n")
cat("     abundance = your_ko_abundance,\n")
cat("     metadata = your_metadata,\n")
cat("     group = 'treatment',\n")
cat("     pathway_type = 'GO',\n")
cat("     go_category = 'BP'\n")
cat("   )\n")

cat("\n   # Analyze all GO categories\n")
cat("   all_categories <- c('BP', 'MF', 'CC')\n")
cat("   results_list <- list()\n")
cat("   for (cat in all_categories) {\n")
cat("     results_list[[cat]] <- pathway_gsea(\n")
cat("       abundance = your_abundance,\n") 
cat("       metadata = your_metadata,\n")
cat("       group = 'condition',\n")
cat("       pathway_type = 'GO',\n")
cat("       go_category = cat\n")
cat("     )\n")
cat("   }\n")

cat("\n   # Add pathway annotations\n")
cat("   annotated_results <- gsea_pathway_annotation(\n")
cat("     gsea_results = gsea_results,\n")
cat("     pathway_type = 'GO'\n")
cat("   )\n")

cat("\n=== Demo completed! ===\n")
cat("\nKey Features Implemented:\n")
cat("✓ GO term support for all three categories (BP, MF, CC)\n")
cat("✓ KO to GO term mapping\n") 
cat("✓ Category-specific analysis\n")
cat("✓ Integration with existing GSEA workflow\n")
cat("✓ Pathway name annotation\n")
cat("✓ Compatible with fgsea and clusterProfiler methods\n")

# Show current GO mapping stats
go_mapping <- create_basic_go_mapping()
bp_count <- sum(go_mapping$category == "BP")
mf_count <- sum(go_mapping$category == "MF") 
cc_count <- sum(go_mapping$category == "CC")

cat(sprintf("\nCurrent GO Database: %d total terms\n", nrow(go_mapping)))
cat(sprintf("  - Biological Process: %d terms\n", bp_count))
cat(sprintf("  - Molecular Function: %d terms\n", mf_count))
cat(sprintf("  - Cellular Component: %d terms\n", cc_count))