# Script to create example data for GSEA functionality

# Load required packages
library(dplyr)
library(tibble)

# Create example GSEA results
example_gsea_results <- data.frame(
  pathway_id = paste0("path:ko", sprintf("%05d", 1:20)),
  pathway_name = paste("Pathway", 1:20),
  size = sample(10:100, 20, replace = TRUE),
  ES = runif(20, -0.8, 0.8),
  NES = runif(20, -2, 2),
  pvalue = runif(20, 0, 0.1),
  p.adjust = runif(20, 0, 0.2),
  leading_edge = replicate(20, paste(paste0("K", sprintf("%05d", sample(1:1000, 5))), collapse = ";")),
  method = rep("fgsea", 20),
  stringsAsFactors = FALSE
)

# Add pathway classes
example_gsea_results$pathway_class <- sample(
  c("Metabolism", "Genetic Information Processing", "Environmental Information Processing", 
    "Cellular Processes", "Organismal Systems", "Human Diseases"),
  20, replace = TRUE
)

# Save the example data
usethis::use_data(example_gsea_results, overwrite = TRUE)

# Create example gene sets for KEGG pathways
ko_to_pathway_example <- list()
for (i in 1:20) {
  pathway_id <- paste0("path:ko", sprintf("%05d", i))
  ko_ids <- paste0("K", sprintf("%05d", sample(1:1000, sample(10:30, 1))))
  ko_to_pathway_example[[pathway_id]] <- ko_ids
}

# Save the example gene sets
usethis::use_data(ko_to_pathway_example, overwrite = TRUE)

# Print confirmation
cat("Example data created and saved.\n")
