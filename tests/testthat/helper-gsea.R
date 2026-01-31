# Shared helper: create mock GSEA results for testing
# pathway_name defaults to pathway_id (unannotated); callers can override
create_gsea_test_results <- function(n_pathways = 10) {
  set.seed(123)
  data.frame(
    pathway_id = paste0("path:ko", sprintf("%05d", 1:n_pathways)),
    pathway_name = paste0("path:ko", sprintf("%05d", 1:n_pathways)),
    size = sample(10:100, n_pathways, replace = TRUE),
    ES = runif(n_pathways, -0.8, 0.8),
    NES = runif(n_pathways, -2, 2),
    pvalue = runif(n_pathways, 0, 0.1),
    p.adjust = runif(n_pathways, 0, 0.2),
    leading_edge = replicate(n_pathways, paste(paste0("K", sprintf("%05d", sample(1:1000, 5))), collapse = ";")),
    method = rep("fgsea", n_pathways),
    stringsAsFactors = FALSE
  )
}
