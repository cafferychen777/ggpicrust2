#' Clean implementation of prepare_gene_sets following Linus principles
prepare_gene_sets <- function(pathway_type = "KEGG", organism = "ko", go_category = "BP") {
  
  # Load and process pathway data based on type - eliminates special cases
  gene_sets <- switch(pathway_type,
    "KEGG" = load_kegg_gene_sets(organism),
    "MetaCyc" = load_metacyc_gene_sets(),
    "GO" = load_go_gene_sets(go_category),
    stop("pathway_type must be one of 'KEGG', 'MetaCyc', or 'GO'")
  )
  
  # Universal validation - no special cases needed
  if (validate_pathway_data(gene_sets, pathway_type)) {
    message(sprintf("%s gene sets prepared successfully", pathway_type))
  }
  
  return(gene_sets)
}