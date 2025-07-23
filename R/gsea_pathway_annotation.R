#' Annotate GSEA results with pathway information
#'
#' This function adds pathway annotations to GSEA results, including pathway names,
#' descriptions, and classifications.
#'
#' @param gsea_results A data frame containing GSEA results from the pathway_gsea function
#' @param pathway_type A character string specifying the pathway type: "KEGG", "MetaCyc", or "GO"
#'
#' @return A data frame with annotated GSEA results
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(ko_abundance)
#' data(metadata)
#' 
#' # Prepare abundance data
#' abundance_data <- as.data.frame(ko_abundance)
#' rownames(abundance_data) <- abundance_data[, "#NAME"]
#' abundance_data <- abundance_data[, -1]
#' 
#' # Run GSEA analysis
#' gsea_results <- pathway_gsea(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment",
#'   pathway_type = "KEGG",
#'   method = "fgsea"
#' )
#' 
#' # Annotate results
#' annotated_results <- gsea_pathway_annotation(
#'   gsea_results = gsea_results,
#'   pathway_type = "KEGG"
#' )
#' }
gsea_pathway_annotation <- function(gsea_results, 
                                   pathway_type = "KEGG") {
  
  # Input validation
  if (!is.data.frame(gsea_results)) {
    stop("'gsea_results' must be a data frame")
  }
  
  if (!pathway_type %in% c("KEGG", "MetaCyc", "GO")) {
    stop("pathway_type must be one of 'KEGG', 'MetaCyc', or 'GO'")
  }
  
  # Check if required columns exist
  if (!"pathway_id" %in% colnames(gsea_results)) {
    stop("GSEA results missing required column: pathway_id")
  }
  
  # Load appropriate reference data based on pathway_type
  if (pathway_type == "KEGG") {
    # Load KEGG pathway reference data
    if (!exists("kegg_reference")) {
      kegg_path <- system.file("extdata", "ko_pathway_reference.RData", package = "ggpicrust2")
      load(kegg_path, envir = environment())
    }
    
    # Convert to data frame
    kegg_reference <- as.data.frame(ko_pathway_reference)
    
    # Merge with GSEA results
    annotated_results <- merge(
      gsea_results,
      kegg_reference,
      by.x = "pathway_id",
      by.y = "pathway_id",
      all.x = TRUE
    )
    
    # If pathway_name is missing, use pathway_id
    if ("pathway_name" %in% colnames(annotated_results)) {
      annotated_results$pathway_name[is.na(annotated_results$pathway_name)] <- 
        annotated_results$pathway_id[is.na(annotated_results$pathway_name)]
    } else {
      annotated_results$pathway_name <- annotated_results$pathway_id
    }
    
  } else if (pathway_type == "MetaCyc") {
    # Load MetaCyc pathway reference data
    if (!exists("metacyc_reference")) {
      data("metacyc_reference", package = "ggpicrust2", envir = environment())
    }
    
    # Convert to data frame
    metacyc_reference <- as.data.frame(metacyc_reference)
    
    # Merge with GSEA results
    annotated_results <- merge(
      gsea_results,
      metacyc_reference,
      by.x = "pathway_id",
      by.y = "pathway",
      all.x = TRUE
    )
    
    # If pathway_name is missing, use pathway_id
    if ("pathway_name" %in% colnames(annotated_results)) {
      annotated_results$pathway_name[is.na(annotated_results$pathway_name)] <- 
        annotated_results$pathway_id[is.na(annotated_results$pathway_name)]
    } else {
      annotated_results$pathway_name <- annotated_results$pathway_id
    }
    
  } else if (pathway_type == "GO") {
    # For GO terms, we would need additional reference data
    # For now, we'll return the original results with a warning
    warning("GO pathway annotation not yet implemented")
    annotated_results <- gsea_results
  }
  
  return(annotated_results)
}
