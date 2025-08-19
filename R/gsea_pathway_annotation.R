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
      # Load directly from extdata file since data() doesn't work for this dataset
      kegg_ref_path <- system.file("extdata", "kegg_reference.RData", package = "ggpicrust2")
      if (file.exists(kegg_ref_path)) {
        # Load the file into current environment
        load(kegg_ref_path)
      } else {
        stop("kegg_reference data file not found")
      }
    }
    
    # Convert to data frame
    kegg_reference <- as.data.frame(kegg_reference)
    
    # Merge with GSEA results
    annotated_results <- merge(
      gsea_results,
      kegg_reference,
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
    
  } else if (pathway_type == "MetaCyc") {
    # Load MetaCyc pathway reference data
    if (!exists("MetaCyc_reference")) {
      # Load directly from extdata file
      metacyc_ref_path <- system.file("extdata", "MetaCyc_reference.RData", package = "ggpicrust2")
      if (file.exists(metacyc_ref_path)) {
        # Load the file into current environment
        load(metacyc_ref_path)
      } else {
        stop("MetaCyc_reference data file not found")
      }
    }
    
    # Convert to data frame
    MetaCyc_reference <- as.data.frame(MetaCyc_reference)
    
    # Rename columns for consistency
    colnames(MetaCyc_reference) <- c("pathway", "description")
    
    # Merge with GSEA results
    annotated_results <- merge(
      gsea_results,
      MetaCyc_reference,
      by.x = "pathway_id",
      by.y = "pathway",
      all.x = TRUE
    )
    
    # Update pathway names using the description from reference data
    if ("description" %in% colnames(annotated_results)) {
      # Use description as pathway_name, fall back to pathway_id if missing
      annotated_results$pathway_name <- ifelse(
        is.na(annotated_results$description) | annotated_results$description == "",
        annotated_results$pathway_id,
        annotated_results$description
      )
      # Remove the description column to avoid confusion
      annotated_results$description <- NULL
    } else {
      annotated_results$pathway_name <- annotated_results$pathway_id
    }
    
  } else if (pathway_type == "GO") {
    # Load GO reference data
    tryCatch({
      data("ko_to_go_reference", package = "ggpicrust2", envir = environment())
    }, error = function(e) {
      # Create basic GO mapping if reference data doesn't exist
      message("Creating basic GO mapping for annotation")
    })
    
    # Always create mapping if it doesn't exist
    if (!exists("ko_to_go_reference", envir = environment())) {
      # Source the pathway_gsea.R to get the create_basic_go_mapping function
      tryCatch({
        source_file <- system.file("R", "pathway_gsea.R", package = "ggpicrust2")
        if (file.exists(source_file)) {
          source(source_file, local = TRUE)
        } else {
          # Try relative path
          source("R/pathway_gsea.R", local = TRUE)
        }
        ko_to_go_reference <- create_basic_go_mapping()
      }, error = function(e) {
        # Fallback: create a minimal mapping inline
        ko_to_go_reference <- data.frame(
          go_id = c("GO:0006096", "GO:0006099", "GO:0006631", "GO:0006520"),
          go_name = c("Glycolytic process", "Tricarboxylic acid cycle", "Fatty acid metabolic process", "Cellular amino acid metabolic process"),
          category = c("BP", "BP", "BP", "BP"),
          ko_members = c("K00134;K01810", "K01902;K01903", "K00059;K00625", "K01915;K00928"),
          stringsAsFactors = FALSE
        )
      })
    }
    
    # Convert to data frame
    go_reference <- as.data.frame(ko_to_go_reference)
    
    # Create lookup data frame for merging
    go_lookup <- data.frame(
      pathway_id = go_reference$go_id,
      pathway_name = go_reference$go_name,
      category = if("category" %in% colnames(go_reference)) go_reference$category else "BP",
      stringsAsFactors = FALSE
    )
    
    # Merge with GSEA results
    annotated_results <- merge(
      gsea_results,
      go_lookup,
      by = "pathway_id",
      all.x = TRUE
    )
    
    # If pathway_name is missing, use pathway_id
    if ("pathway_name.y" %in% colnames(annotated_results)) {
      annotated_results$pathway_name <- ifelse(
        is.na(annotated_results$pathway_name.y),
        annotated_results$pathway_id,
        annotated_results$pathway_name.y
      )
      # Remove duplicate columns
      annotated_results$pathway_name.x <- NULL
      annotated_results$pathway_name.y <- NULL
    } else if (!"pathway_name" %in% colnames(annotated_results)) {
      annotated_results$pathway_name <- annotated_results$pathway_id
    }
  }
  
  return(annotated_results)
}
