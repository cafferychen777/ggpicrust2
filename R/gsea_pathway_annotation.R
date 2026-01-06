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
#' # Run GSEA analysis (using camera method - recommended)
#' gsea_results <- pathway_gsea(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment",
#'   pathway_type = "KEGG",
#'   method = "camera"
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

  valid_types <- c("KEGG", "MetaCyc", "GO")
  if (length(pathway_type) != 1 || !pathway_type %in% valid_types) {
    stop(sprintf("pathway_type must be one of: %s", paste(valid_types, collapse = ", ")))
  }

  if (!"pathway_id" %in% colnames(gsea_results)) {
    stop("GSEA results missing required column: pathway_id")
  }

  # Annotate based on pathway type
  if (pathway_type == "KEGG") {
    annotated_results <- annotate_kegg_gsea(gsea_results)
  } else if (pathway_type == "MetaCyc") {
    annotated_results <- annotate_metacyc_gsea(gsea_results)
  } else if (pathway_type == "GO") {
    annotated_results <- annotate_go_gsea(gsea_results)
  }

  return(annotated_results)
}

#' Annotate GSEA results with KEGG pathway information
#' @param gsea_results GSEA results data frame
#' @return Annotated results
#' @noRd
annotate_kegg_gsea <- function(gsea_results) {
  # Load KEGG pathway reference using unified loader
  kegg_ref <- load_reference_data("KEGG")

  # Remove pathway_name from gsea_results if present (will be replaced by reference)
  if ("pathway_name" %in% colnames(gsea_results)) {
    gsea_results$pathway_name <- NULL
  }

  # Merge with GSEA results
  annotated_results <- merge(
    gsea_results,
    kegg_ref,
    by.x = "pathway_id",
    by.y = "pathway",
    all.x = TRUE
  )

  # Fill missing pathway names with pathway_id
  if ("pathway_name" %in% colnames(annotated_results)) {
    annotated_results$pathway_name[is.na(annotated_results$pathway_name)] <-
      annotated_results$pathway_id[is.na(annotated_results$pathway_name)]
  } else {
    annotated_results$pathway_name <- annotated_results$pathway_id
  }

  annotated_results
}

#' Annotate GSEA results with MetaCyc pathway information
#' @param gsea_results GSEA results data frame
#' @return Annotated results
#' @noRd
annotate_metacyc_gsea <- function(gsea_results) {
  # Load MetaCyc reference using unified loader
  metacyc_ref <- load_reference_data("MetaCyc")

  # Remove pathway_name from gsea_results if present (will be replaced by reference)
  if ("pathway_name" %in% colnames(gsea_results)) {
    gsea_results$pathway_name <- NULL
  }

  # Merge with GSEA results
  annotated_results <- merge(
    gsea_results,
    metacyc_ref,
    by.x = "pathway_id",
    by.y = "id",
    all.x = TRUE
  )

  # Use description as pathway_name
  annotated_results$pathway_name <- ifelse(
    is.na(annotated_results$description) | annotated_results$description == "",
    annotated_results$pathway_id,
    annotated_results$description
  )
  annotated_results$description <- NULL

  annotated_results
}

#' Annotate GSEA results with GO term information
#' @param gsea_results GSEA results data frame
#' @return Annotated results
#' @noRd
annotate_go_gsea <- function(gsea_results) {
  # Load GO reference using unified loader
  go_ref <- load_reference_data("ko_to_go")

  # Remove pathway_name from gsea_results if present (will be replaced by reference)
  if ("pathway_name" %in% colnames(gsea_results)) {
    gsea_results$pathway_name <- NULL
  }

  # Create lookup data frame for merging (unique GO IDs with names)
  go_lookup <- unique(data.frame(
    pathway_id = go_ref$go_id,
    pathway_name = go_ref$go_name,
    stringsAsFactors = FALSE
  ))

  # Merge with GSEA results
  annotated_results <- merge(
    gsea_results,
    go_lookup,
    by = "pathway_id",
    all.x = TRUE
  )

  # Fill missing pathway names with pathway_id
  if ("pathway_name" %in% colnames(annotated_results)) {
    annotated_results$pathway_name[is.na(annotated_results$pathway_name)] <-
      annotated_results$pathway_id[is.na(annotated_results$pathway_name)]
  } else {
    annotated_results$pathway_name <- annotated_results$pathway_id
  }

  annotated_results
}
