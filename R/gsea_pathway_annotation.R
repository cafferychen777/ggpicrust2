#' Annotate GSEA results with pathway information
#'
#' This function adds pathway annotations to GSEA results, including pathway names,
#' descriptions, and classifications.
#'
#' @param gsea_results A data frame containing GSEA results from the pathway_gsea function
#' @param pathway_type A character string specifying the pathway type: "KEGG", "MetaCyc", or "GO"
#'
#' @return A data frame with annotated GSEA results
#'
#' @details
#' The \code{pathway_id} column must contain non-empty values without
#' \code{NA}. Unknown but non-empty pathway IDs are retained as their own
#' display names when no reference annotation is found.
#'
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
  if (!is.character(pathway_type) || length(pathway_type) != 1 ||
      is.na(pathway_type) || !pathway_type %in% valid_types) {
    stop(sprintf("pathway_type must be one of: %s", paste(valid_types, collapse = ", ")))
  }

  if (!"pathway_id" %in% colnames(gsea_results)) {
    stop("GSEA results missing required column: pathway_id")
  }
  gsea_results$pathway_id <- validate_nonempty_character_column(
    gsea_results$pathway_id,
    "pathway_id",
    "gsea_results"
  )

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
  annotate_gsea_from_lookup(
    gsea_results,
    reference_ids = kegg_ref$pathway,
    reference_names = kegg_ref$pathway_name,
    reference_name = "KEGG reference"
  )
}

#' Annotate GSEA results with MetaCyc pathway information
#' @param gsea_results GSEA results data frame
#' @return Annotated results
#' @noRd
annotate_metacyc_gsea <- function(gsea_results) {
  # Load MetaCyc reference using unified loader
  metacyc_ref <- load_reference_data("MetaCyc")
  annotate_gsea_from_lookup(
    gsea_results,
    reference_ids = metacyc_ref$id,
    reference_names = metacyc_ref$description,
    reference_name = "MetaCyc reference"
  )
}

#' Annotate GSEA results with GO term information
#' @param gsea_results GSEA results data frame
#' @return Annotated results
#' @noRd
annotate_go_gsea <- function(gsea_results) {
  # Load GO reference using unified loader
  go_ref <- load_reference_data("ko_to_go")
  annotate_gsea_from_lookup(
    gsea_results,
    reference_ids = go_ref$go_id,
    reference_names = go_ref$go_name,
    reference_name = "GO reference"
  )
}

#' Annotate GSEA rows through a strict many-to-one lookup
#'
#' @noRd
annotate_gsea_from_lookup <- function(gsea_results,
                                      reference_ids,
                                      reference_names,
                                      reference_name) {
  reference_ids <- as.character(reference_ids)
  reference_names <- as.character(reference_names)
  if (length(reference_ids) != length(reference_names)) {
    stop(reference_name, " has mismatched ID and name lengths.",
         call. = FALSE)
  }

  valid_ids <- !is.na(reference_ids) & nzchar(trimws(reference_ids))
  reference_ids <- reference_ids[valid_ids]
  reference_names <- reference_names[valid_ids]
  if (anyDuplicated(reference_ids)) {
    duplicated_ids <- unique(reference_ids[duplicated(reference_ids)])
    stop(
      reference_name,
      " contains duplicated pathway IDs: ",
      paste(utils::head(duplicated_ids, 5), collapse = ", "),
      ". Each pathway ID must map to exactly one annotation row.",
      call. = FALSE
    )
  }

  lookup_index <- match(gsea_results$pathway_id, reference_ids)
  pathway_names <- reference_names[lookup_index]
  missing_names <- is.na(pathway_names) | !nzchar(trimws(pathway_names))
  pathway_names[missing_names] <- gsea_results$pathway_id[missing_names]

  gsea_results$pathway_name <- pathway_names
  rownames(gsea_results) <- NULL
  gsea_results
}
