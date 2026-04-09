# =============================================================================
# Taxa-Function Contribution Analysis
# =============================================================================
# Parse and aggregate PICRUSt2 per-sequence contribution data to identify
# which taxa drive differentially abundant pathways.

#' Read PICRUSt2 contribution file
#'
#' Parses the \code{pred_metagenome_contrib.tsv} file produced by PICRUSt2's
#' \code{--per_sequence_contrib} flag.
#'
#' @param file Path to the contribution file (.tsv, .txt, or .csv).
#' @param data A data.frame already loaded from the contribution file.
#'   If both \code{file} and \code{data} are provided, \code{data} is used.
#'
#' @return A data.frame with columns: \code{sample}, \code{function_id},
#'   \code{taxon}, \code{taxon_function_abun}, \code{norm_taxon_function_contrib},
#'   and any additional columns from the original file.
#'
#' @details
#' The contribution file records how much each ASV/OTU contributes to the
#' predicted abundance of each gene family in each sample. The
#' \code{norm_taxon_function_contrib} column gives the proportion of a
#' function's abundance attributable to each taxon.
#'
#' @examples
#' \donttest{
#' # From a data.frame
#' contrib_df <- data.frame(
#'   sample = rep(c("S1", "S2"), each = 4),
#'   `function` = rep(c("K00001", "K00002"), 4),
#'   taxon = rep(c("ASV1", "ASV2"), each = 2, times = 2),
#'   taxon_function_abun = runif(8),
#'   norm_taxon_function_contrib = runif(8),
#'   check.names = FALSE
#' )
#' result <- read_contrib_file(data = contrib_df)
#' head(result)
#' }
#'
#' @export
read_contrib_file <- function(file = NULL, data = NULL) {
  if (is.null(file) && is.null(data)) {
    stop("Please provide either a file path or a data.frame.")
  }
  if (!is.null(file) && !is.null(data)) {
    warning("Both file and data provided. Using data and ignoring file.")
  }

  if (!is.null(data)) {
    if (!is.data.frame(data)) stop("'data' must be a data.frame")
    contrib <- data
  } else {
    contrib <- read_abundance_file(file)
  }

  # Validate required columns
  function_col <- if ("function_id" %in% colnames(contrib)) {
    "function_id"
  } else if ("function" %in% colnames(contrib)) {
    "function"
  } else {
    NULL
  }

  required <- c("sample", "taxon",
                "taxon_function_abun", "norm_taxon_function_contrib")
  missing <- setdiff(required, colnames(contrib))
  if (is.null(function_col)) {
    missing <- c(missing, "function")
  }
  if (length(missing) > 0) {
    stop(sprintf(
      "Missing required columns: %s. Expected PICRUSt2 contrib format.",
      paste(missing, collapse = ", ")
    ))
  }

  if (identical(function_col, "function")) {
    # Rename 'function' -> 'function_id' (reserved word in R)
    colnames(contrib)[colnames(contrib) == "function"] <- "function_id"
  }

  # Strip ko: prefix
  contrib$function_id <- gsub("^ko:", "", contrib$function_id)

  contrib <- as.data.frame(contrib)
  contrib
}


#' Read PICRUSt2 stratified abundance file
#'
#' Parses the \code{pred_metagenome_strat.tsv} file produced by PICRUSt2's
#' \code{--strat_out} flag and converts the wide format to tidy long format.
#'
#' @param file Path to the stratified file (.tsv, .txt, or .csv).
#' @param data A data.frame already loaded from the stratified file.
#'   If both \code{file} and \code{data} are provided, \code{data} is used.
#'
#' @return A tidy data.frame with columns: \code{function_id}, \code{taxon},
#'   \code{sample}, \code{abundance}.
#'
#' @details
#' The stratified file has function IDs in the first column, sequence/taxon
#' IDs in the second column, and sample abundances in the remaining columns
#' (wide format). This function pivots to long format for downstream analysis.
#'
#' @examples
#' \donttest{
#' # From a data.frame
#' strat_df <- data.frame(
#'   `function` = c("K00001", "K00001", "K00002"),
#'   sequence = c("ASV1", "ASV2", "ASV1"),
#'   S1 = c(10, 5, 8),
#'   S2 = c(12, 3, 7),
#'   check.names = FALSE
#' )
#' result <- read_strat_file(data = strat_df)
#' head(result)
#' }
#'
#' @export
read_strat_file <- function(file = NULL, data = NULL) {
  if (is.null(file) && is.null(data)) {
    stop("Please provide either a file path or a data.frame.")
  }
  if (!is.null(file) && !is.null(data)) {
    warning("Both file and data provided. Using data and ignoring file.")
  }

  if (!is.null(data)) {
    if (!is.data.frame(data)) stop("'data' must be a data.frame")
    strat <- data
  } else {
    strat <- read_abundance_file(file)
  }

  # Validate minimum structure: function column + sequence column + >= 1 sample
  if (ncol(strat) < 3) {
    stop("Stratified file must have at least 3 columns (function, sequence, samples).")
  }

  function_col <- if ("function_id" %in% colnames(strat)) {
    "function_id"
  } else if ("function" %in% colnames(strat)) {
    "function"
  } else {
    NULL
  }
  taxon_col <- if ("taxon" %in% colnames(strat)) {
    "taxon"
  } else if ("sequence" %in% colnames(strat)) {
    "sequence"
  } else {
    NULL
  }

  if (is.null(function_col)) {
    stop("Expected a 'function' column in the stratified file.")
  }
  if (is.null(taxon_col)) {
    stop("Expected a 'sequence' column in the stratified file.")
  }

  if (identical(function_col, "function")) {
    colnames(strat)[colnames(strat) == "function"] <- "function_id"
  }
  if (identical(taxon_col, "sequence")) {
    colnames(strat)[colnames(strat) == "sequence"] <- "taxon"
  }

  # Pivot sample columns to long format
  sample_cols <- setdiff(colnames(strat), c("function_id", "taxon"))
  id_cols <- c("function_id", "taxon")
  strat <- tidyr::pivot_longer(
    strat,
    cols = -all_of(id_cols),
    names_to = "sample",
    values_to = "abundance"
  )

  # Strip ko: prefix
  strat$function_id <- gsub("^ko:", "", strat$function_id)

  as.data.frame(strat)
}


#' Aggregate taxa contributions for visualization
#'
#' Core aggregation function that bridges PICRUSt2 contribution data with
#' differential abundance analysis results. Optionally maps ASV/OTU IDs to
#' taxonomic names and filters to significant pathways.
#'
#' @param contrib_data A data.frame from \code{\link{read_contrib_file}} or
#'   \code{\link{read_strat_file}}.
#' @param taxonomy Optional data.frame mapping taxon IDs to taxonomy. Supports
#'   QIIME2 format (semicolon-delimited taxonomy strings) or DADA2 format
#'   (separate columns for each rank).
#' @param tax_level Character. Taxonomic rank for aggregation. One of
#'   \code{"Kingdom"}, \code{"Phylum"}, \code{"Class"}, \code{"Order"},
#'   \code{"Family"}, \code{"Genus"}, \code{"Species"}. Default \code{"Genus"}.
#' @param top_n Integer. Number of top taxa to keep; remaining are lumped as
#'   "Other". Default 10.
#' @param daa_results_df Optional data.frame from \code{\link{pathway_daa}},
#'   used to filter contributions to significant pathways.
#' @param pathway_ids Optional character vector of pathway IDs to filter.
#'   Alternative to \code{daa_results_df}.
#' @param p_threshold Numeric. Significance cutoff when using
#'   \code{daa_results_df}. Default 0.05.
#'
#' @return A tidy data.frame with columns: \code{sample}, \code{function_id},
#'   \code{taxon_label}, \code{contribution}.
#'
#' @details
#' When \code{daa_results_df} is provided, the function:
#' \enumerate{
#'   \item Extracts significant pathway IDs from the DAA results
#'   \item Maps pathway IDs to their constituent KO IDs using the internal
#'     ko_to_kegg reference
#'   \item Filters contribution data to only matching KO IDs
#' }
#'
#' Taxonomy can be provided in two formats:
#' \itemize{
#'   \item QIIME2: A column named \code{Taxon} or \code{taxonomy} containing
#'     semicolon-delimited strings (e.g., "k__Bacteria;p__Firmicutes;...")
#'   \item DADA2: Separate columns for each rank (\code{Kingdom}, \code{Phylum}, etc.)
#' }
#'
#' @examples
#' \donttest{
#' # Basic usage with synthetic data
#' contrib <- data.frame(
#'   sample = rep(c("S1", "S2"), each = 6),
#'   function_id = rep(c("K00001", "K00002", "K00003"), 4),
#'   taxon = rep(c("ASV1", "ASV2"), each = 3, times = 2),
#'   taxon_function_abun = runif(12),
#'   norm_taxon_function_contrib = runif(12)
#' )
#' agg <- aggregate_taxa_contributions(contrib, top_n = 2)
#' head(agg)
#' }
#'
#' @export
aggregate_taxa_contributions <- function(contrib_data,
                                         taxonomy = NULL,
                                         tax_level = "Genus",
                                         top_n = 10,
                                         daa_results_df = NULL,
                                         pathway_ids = NULL,
                                         p_threshold = 0.05) {
  # Validate input
  validate_dataframe(contrib_data, param_name = "contrib_data")

  if (!is.null(daa_results_df) && !is.null(pathway_ids)) {
    warning("Both daa_results_df and pathway_ids provided. Using daa_results_df.")
  }

  # Detect input format (contrib vs strat)
  if ("norm_taxon_function_contrib" %in% colnames(contrib_data)) {
    contrib_col <- "norm_taxon_function_contrib"
  } else if ("abundance" %in% colnames(contrib_data)) {
    contrib_col <- "abundance"
  } else {
    stop("contrib_data must contain 'norm_taxon_function_contrib' or 'abundance' column.")
  }

  # Ensure required columns exist
  required <- c("sample", "function_id", "taxon", contrib_col)
  missing <- setdiff(required, colnames(contrib_data))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
  }

  # Filter to significant pathways via DAA results
  if (!is.null(daa_results_df)) {
    sig_pathways <- daa_results_df$feature[daa_results_df$p_adjust < p_threshold]
    sig_pathways <- unique(sig_pathways[!is.na(sig_pathways)])
    if (length(sig_pathways) == 0) {
      stop("No significant pathways found at p_threshold = ", p_threshold)
    }
    # Map pathway IDs to KO IDs
    ko_to_kegg <- load_reference_data("ko_to_kegg")
    pathway_kos <- ko_to_kegg$ko_id[ko_to_kegg$pathway_id %in% sig_pathways]
    contrib_data <- contrib_data[contrib_data$function_id %in% pathway_kos, ]
  } else if (!is.null(pathway_ids)) {
    # Map pathway IDs to KO IDs
    ko_to_kegg <- load_reference_data("ko_to_kegg")
    pathway_kos <- ko_to_kegg$ko_id[ko_to_kegg$pathway_id %in% pathway_ids]
    contrib_data <- contrib_data[contrib_data$function_id %in% pathway_kos, ]
  }

  if (nrow(contrib_data) == 0) {
    stop("No data remaining after filtering. Check that function IDs match between contrib_data and pathway/DAA results.")
  }

  # Apply taxonomy labels
  if (!is.null(taxonomy)) {
    contrib_data <- apply_taxonomy(contrib_data, taxonomy, tax_level)
  } else {
    contrib_data$taxon_label <- contrib_data$taxon
  }

  # Aggregate by sample, function, taxon_label
  agg <- stats::aggregate(
    stats::as.formula(paste(contrib_col, "~ sample + function_id + taxon_label")),
    data = contrib_data,
    FUN = sum
  )
  colnames(agg)[colnames(agg) == contrib_col] <- "contribution"

  # Rank taxa and keep top_n
  taxa_means <- stats::aggregate(contribution ~ taxon_label, data = agg, FUN = mean)
  taxa_means <- taxa_means[order(-taxa_means$contribution), ]
  top_taxa <- utils::head(taxa_means$taxon_label, top_n)

  agg$taxon_label <- ifelse(agg$taxon_label %in% top_taxa, agg$taxon_label, "Other")

  # Re-aggregate after lumping
  agg <- stats::aggregate(
    contribution ~ sample + function_id + taxon_label,
    data = agg,
    FUN = sum
  )

  agg
}


# =============================================================================
# Internal helpers
# =============================================================================

#' Parse a QIIME2-style taxonomy string
#'
#' @param tax_string Character. Semicolon-delimited taxonomy string.
#' @param level Character. Rank to extract.
#' @return Character. The taxon name at the requested level, or NA.
#' @noRd
parse_taxonomy_string <- function(tax_string, level) {
  prefix_map <- c(
    Kingdom = "k__", Phylum = "p__", Class = "c__", Order = "o__",
    Family = "f__", Genus = "g__", Species = "s__"
  )
  prefix <- prefix_map[[level]]
  if (is.null(prefix)) return(NA_character_)

  parts <- strsplit(tax_string, ";\\s*")[[1]]
  match <- parts[startsWith(parts, prefix)]
  if (length(match) == 0) return(NA_character_)

  name <- sub(paste0("^", prefix), "", match[1])
  name <- trimws(name)
  if (nchar(name) == 0 || name == "__") return(NA_character_)
  name
}

#' Apply taxonomy mapping to contribution data
#'
#' @param contrib_data Data.frame with a \code{taxon} column.
#' @param taxonomy Data.frame with taxonomy information.
#' @param tax_level Character. Rank for labeling.
#' @return Data.frame with added \code{taxon_label} column.
#' @noRd
apply_taxonomy <- function(contrib_data, taxonomy, tax_level) {
  valid_levels <- c("Kingdom", "Phylum", "Class", "Order",
                    "Family", "Genus", "Species")
  if (!tax_level %in% valid_levels) {
    stop(sprintf(
      "Invalid tax_level '%s'. Must be one of: %s",
      tax_level, paste(valid_levels, collapse = ", ")
    ))
  }

  # Detect taxonomy format
  tax_cols <- colnames(taxonomy)
  id_col <- intersect(
    c("Feature.ID", "Feature ID", "feature_id", "ASV", "OTU", "#OTU ID"),
    tax_cols
  )

  # Also check first column as ID
  if (length(id_col) == 0) {
    id_col <- tax_cols[1]
  } else {
    id_col <- id_col[1]
  }

  # QIIME2 format: has a Taxon or taxonomy column with semicolon-delimited strings
  qiime2_col <- intersect(c("Taxon", "taxonomy", "Taxonomy"), tax_cols)

  if (length(qiime2_col) > 0) {
    # QIIME2 format
    qiime2_col <- qiime2_col[1]
    tax_labels <- vapply(
      taxonomy[[qiime2_col]],
      parse_taxonomy_string,
      character(1),
      level = tax_level
    )
    tax_map <- data.frame(
      taxon = taxonomy[[id_col]],
      taxon_label = tax_labels,
      stringsAsFactors = FALSE
    )
  } else if (tax_level %in% tax_cols) {
    # DADA2 format: separate columns per rank
    tax_map <- data.frame(
      taxon = taxonomy[[id_col]],
      taxon_label = taxonomy[[tax_level]],
      stringsAsFactors = FALSE
    )
  } else {
    warning(sprintf(
      "Could not find '%s' in taxonomy. Using raw taxon IDs.", tax_level
    ))
    contrib_data$taxon_label <- contrib_data$taxon
    return(contrib_data)
  }

  tax_map <- tax_map[!is.na(tax_map$taxon), , drop = FALSE]
  duplicate_taxa <- unique(tax_map$taxon[duplicated(tax_map$taxon)])
  if (length(duplicate_taxa) > 0) {
    warning(
      sprintf(
        "Duplicate taxon IDs found in taxonomy. Using the first match for %d taxa.",
        length(duplicate_taxa)
      ),
      call. = FALSE
    )
    tax_map <- tax_map[!duplicated(tax_map$taxon), , drop = FALSE]
  }

  match_idx <- match(contrib_data$taxon, tax_map$taxon)
  contrib_data$taxon_label <- tax_map$taxon_label[match_idx]
  contrib_data$taxon_label[
    is.na(contrib_data$taxon_label) | trimws(contrib_data$taxon_label) == ""
  ] <- "Unclassified"

  contrib_data
}
