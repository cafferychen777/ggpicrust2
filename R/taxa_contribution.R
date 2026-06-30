# =============================================================================
# Taxa-Function Contribution Analysis
# =============================================================================
# Parse and aggregate PICRUSt2 per-sequence contribution data to identify
# which taxa drive differentially abundant pathways.

#' Read PICRUSt2 contribution file
#'
#' Parses PICRUSt2 contribution files such as
#' \code{pred_metagenome_contrib.tsv}. It also accepts the long contribution
#' schema used by \code{path_abun_contrib.tsv}; for clarity, use
#' \code{\link{read_pathway_contrib_file}} when reading pathway-level
#' contribution output.
#'
#' @param file Path to the contribution file (.tsv, .txt, .csv, or gzipped variants).
#' @param data A data.frame already loaded from the contribution file.
#'   If both \code{file} and \code{data} are provided, \code{data} is used.
#' @param type Character. Contribution file level. One of \code{"auto"},
#'   \code{"gene_family"}, or \code{"pathway"}. Default \code{"auto"}.
#'
#' @return A data.frame with columns: \code{sample}, \code{function_id},
#'   \code{taxon}, contribution abundance columns from the original file, and
#'   \code{feature_level}.
#'
#' @details
#' The contribution file records how much each ASV/OTU contributes to the
#' predicted abundance of each gene family or pathway in each sample.
#' PICRUSt2 versions differ in whether they include
#' \code{norm_taxon_function_contrib}; when that column is absent downstream
#' aggregation can use \code{taxon_function_abun} or
#' \code{taxon_rel_function_abun}.
#' Contribution tables must describe one feature level at a time. Mixed
#' gene-family identifiers such as KOs/ECs and pathway identifiers are rejected,
#' because direct pathway matching and KEGG pathway-to-KO expansion have
#' different biological meanings.
#' Identifier columns (\code{sample}, \code{function_id}/\code{function}, and
#' \code{taxon}) must contain non-empty values without \code{NA}; otherwise
#' downstream aggregation would silently drop or mislabel contributions.
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
read_contrib_file <- function(file = NULL, data = NULL,
                              type = c("auto", "gene_family", "pathway")) {
  type <- match.arg(type)
  if (is.null(file) && is.null(data)) {
    stop("Please provide either a file path or a data.frame.")
  }
  if (!is.null(file) && !is.null(data)) {
    warning("Both file and data provided. Using data and ignoring file.")
  }

  if (!is.null(data)) {
    validate_dataframe(data, param_name = "data")
    contrib <- data
  } else {
    contrib <- read_abundance_file(file)
  }

  normalize_contrib_table(contrib, type = type)
}

#' Read PICRUSt2 pathway-level contribution file
#'
#' Parses \code{path_abun_contrib.tsv} output from PICRUSt2's pathway
#' pipeline and returns the same normalized contribution schema used by
#' \code{\link{aggregate_taxa_contributions}}.
#'
#' @param file Path to the contribution file (.tsv, .txt, .csv, or gzipped variants).
#' @param data A data.frame already loaded from the contribution file.
#'   If both \code{file} and \code{data} are provided, \code{data} is used.
#'
#' @return A normalized data.frame with pathway IDs in \code{function_id}.
#'
#' @examples
#' \donttest{
#' # From a data.frame
#' path_contrib_df <- data.frame(
#'   sample = rep(c("S1", "S2"), each = 2),
#'   `function` = rep(c("PWY-1234", "PWY-5678"), times = 2),
#'   taxon = c("ASV1", "ASV2", "ASV1", "ASV2"),
#'   taxon_function_abun = c(10, 5, 12, 4),
#'   check.names = FALSE
#' )
#' path_contrib <- read_pathway_contrib_file(data = path_contrib_df)
#' head(path_contrib)
#' }
#'
#' @export
read_pathway_contrib_file <- function(file = NULL, data = NULL) {
  read_contrib_file(file = file, data = data, type = "pathway")
}


#' Read PICRUSt2 stratified abundance file
#'
#' Parses the \code{pred_metagenome_strat.tsv} file produced by PICRUSt2's
#' \code{--strat_out} flag and converts the wide format to tidy long format.
#'
#' @param file Path to the stratified file (.tsv, .txt, .csv, or gzipped variants).
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
#' Function IDs, taxon IDs, and sample column names must be non-empty; sample
#' column names must also be unique.
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
    validate_dataframe(data, param_name = "data")
    strat <- data
  } else {
    strat <- read_abundance_file(file)
  }

  # Validate minimum structure: function column + sequence column + >= 1 sample
  if (ncol(strat) < 3) {
    stop("Stratified file must have at least 3 columns (function, sequence, samples).")
  }

  function_col <- resolve_function_id_column(strat)
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
  sample_cols <- colnames(strat)[!(colnames(strat) %in% c("function_id", "taxon"))]
  validate_contrib_key_columns(strat, c("function_id", "taxon"), "stratified file")
  validate_strat_sample_columns(sample_cols)
  id_cols <- c("function_id", "taxon")
  strat <- tidyr::pivot_longer(
    strat,
    cols = -all_of(id_cols),
    names_to = "sample",
    values_to = "abundance"
  )

  strat$function_id <- normalize_contrib_feature_ids(strat$function_id)
  validate_contrib_key_columns(strat, c("sample", "function_id", "taxon"),
                               "stratified file")

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
#' @param top_n Integer. Number of top taxa to keep by total contribution;
#'   remaining are lumped as "Other". Default 10.
#' @param daa_results_df Optional data.frame from \code{\link{pathway_daa}},
#'   used to filter contributions to significant pathways.
#' @param pathway_ids Optional character vector of pathway IDs to filter.
#'   Alternative to \code{daa_results_df}.
#' @param p_threshold Numeric. Significance cutoff when using
#'   \code{daa_results_df}. Default 0.05. Must be in the range (0, 1].
#' @param contribution_col Character. Column to aggregate. Use \code{"auto"}
#'   to select the first available column from
#'   \code{norm_taxon_function_contrib}, \code{taxon_function_abun},
#'   \code{taxon_rel_function_abun}, and \code{abundance}.
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
#' Contribution table identifier columns (\code{sample}, \code{function_id},
#' and \code{taxon}) and requested pathway/function IDs must be non-empty and
#' non-missing, because R's aggregation functions otherwise omit missing keys
#' without preserving contribution totals. Contribution tables must not mix
#' gene-family-level and pathway-level identifiers in the same aggregation.
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
                                         p_threshold = 0.05,
                                         contribution_col = "auto") {
  # Validate input
  validate_dataframe(contrib_data, param_name = "contrib_data")
  validate_count_parameter(top_n, "top_n")
  validate_probability_threshold(p_threshold, "p_threshold")

  if (!is.null(daa_results_df) && !is.null(pathway_ids)) {
    warning("Both daa_results_df and pathway_ids provided. Using daa_results_df.")
  }

  contrib_col <- choose_contribution_column(contrib_data, contribution_col)
  if (!is.numeric(contrib_data[[contrib_col]])) {
    contrib_data[[contrib_col]] <- normalize_contribution_values(
      contrib_data[[contrib_col]],
      contrib_col
    )
  }
  validate_contribution_values(contrib_data[[contrib_col]], contrib_col)

  # Ensure required columns exist
  required <- c("sample", "function_id", "taxon", contrib_col)
  missing <- setdiff(required, colnames(contrib_data))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")))
  }
  validate_contrib_key_columns(contrib_data,
                               c("sample", "function_id", "taxon"),
                               "contrib_data")
  validate_contrib_feature_level_column(contrib_data, "contrib_data")

  # KEGG pathway IDs are expanded to KOs only for KO-level contribution data.
  if (!is.null(daa_results_df)) {
    validate_dataframe(daa_results_df,
                       required_cols = c("feature", "p_adjust"),
                       param_name = "daa_results_df")
    validate_probability_values(daa_results_df$p_adjust,
                                "p_adjust",
                                "daa_results_df")
    daa_features <- validate_nonempty_character_column(
      daa_results_df$feature,
      "feature",
      "daa_results_df"
    )
    sig_features <- daa_features[daa_results_df$p_adjust < p_threshold]
    sig_features <- unique(sig_features[!is.na(sig_features)])
    if (length(sig_features) == 0) {
      stop("No significant features found at p_threshold = ", p_threshold)
    }
    contrib_data <- filter_contrib_by_features(contrib_data, sig_features)
  } else if (!is.null(pathway_ids)) {
    contrib_data <- filter_contrib_by_features(contrib_data, pathway_ids)
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
    build_two_sided_formula(
      contrib_col,
      c("sample", "function_id", "taxon_label")
    ),
    data = contrib_data,
    FUN = sum
  )
  colnames(agg)[colnames(agg) == contrib_col] <- "contribution"

  # Rank taxa by total contribution mass, not by mean over rows where the
  # taxon happens to appear. PICRUSt2 contribution tables are sparse in
  # practice (zero-contribution taxon/function/sample combinations are often
  # absent), so a mean over observed rows can over-rank a taxon with one
  # large contribution and under-rank a taxon that explains more total
  # functional abundance across samples/functions.
  taxa_totals <- stats::aggregate(contribution ~ taxon_label, data = agg, FUN = sum)
  taxa_totals <- taxa_totals[order(-taxa_totals$contribution), ]
  top_taxa <- utils::head(taxa_totals$taxon_label, top_n)

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

#' Filter contribution data by a vector of DAA features (KO or pathway IDs)
#'
#' Dispatches on ID shape: KO-like (`K\d{5}`) intersects directly with
#' `function_id`; pathway-like (`ko\d{5}`, `map\d{5}`, `ec\d{5}`) is first
#' expanded to its constituent KOs via the ko_to_kegg reference. Mixed
#' inputs are handled by combining both resolved sets.
#'
#' @noRd
filter_contrib_by_features <- function(contrib_data, features) {
  features <- unique(validate_nonempty_character_column(
    features,
    "features",
    "feature filter"
  ))
  function_ids <- as.character(contrib_data$function_id)

  direct_matches <- intersect(function_ids, features)
  matched_ids <- direct_matches

  unresolved <- setdiff(features, direct_matches)
  if (length(unresolved) > 0 && is_ko_level_contrib(contrib_data, function_ids)) {
    kegg_pathway_features <- unresolved[is_kegg_pathway_id(unresolved)]
    if (length(kegg_pathway_features) > 0) {
      ko_to_kegg <- load_reference_data("ko_to_kegg")
      pathway_kos <- ko_to_kegg$ko_id[
        ko_to_kegg$pathway_id %in% kegg_pathway_features
      ]
      matched_ids <- c(matched_ids, pathway_kos)
    }
  }

  matched_ids <- unique(matched_ids)

  if (length(matched_ids) == 0) {
    stop(
      "No contribution features matched the requested IDs. Direct matching is ",
      "used for pathway-level contribution data; KEGG pathway IDs are expanded ",
      "to KO IDs only when contribution data is KO-level. Requested examples: '",
      paste(head(features, 2), collapse = "', '"), "'. Contribution examples: '",
      paste(head(unique(function_ids), 2), collapse = "', '"), "'."
    )
  }

  filtered <- contrib_data[function_ids %in% matched_ids, ]
  if (nrow(filtered) == 0) {
    stop(
      "No data remaining after filtering. Check that function IDs match ",
      "between contrib_data (e.g. '",
      paste(head(unique(contrib_data$function_id), 2), collapse = "', '"),
      "') and the requested feature IDs (e.g. '",
      paste(head(features, 2), collapse = "', '"), "')."
    )
  }
  filtered
}

#' Contribution feature-level labels
#' @noRd
contribution_feature_levels <- function() {
  c(auto = "auto", gene_family = "gene_family", pathway = "pathway")
}

#' Normalize PICRUSt2 contribution table columns
#' @noRd
normalize_contrib_table <- function(contrib, type = contribution_feature_levels()) {
  type <- match.arg(type)
  contrib <- as.data.frame(contrib)

  function_col <- resolve_function_id_column(contrib)

  missing <- setdiff(c("sample", "taxon"), colnames(contrib))
  if (is.null(function_col)) {
    missing <- c(missing, "function")
  }

  if (!any(contribution_metric_candidates() %in% colnames(contrib))) {
    missing <- c(missing, paste(contribution_metric_candidates(), collapse = " or "))
  }

  if (length(missing) > 0) {
    stop(sprintf(
      "Missing required columns: %s. Expected PICRUSt2 contribution format.",
      paste(unique(missing), collapse = ", ")
    ))
  }

  if (identical(function_col, "function")) {
    colnames(contrib)[colnames(contrib) == "function"] <- "function_id"
  }

  contrib$function_id <- normalize_contrib_feature_ids(contrib$function_id)
  validate_contrib_key_columns(contrib,
                               c("sample", "function_id", "taxon"),
                               "contribution file")
  contrib$feature_level <- resolve_contrib_feature_level(
    contrib$function_id,
    requested_level = type,
    context = "contribution file"
  )

  contrib
}

#' Validate contribution table identifier columns
#' @noRd
validate_contrib_key_columns <- function(data, columns, context = "contrib_data") {
  for (column in columns) {
    validate_nonempty_character_column(data[[column]], column, context)
  }
  invisible(TRUE)
}

#' Validate sample columns in wide stratified contribution input
#' @noRd
validate_strat_sample_columns <- function(sample_cols) {
  sample_cols <- as.character(sample_cols)
  if (length(sample_cols) == 0) {
    stop("Stratified file must contain at least one sample column.",
         call. = FALSE)
  }
  bad <- is.na(sample_cols) | !nzchar(trimws(sample_cols))
  if (any(bad)) {
    stop("Stratified file sample column names must be non-empty.",
         call. = FALSE)
  }
  if (anyDuplicated(sample_cols)) {
    duplicated_samples <- unique(sample_cols[duplicated(sample_cols)])
    stop("Stratified file sample column names are duplicated: ",
         paste(utils::head(duplicated_samples, 5), collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Resolve a contribution function ID column
#' @noRd
resolve_function_id_column <- function(data) {
  if ("function_id" %in% colnames(data)) {
    "function_id"
  } else if ("function" %in% colnames(data)) {
    "function"
  } else {
    NULL
  }
}

#' Normalize contribution feature IDs without corrupting KEGG pathway IDs
#' @noRd
normalize_contrib_feature_ids <- function(ids) {
  ids <- as.character(ids)
  sub("^ko:(K\\d{5})$", "\\1", ids)
}

#' Detect whether contribution IDs represent gene families or pathways
#' @noRd
detect_contrib_feature_level <- function(ids) {
  resolve_contrib_feature_level(ids, requested_level = "auto")
}

#' Resolve and validate a contribution feature level
#' @noRd
resolve_contrib_feature_level <- function(ids,
                                          requested_level = "auto",
                                          context = "contrib_data") {
  requested_level <- match.arg(requested_level, contribution_feature_levels())
  ids <- as.character(ids)
  feature_levels <- contribution_feature_levels()
  pathway_like <- is_pathway_id(ids)
  gene_family_like <- is_gene_family_id(ids)

  if (requested_level == feature_levels[["auto"]]) {
    if (any(pathway_like) && any(!pathway_like)) {
      stop_mixed_contrib_feature_levels(ids, pathway_like, !pathway_like, context)
    }
    if (any(pathway_like)) {
      return(feature_levels[["pathway"]])
    }
    return(feature_levels[["gene_family"]])
  }

  if (requested_level == feature_levels[["gene_family"]]) {
    if (any(pathway_like)) {
      stop_feature_level_conflict(
        ids[pathway_like],
        requested_level,
        "pathway",
        context
      )
    }
    return(feature_levels[["gene_family"]])
  }

  if (any(gene_family_like)) {
    stop_feature_level_conflict(
      ids[gene_family_like],
      requested_level,
      "gene-family",
      context
    )
  }
  feature_levels[["pathway"]]
}

#' Validate an optional contribution feature_level column
#' @noRd
validate_contrib_feature_level_column <- function(contrib_data,
                                                  context = "contrib_data") {
  ids <- as.character(contrib_data$function_id)
  if (!"feature_level" %in% colnames(contrib_data)) {
    resolve_contrib_feature_level(ids, requested_level = "auto",
                                  context = context)
    return(invisible(TRUE))
  }

  levels <- as.character(contrib_data$feature_level)
  bad <- is.na(levels) | !nzchar(trimws(levels))
  if (any(bad)) {
    stop("Column 'feature_level' in ", context,
         " must contain non-empty values.", call. = FALSE)
  }

  allowed <- contribution_feature_levels()[c("gene_family", "pathway")]
  invalid <- unique(levels[!levels %in% allowed])
  if (length(invalid) > 0) {
    stop("Column 'feature_level' in ", context,
         " must contain only 'gene_family' or 'pathway'. Invalid values: ",
         paste(utils::head(invalid, 5), collapse = ", "), ".",
         call. = FALSE)
  }

  unique_levels <- unique(levels)
  if (length(unique_levels) > 1) {
    stop(context,
         " mixes contribution feature levels in column 'feature_level': ",
         paste(unique_levels, collapse = ", "),
         ". Aggregate one PICRUSt2 output level at a time.",
         call. = FALSE)
  }

  resolve_contrib_feature_level(ids,
                                requested_level = unique_levels[1],
                                context = context)
  invisible(TRUE)
}

#' Stop on mixed contribution feature levels
#' @noRd
stop_mixed_contrib_feature_levels <- function(ids, pathway_like,
                                              gene_family_like, context) {
  pathway_examples <- unique(ids[pathway_like])
  gene_family_examples <- unique(ids[gene_family_like])
  stop(
    context,
    " mixes pathway-level and gene-family-level identifiers in 'function_id'. ",
    "Aggregate one PICRUSt2 output level at a time. Pathway examples: '",
    paste(utils::head(pathway_examples, 3), collapse = "', '"),
    "'; gene-family examples: '",
    paste(utils::head(gene_family_examples, 3), collapse = "', '"), "'.",
    call. = FALSE
  )
}

#' Stop on explicit contribution feature-level conflicts
#' @noRd
stop_feature_level_conflict <- function(ids, requested_level, detected_label,
                                        context) {
  examples <- unique(as.character(ids))
  stop(
    context,
    " was declared as '",
    requested_level,
    "' but contains ",
    detected_label,
    "-level identifiers in 'function_id': '",
    paste(utils::head(examples, 3), collapse = "', '"),
    "'. Use the matching reader or aggregate one PICRUSt2 output level at a time.",
    call. = FALSE
  )
}

#' Check whether contribution data is KO-level
#' @noRd
is_ko_level_contrib <- function(contrib_data, function_ids) {
  if ("feature_level" %in% colnames(contrib_data)) {
    gene_family_level <- contribution_feature_levels()[["gene_family"]]
    validate_contrib_feature_level_column(contrib_data, "contrib_data")
    return(identical(unique(as.character(contrib_data$feature_level)),
                     gene_family_level))
  }
  identical(
    resolve_contrib_feature_level(function_ids, requested_level = "auto"),
    contribution_feature_levels()[["gene_family"]]
  )
}

#' Choose the contribution column used for aggregation
#' @noRd
choose_contribution_column <- function(contrib_data, contribution_col = "auto") {
  if (!is.character(contribution_col) || length(contribution_col) != 1) {
    stop("'contribution_col' must be a single column name or 'auto'.")
  }

  if (!identical(contribution_col, "auto")) {
    if (!contribution_col %in% colnames(contrib_data)) {
      stop("Column '", contribution_col, "' not found in contrib_data.")
    }
    return(contribution_col)
  }

  candidates <- contribution_metric_candidates()
  available <- candidates[candidates %in% colnames(contrib_data)]
  if (length(available) == 0) {
    stop(
      "contrib_data must contain one usable contribution column: ",
      paste(candidates, collapse = ", "), "."
    )
  }

  available[1]
}

#' Contribution metric columns in preference order
#' @noRd
contribution_metric_candidates <- function() {
  c(
    "norm_taxon_function_contrib",
    "taxon_function_abun",
    "taxon_rel_function_abun",
    "abundance"
  )
}

#' Normalize contribution metric values
#' @noRd
normalize_contribution_values <- function(values, contribution_col) {
  if (is.numeric(values)) {
    return(values)
  }

  converted <- suppressWarnings(as.numeric(values))
  if (any(is.na(converted) & !is.na(values))) {
    stop("Contribution column '", contribution_col, "' must be numeric.")
  }
  converted
}

#' Validate contribution metric values
#' @noRd
validate_contribution_values <- function(values, contribution_col) {
  if (!is.numeric(values)) {
    stop("Contribution column '", contribution_col, "' must be numeric.",
         call. = FALSE)
  }

  bad <- is.na(values) | !is.finite(values) | values < 0
  if (any(bad)) {
    stop(
      "Contribution column '", contribution_col,
      "' must contain only finite non-negative values. ",
      "Choose a different contribution_col if this PICRUSt2 output column ",
      "contains missing or invalid values.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Check KO identifiers
#' @noRd
is_ko_id <- function(ids) {
  grepl("^K\\d{5}$", as.character(ids))
}

#' Check EC identifiers
#' @noRd
is_ec_id <- function(ids) {
  grepl("^(EC[:_])?\\d+\\.(\\d+|-)\\.(\\d+|-)\\.(\\d+|-)$",
        as.character(ids))
}

#' Check gene-family identifiers supported by PICRUSt2 outputs
#' @noRd
is_gene_family_id <- function(ids) {
  ids <- as.character(ids)
  is_ko_id(ids) | is_ec_id(ids)
}

#' Check KEGG pathway identifiers
#' @noRd
is_kegg_pathway_id <- function(ids) {
  grepl("^(ko|map)\\d{5}$", as.character(ids))
}

#' Check pathway-like identifiers supported by PICRUSt2 outputs
#' @noRd
is_pathway_id <- function(ids) {
  ids <- as.character(ids)
  is_kegg_pathway_id(ids) |
    grepl("^(PWY|PWYG|PWY0|PWY\\d|[A-Z0-9]+-PWY)", ids) |
    grepl("-PWY$", ids)
}

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
