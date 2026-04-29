# Note: read_abundance_file() is defined in data_utils.R
# Note: load_reference_data() is defined in data_utils.R

#' Cache manager for KEGG annotations
#' @noRd
kegg_cache <- new.env(parent = emptyenv())

#' Get KEGG annotation with caching
#' @param ko_id KO identifier
#' @param organism Organism code (e.g., 'hsa' for human, 'eco' for E. coli). If NULL, the generic KO entry is retrieved.
#' @return KEGG annotation
#' @noRd
get_kegg_with_cache <- function(ko_id, organism = NULL) {
  # Cache key includes organism so organism-specific pathway IDs do not
  # collide with the generic entries.
  cache_key <- if (!is.null(organism)) {
    paste0(organism, ":", ko_id)
  } else {
    ko_id
  }

  if (!exists(cache_key, envir = kegg_cache)) {
    # Throttle to avoid hammering the KEGG API
    Sys.sleep(0.1)

    # Always fetch the generic KO entry. The previous organism-specific
    # branch picked the *first* gene returned by keggLink() as a
    # "representative" and fetched that gene's record instead, but gene
    # order from keggLink is not semantically meaningful: isozymes or
    # paralogs can belong to different pathways, so the same KO could
    # yield different annotations on different KEGG builds. The generic
    # KO entry is the authoritative KO-level annotation; organism-specific
    # pathway maps share the same numeric IDs as the generic `ko`/`map`
    # pathways and only differ by prefix (e.g. ko00010 <-> hsa00010), so
    # we rewrite the prefix below rather than drilling into individual
    # genes.
    result <- with_retry({
      KEGGREST::keggGet(ko_id)
    })

    if (inherits(result, "kegg_error")) {
      result$ko_id <- ko_id
      return(result)
    }

    if (!is.null(result) && !is.null(organism)) {
      result <- rewrite_kegg_pathway_organism(result, organism)
    }

    if (!is.null(result)) {
      assign(cache_key, result, envir = kegg_cache)
    }
  }

  if (exists(cache_key, envir = kegg_cache)) {
    return(get(cache_key, envir = kegg_cache, inherits = FALSE))
  }

  NULL
}

#' Rewrite KO pathway IDs to organism-specific prefix
#'
#' KO entries reference pathways by their generic `ko`-prefixed IDs
#' (e.g. `ko00010`). KEGG mirrors every generic pathway under the
#' organism code using the same numeric suffix (e.g. `hsa00010`).
#' Substituting the prefix in `PATHWAY` and `PATHWAY_MAP` is equivalent
#' to the organism-specific pathway projection while keeping the KO-level
#' name/description intact, and avoids the fragile first-gene heuristic.
#' @noRd
rewrite_kegg_pathway_organism <- function(result, organism) {
  if (is.null(result) || length(result) == 0) {
    return(result)
  }
  rewrite <- function(ids) {
    if (is.null(ids)) return(ids)
    # Replace leading `ko` (KO pathway prefix) with the organism code.
    sub("^ko", organism, ids)
  }
  for (i in seq_along(result)) {
    entry <- result[[i]]
    if (!is.list(entry)) next
    for (fld in c("PATHWAY", "PATHWAY_MAP")) {
      if (!is.null(entry[[fld]]) && length(entry[[fld]]) > 0) {
        names(entry[[fld]]) <- rewrite(names(entry[[fld]]))
      }
    }
    result[[i]] <- entry
  }
  result
}

#' Retry mechanism for KEGG queries
#' @param expr Expression or function to retry
#' @param max_attempts Maximum number of retry attempts
#' @return Result or error
#' @noRd
with_retry <- function(expr, max_attempts = getOption("ggpicrust2.max_retries", 3)) {
  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch({
      if (is.function(expr)) expr() else if (is.expression(expr)) eval(expr) else expr
    }, error = function(e) {
      structure(list(
        status = if (grepl("HTTP 404", e$message)) "not_found" else "failed",
        message = e$message,
        ko_id = NULL
      ), class = "kegg_error")
    })

    # Success
    if (!inherits(result, "kegg_error")) return(result)

    # 404 is permanent — retrying won't help
    if (result$status == "not_found") return(result)

    # Transient failure: log and back off before next attempt
    if (attempt < max_attempts) {
      log_message(sprintf("Attempt %d/%d failed: %s, retrying...",
                          attempt, max_attempts, result$message), "WARN")
      Sys.sleep(min(2^attempt, 30))
    }
  }

  # All retries exhausted: return last error (always a kegg_error, never NULL)
  result
}

#' Enhanced logging system
#' @noRd
log_message <- function(msg, level = "INFO") {
  if (getOption("ggpicrust2.verbose", default = TRUE)) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(sprintf("[%s] %s: %s", timestamp, level, msg))
  }
}

#' Create enhanced progress bar
#' @noRd
create_progress_bar <- function(total, format = NULL) {
  if (is.null(format)) {
    format <- paste0(
      "  [:bar] :percent | Elapsed: :elapsed | ETA: :eta\n",
      "  :current/:total (:rate/sec) | :what"
    )
  }
  
  progress::progress_bar$new(
    format = format,
    total = total,
    clear = FALSE,
    width = 80
  )
}

#' Process KEGG annotations with improved error handling and caching
#' @param df Data frame with features to annotate
#' @param organism KEGG organism code (e.g., 'hsa' for human)
#' @return Annotated data frame
#' @noRd
process_kegg_annotations <- function(df, organism = NULL, p_adjust_threshold = 0.05) {
  if (nrow(df) == 0) {
    stop("Empty data frame provided for KEGG annotation")
  }

  # Filter for significant pathways
  filtered_df <- df[which(df$p_adjust < p_adjust_threshold), ]

  # Handle case when no significant pathways are found
  if (nrow(filtered_df) == 0) {
    min_p <- min(df$p_adjust, na.rm = TRUE)

    warning(
      "\n================================================================================\n",
      "NO SIGNIFICANT PATHWAYS FOUND\n",
      "================================================================================\n\n",
      sprintf("pathway_annotation() filters results by p_adjust < %g threshold.\n", p_adjust_threshold),
      "Your data contains no pathways meeting this criterion.\n\n",
      "Statistics:\n",
      "  Total features in input:       ", nrow(df), "\n",
      sprintf("  Significant features (p<%g): 0\n", p_adjust_threshold),
      "  Minimum p_adjust value:        ", sprintf("%.6f", min_p), "\n\n",
      "Possible reasons:\n",
      "  1. Insufficient statistical power (small sample size)\n",
      "  2. Low biological effect size (small differences between groups)\n",
      "  3. High biological variability in your samples\n",
      "  4. DAA method may not be suitable for your data\n\n",
      "Recommendations:\n",
      "  1. Check input data quality and normalization\n",
      "  2. Verify adequate sample size in each group\n",
      "  3. Try alternative DAA methods (LinDA, DESeq2, edgeR)\n",
      "  4. For annotation only, use ko_to_kegg = FALSE (uses local reference)\n\n",
      "Returning data with NA annotation columns.\n",
      "================================================================================\n",
      call. = FALSE
    )

    # Return original data with empty annotation columns for compatibility
    new_cols <- c("pathway_name", "pathway_description", "pathway_class", "pathway_map")
    df[new_cols] <- NA_character_

    log_message("No significant pathways found. Returning data with empty annotation columns.", "WARN")
    return(df)
  }
  
  # Initialize new columns for significant pathways
  new_cols <- c("pathway_name", "pathway_description", "pathway_class", "pathway_map")
  filtered_df[new_cols] <- NA_character_
  
  # Create progress bar
  total_features <- nrow(filtered_df)
  pb <- create_progress_bar(total_features)
  log_message("Starting KEGG annotation process")
  
  # Track statistics
  success_count <- 0
  not_found_count <- 0
  error_count <- 0
  not_found_ids <- character(0)
  error_ids <- character(0)
  
  # Process each feature
  for (i in seq_len(nrow(filtered_df))) {
    ko_id <- filtered_df$feature[i]
    
    # Update progress bar
    pb$tick(tokens = list(what = sprintf("Processing %s", ko_id)))
    
    # Get KEGG annotation
    entry <- tryCatch({
      result <- get_kegg_with_cache(ko_id, organism)
      
      # Handle error objects
      if (inherits(result, "kegg_error")) {
        if (result$status == "not_found") {
          not_found_count <- not_found_count + 1
          not_found_ids <- c(not_found_ids, ko_id)
          log_message(sprintf("KO ID %s not found in KEGG database (HTTP 404)", ko_id), "WARN")
        } else {
          error_count <- error_count + 1
          error_ids <- c(error_ids, ko_id)
          log_message(sprintf("Error querying KEGG for %s: %s", ko_id, result$message), "ERROR")
        }
        next
      }
      
      if (is.null(result)) {
        log_message(sprintf("No KEGG data found for %s", ko_id), "WARN")
        not_found_count <- not_found_count + 1
        not_found_ids <- c(not_found_ids, ko_id)
        next
      }
      
      success_count <- success_count + 1
      result
    }, error = function(e) {
      log_message(sprintf("Error processing %s: %s", ko_id, e$message), "ERROR")
      error_count <- error_count + 1
      error_ids <- c(error_ids, ko_id)
      NULL
    })
    
    # Safely extract data, check if fields exist
    if (!is.null(entry) && length(entry) > 0) {
      filtered_df$pathway_name[i] <- safe_extract(entry[[1]], "NAME", 1)

      # Prefer PATHWAY (KO entries); fall back to DESCRIPTION (pathway entries).
      if ("PATHWAY" %in% names(entry[[1]]) &&
          !is.null(entry[[1]][["PATHWAY"]]) &&
          length(entry[[1]][["PATHWAY"]]) > 0) {
        pathway_names <- as.character(entry[[1]][["PATHWAY"]])
        filtered_df$pathway_description[i] <- paste(pathway_names, collapse = "; ")
      } else if ("DESCRIPTION" %in% names(entry[[1]]) &&
                 !is.null(entry[[1]][["DESCRIPTION"]]) &&
                 length(entry[[1]][["DESCRIPTION"]]) > 0) {
        filtered_df$pathway_description[i] <- paste(as.character(entry[[1]][["DESCRIPTION"]]), collapse = "; ")
      } else {
        filtered_df$pathway_description[i] <- NA_character_
      }

      # Use CLASS field for pathway_class (contains functional classification)
      # KEGG returns CLASS, not BRITE for pathway entries
      if("CLASS" %in% names(entry[[1]]) && !is.null(entry[[1]][["CLASS"]]) && length(entry[[1]][["CLASS"]]) > 0) {
        # Extract CLASS field
        pathway_class <- as.character(entry[[1]][["CLASS"]])
        filtered_df$pathway_class[i] <- paste(pathway_class, collapse = "; ")
      } else if("BRITE" %in% names(entry[[1]]) && !is.null(entry[[1]][["BRITE"]]) && length(entry[[1]][["BRITE"]]) > 0) {
        # Fallback to BRITE if CLASS is not available
        brite_classes <- as.character(entry[[1]][["BRITE"]])
        filtered_df$pathway_class[i] <- paste(head(brite_classes, 3), collapse = "; ")
      } else {
        filtered_df$pathway_class[i] <- NA_character_
      }

      # Use PATHWAY_MAP field for pathway_map
      if("PATHWAY_MAP" %in% names(entry[[1]]) && !is.null(entry[[1]][["PATHWAY_MAP"]]) && length(entry[[1]][["PATHWAY_MAP"]]) > 0) {
        # Extract pathway map IDs (names of the PATHWAY_MAP vector)
        pathway_maps <- names(entry[[1]][["PATHWAY_MAP"]])
        filtered_df$pathway_map[i] <- paste(pathway_maps, collapse = "; ")
      } else if("PATHWAY" %in% names(entry[[1]]) && !is.null(entry[[1]][["PATHWAY"]]) && length(entry[[1]][["PATHWAY"]]) > 0) {
        # Fallback to PATHWAY field if PATHWAY_MAP is not available
        pathway_maps <- names(entry[[1]][["PATHWAY"]])
        filtered_df$pathway_map[i] <- paste(pathway_maps, collapse = "; ")
      } else {
        filtered_df$pathway_map[i] <- NA_character_
      }
    }
  }
  
  log_message("KEGG annotation process completed")
  
  # Provide summary information after processing all KO IDs
  log_message(sprintf("KEGG annotation complete: %d successful, %d not found, %d errors", 
                      success_count, not_found_count, error_count))
  
  if (not_found_count > 0) {
    if (length(not_found_ids) <= 10) {
      log_message(sprintf("KO IDs not found: %s", paste(not_found_ids, collapse = ", ")), "WARN")
    } else {
      log_message(sprintf("First 10 KO IDs not found: %s...", 
                          paste(not_found_ids[1:10], collapse = ", ")), "WARN")
    }
  }
  
  # Only throw an error when all KO IDs fail
  if (success_count == 0) {
    error_msg <- paste0(
      "\n================================================================================\n",
      "FAILED TO RETRIEVE KEGG ANNOTATIONS\n",
      "================================================================================\n\n",
      "Statistics:\n",
      "  Total features queried:  ", total_features, "\n",
      "  Successful queries:      ", success_count, "\n",
      "  Not found (HTTP 404):    ", not_found_count, "\n",
      "  Failed with errors:      ", error_count, "\n\n"
    )

    if (not_found_count > 0 && error_count == 0) {
      ko_list <- if (length(not_found_ids) <= 10) {
        paste("  ", paste(not_found_ids, collapse = ", "))
      } else {
        paste("  ", paste(not_found_ids[1:10], collapse = ", "), "...")
      }

      error_msg <- paste0(error_msg,
        "Problem: All KO IDs were not found in KEGG database (HTTP 404).\n\n",
        "Possible causes:\n",
        "  1. Invalid KO IDs in your data\n",
        "  2. KEGG database structure has changed\n",
        "  3. KO IDs are outdated or deprecated\n\n",
        "KO IDs not found:\n", ko_list, "\n\n",
        "Recommendations:\n",
        "  1. Verify your KO IDs are valid (check KEGG database)\n",
        "  2. Use ko_to_kegg = FALSE for local annotation (no KEGG API needed)\n",
        "  3. Update your PICRUSt2 reference database\n\n"
      )
    } else if (error_count > 0) {
      ko_list <- if (length(error_ids) <= 10) {
        paste("  ", paste(error_ids, collapse = ", "))
      } else {
        paste("  ", paste(error_ids[1:10], collapse = ", "), "...")
      }

      error_msg <- paste0(error_msg,
        "Problem: Network or KEGG API errors.\n\n",
        "Possible causes:\n",
        "  1. No internet connection or firewall blocking KEGG\n",
        "  2. KEGG API is down or rate-limiting your requests\n",
        "  3. Network instability (especially common in China)\n\n",
        "KO IDs with errors:\n", ko_list, "\n\n",
        "Recommendations:\n",
        "  1. Check your internet connection\n",
        "  2. Try again later (KEGG API may be temporarily down)\n",
        "  3. Use VPN if in China or behind restrictive firewall\n",
        "  4. Use ko_to_kegg = FALSE for local annotation (recommended)\n\n"
      )
    } else {
      error_msg <- paste0(error_msg,
        "Problem: Unknown error occurred.\n\n",
        "Recommendations:\n",
        "  1. Check sessionInfo() for package versions\n",
        "  2. Report this issue on GitHub with reproducible example\n",
        "  3. Use ko_to_kegg = FALSE as temporary workaround\n\n"
      )
    }

    error_msg <- paste0(error_msg, "================================================================================\n")

    stop(error_msg, call. = FALSE)
  }
  
  # Check if any annotations were found
  has_annotations <- any(!is.na(filtered_df$pathway_name))
  if (!has_annotations) {
    warning("No valid KEGG annotations found for any features")
  }

  # Merge annotations back onto the full input so the return value keeps
  # every row of `daa_results_df`. Previously we returned only the
  # p_adjust < threshold subset, which silently dropped non-significant
  # features that downstream code (e.g. ggpicrust2()'s plot_result_list$
  # daa_results_df and user post-hoc analyses) expected to still be there.
  # Non-significant rows get NA annotation columns, mirroring the
  # no-significant-pathways branch above.
  new_cols <- c("pathway_name", "pathway_description", "pathway_class", "pathway_map")
  df[new_cols] <- NA_character_
  match_idx <- match(df$feature, filtered_df$feature)
  populated <- !is.na(match_idx)
  for (col in new_cols) {
    df[populated, col] <- filtered_df[[col]][match_idx[populated]]
  }

  df
}

#' Annotate pathways using reference data
#' @param data Data frame to annotate
#' @param pathway_type Type of pathway
#' @param ref_data Reference data
#' @return Annotated data frame
#' @noRd
annotate_pathways <- function(data, pathway_type, ref_data) {
  message("Starting pathway annotation...")
  
  # Check if data is from DAA results (check for feature column with any p-value column)
  is_daa_results <- "feature" %in% colnames(data) &&
    any(c("p_values", "p_adjust", "pvalue", "p.adjust") %in% colnames(data))

  # Extract features.
  # DAA results:  feature IDs live in the `feature` column (one per row).
  # File input:   the abundance data frame produced by read_abundance_file()
  #               followed by add_column(description, .after = 1) has feature
  #               IDs in column 1 (with description now in column 2 and samples
  #               in columns 3+). Either way, features are a per-row attribute.
  if (is_daa_results) {
    message("DAA results data frame is not null. Proceeding...")
    features <- data$feature
  } else {
    features <- data[[1]]
  }

  # Match features with reference data.
  if (pathway_type == "EC") {
    matches <- match(features, ref_data$id)

    unmatched <- is.na(matches)
    if (any(unmatched)) {
      features_with_prefix <- paste0("EC:", features[unmatched])
      matches_with_prefix <- match(features_with_prefix, ref_data$id)
      matches[unmatched] <- matches_with_prefix
    }
  } else if (pathway_type == "KO") {
    matches <- match(features, ref_data$id)

    unmatched <- is.na(matches)
    if (any(unmatched)) {
      features_without_prefix <- sub("^ko:", "", features[unmatched])
      matches_without_prefix <- match(features_without_prefix, ref_data$id)
      matches[unmatched] <- matches_without_prefix
    }
  } else {
    matches <- match(features, ref_data$id)
  }

  # Create description vector (one entry per feature / per row)
  descriptions <- rep(NA_character_, length(features))
  valid_matches <- !is.na(matches)
  if (any(valid_matches)) {
    descriptions[valid_matches] <- ref_data$description[matches[valid_matches]]
  }

  # Assign row-wise; works identically for DAA results and file-mode tables
  # because both have one feature per row.
  data$description <- descriptions

  message("Pathway annotation completed.")
  data
}

#' Pathway information annotation
#'
#' @description
#' This function serves two main purposes:
#' 1. Annotating pathway information from PICRUSt2 output files.
#' 2. Annotating pathway information from the output of `pathway_daa` function, and converting KO abundance to KEGG pathway abundance when `ko_to_kegg` is set to TRUE.
#'
#' **Important**: When `ko_to_kegg = TRUE`, this function automatically filters pathways by
#' `p_adjust < p_adjust_threshold`. If no pathways meet this criterion, the function returns
#' the original data with NA annotation columns and issues a detailed warning message with
#' diagnostic information and recommendations.
#'
#' @param file A character string, the path to the PICRUSt2 output file.
#' @param pathway A character string, the type of pathway to annotate. Options are "KO", "EC", or "MetaCyc".
#' @param daa_results_df A data frame, the output from `pathway_daa` function. When `ko_to_kegg = TRUE`,
#'   must contain columns: feature, p_values, p_adjust, and method.
#' @param ko_to_kegg A logical, decide if convert KO abundance to KEGG pathway abundance. Default is FALSE.
#'   Set to TRUE when using the function for the second use case. When TRUE, queries KEGG database for
#'   pathway annotations (requires internet connection) and filters for significant pathways.
#' @param organism A character string specifying the KEGG organism code (e.g., 'hsa' for human, 'eco' for E. coli).
#'   Default is NULL, which retrieves generic KO information not specific to any organism. Only used when ko_to_kegg is TRUE.
#' @param p_adjust_threshold A numeric value specifying the significance threshold for filtering
#'   pathways when `ko_to_kegg = TRUE`. Only pathways with `p_adjust < p_adjust_threshold` will be
#'   annotated via KEGG API. Default is 0.05. Ignored when `ko_to_kegg = FALSE`.
#'
#' @return A data frame with annotated pathway information.
#'
#' If using the function for the first use case (file input), the output data frame will include:
#' \itemize{
#'   \item \code{id}: The pathway ID.
#'   \item \code{description}: The description of the pathway.
#'   \item \code{sample1, sample2, ...}: Abundance values for each sample.
#' }
#'
#' If \code{ko_to_kegg} is set to TRUE, the output data frame will also include:
#' \itemize{
#'   \item \code{pathway_name}: The name of the KEGG pathway.
#'   \item \code{pathway_description}: The description of the KEGG pathway.
#'   \item \code{pathway_class}: The class of the KEGG pathway.
#'   \item \code{pathway_map}: The KEGG pathway map ID.
#' }
#'
#' **Note**: When \code{ko_to_kegg = TRUE}, only pathways with
#' \code{p_adjust < p_adjust_threshold} are processed. If no pathways meet this criterion,
#' all annotation columns will be NA, and a detailed warning message will be issued with
#' diagnostic information.
#'
#' When \code{ko_to_kegg} is TRUE, the function queries the KEGG database for pathway information.
#' By default (organism = NULL), it retrieves generic KO information that is not specific to any organism.
#' If you are interested in organism-specific pathway information, you can specify the KEGG organism code
#' using the \code{organism} parameter.
#'
#' @examples
#' \dontrun{
#' # Example 1: Annotate pathways from PICRUSt2 output file
#' pathway_annotation(file = "path/to/picrust2_output.tsv",
#'                               pathway = "KO")
#' }
#' 
#' \dontrun{
#' # Example 2: Annotate pathways from pathway_daa output
#' # Assuming you have daa_results from pathway_daa function
#' daa_results <- pathway_daa(abundance, metadata, group = "Group")
#' annotated_results <- pathway_annotation(pathway = "KO",
#'                               daa_results_df = daa_results,
#'                               ko_to_kegg = FALSE)
#' }
#' @export
pathway_annotation <- function(file = NULL,
                             pathway = NULL,
                             daa_results_df = NULL,
                             ko_to_kegg = FALSE,
                             organism = NULL,
                             p_adjust_threshold = 0.05) {
  
  # Input validation
  if (is.null(file) && is.null(daa_results_df)) {
    stop("Please input the picrust2 output or results of pathway_daa daa_results_df")
  }

  if (!is.null(daa_results_df) && nrow(daa_results_df) == 0) {
    stop("Input data frame is empty")
  }

  # ko_to_kegg is only meaningful for KO-level features: the KEGG REST API
  # we query (KEGGREST::keggGet) expects KO IDs, and the whole branch is
  # documented as "convert KO abundance to KEGG pathway abundance". If the
  # caller also passes pathway = "EC" or "MetaCyc", they have given
  # contradictory information. The previous code silently ignored `pathway`
  # and ran the KEGG path anyway, which surfaced later as either cryptic
  # HTTP 404s or mis-shaped enzyme records. Fail fast with an actionable
  # message instead.
  if (isTRUE(ko_to_kegg) && !is.null(pathway) && !identical(pathway, "KO")) {
    stop(sprintf(
      "ko_to_kegg = TRUE requires pathway = 'KO' (or NULL), got pathway = '%s'. KEGG annotation operates on KO IDs; if your features are %s IDs, set ko_to_kegg = FALSE.",
      pathway, pathway
    ), call. = FALSE)
  }
  
  # Process file input
  if (!is.null(file)) {
    # File mode uses local reference data — KEGG-specific parameters have
    # no effect here. Warn so the user doesn't silently get a different
    # annotation source than they asked for.
    if (isTRUE(ko_to_kegg)) {
      warning(
        "ko_to_kegg is ignored in file mode. ",
        "File-based annotation uses local reference data and does not query KEGG. ",
        "To get KEGG API annotations, first run pathway_daa(), then call ",
        "pathway_annotation(daa_results_df = ..., ko_to_kegg = TRUE).",
        call. = FALSE
      )
    }
    if (!is.null(organism)) {
      warning(
        "organism is ignored in file mode (only used with ko_to_kegg = TRUE on DAA results).",
        call. = FALSE
      )
    }

    abundance <- read_abundance_file(file)
    abundance <- tibble::add_column(abundance, description = NA, .after = 1)
    ref_data <- load_reference_data(pathway)
    return(annotate_pathways(abundance, pathway, ref_data))
  }
  
  # Process DAA results
  if (!is.null(daa_results_df)) {
    if (!ko_to_kegg) {
      ref_data <- load_reference_data(pathway)
      return(annotate_pathways(daa_results_df, pathway, ref_data))
    } else {
      message("KO to KEGG is set to TRUE. Proceeding with KEGG pathway annotations...")
      if (!is.null(organism)) {
        message("Using organism code: ", organism, " for species-specific pathway information.")
      } else {
        message("No organism specified. Using generic KO information across all species.")
      }
      message("Connecting to KEGG database...")
      return(process_kegg_annotations(daa_results_df, organism, p_adjust_threshold))
    }
  }
}

#' Safely Extract Elements from a List
#'
#' Safely extracts elements from a list, returning NA if the extraction fails
#'
#' @param list A list object from which to extract elements
#' @param field The name of the field to extract from the list
#' @param index The index position to extract from the field. Default is 1
#'
#' @return The extracted element if successful, NA if extraction fails
#'
#' @examples
#' # Create a sample list
#' my_list <- list(
#'   a = list(x = 1:3),
#'   b = list(y = 4:6)
#' )
#'
#' # Extract existing element
#' safe_extract(my_list, "a", 1)
#'
#' # Extract non-existing element (returns NA)
#' safe_extract(my_list, "c", 1)
#' @export
safe_extract <- function(list, field, index = 1) {
  tryCatch({
    if (is.null(list) || !field %in% names(list) || is.null(list[[field]]) || length(list[[field]]) == 0) {
      NA_character_
    } else {
      as.character(list[[field]][index])
    }
  }, error = function(e) {
    message(paste("Error in safe_extract:", e$message))
    NA_character_
  })
}
