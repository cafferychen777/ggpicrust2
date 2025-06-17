#' Read and validate input file
#' @param file_path Path to the input file
#' @return Abundance data frame
#' @noRd
read_abundance_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  file_ext <- tolower(tools::file_ext(file_path))
  valid_formats <- c("txt", "tsv", "csv")
  
  if (!file_ext %in% valid_formats) {
    stop(
      "Invalid file format. Please input file in .tsv, .txt or .csv format. ",
      "The best input file format is the output file from PICRUSt2, ",
      "specifically 'pred_metagenome_unstrat.tsv'."
    )
  }
  
  delimiter <- if (file_ext == "csv") "," else "\t"
  
  tryCatch({
    abundance <- readr::read_delim(
      file_path,
      delim = delimiter,
      show_col_types = FALSE,
      progress = FALSE
    )
    
    if (ncol(abundance) < 2) {
      stop("Input file must contain at least two columns")
    }
    
    abundance %>% tibble::add_column(description = NA, .after = 1)
  }, 
  error = function(e) {
    stop("Error reading file: ", e$message)
  })
}

#' Load reference data for pathway annotation
#' @param pathway_type One of "KO", "EC", or "MetaCyc"
#' @return Reference data list
#' @noRd
load_reference_data <- function(pathway_type) {
  if (!pathway_type %in% c("KO", "EC", "MetaCyc")) {
    stop("Invalid pathway option. Please provide one of the following options: 'KO', 'EC', 'MetaCyc'.")
  }

  ref_file <- sprintf("%s_reference.RData", pathway_type)

  # First try to load from data/ (lazy loading)
  ref_data_name <- paste0(pathway_type, "_reference")
  if (exists(ref_data_name, envir = asNamespace("ggpicrust2"))) {
    ref_data <- get(ref_data_name, envir = asNamespace("ggpicrust2"))

    # FIX: Standardize MetaCyc column names
    if (pathway_type == "MetaCyc" && all(c("X1", "X2") %in% colnames(ref_data))) {
      colnames(ref_data) <- c("id", "description")
    }

    return(ref_data)
  }

  # If not found in data/, try inst/extdata/
  ref_path <- system.file("extdata", ref_file, package = "ggpicrust2", mustWork = FALSE)

  if (file.exists(ref_path)) {
    load(ref_path)
    ref_data <- get(ref_data_name)

    # FIX: Standardize MetaCyc column names
    if (pathway_type == "MetaCyc" && all(c("X1", "X2") %in% colnames(ref_data))) {
      colnames(ref_data) <- c("id", "description")
    }

    return(ref_data)
  }

  # If still not found, try loading from the package's installed location
  ref_path <- system.file("inst/extdata", ref_file, package = "ggpicrust2", mustWork = FALSE)
  if (file.exists(ref_path)) {
    load(ref_path)
    ref_data <- get(ref_data_name)

    # FIX: Standardize MetaCyc column names
    if (pathway_type == "MetaCyc" && all(c("X1", "X2") %in% colnames(ref_data))) {
      colnames(ref_data) <- c("id", "description")
    }

    return(ref_data)
  }

  # If we reach here, the file was not found in any location
  stop(sprintf("Reference data file '%s' not found in any standard location.\nPlease ensure the package was installed correctly.", ref_file))

  # This code is unreachable, just for reference
  ref_data
}

#' Cache manager for KEGG annotations
#' @noRd
kegg_cache <- new.env(parent = emptyenv())

#' Get KEGG annotation with caching
#' @param ko_id KO identifier
#' @param organism Organism code (e.g., 'hsa' for human, 'eco' for E. coli). If NULL, the generic KO entry is retrieved.
#' @return KEGG annotation
#' @noRd
get_kegg_with_cache <- function(ko_id, organism = NULL) {
  # Create a cache key that includes the organism if specified
  cache_key <- if (!is.null(organism)) {
    paste0(organism, ":", ko_id)
  } else {
    ko_id
  }
  
  if (!exists(cache_key, envir = kegg_cache)) {
    # Add request interval
    Sys.sleep(0.1)  # 100ms delay to avoid too frequent requests
    
    # Different query strategy based on whether organism is specified
    if (!is.null(organism)) {
      # For organism-specific queries, we need a two-step process:
      # 1. Find organism-specific genes linked to this KO
      # 2. Get pathway information for those genes
      
      # Step 1: Find organism-specific genes
      organism_genes <- tryCatch({
        # Use keggLink to find genes in the specified organism that are linked to this KO
        links <- KEGGREST::keggLink(organism, ko_id)
        links
      }, error = function(e) {
        log_message(sprintf("Error finding %s genes for %s: %s", organism, ko_id, e$message), "WARN")
        NULL
      })
      
      # If we found organism-specific genes
      if (!is.null(organism_genes) && length(organism_genes) > 0) {
        # Step 2: Get pathway information for the first gene (as a representative)
        # We use the first gene as most genes linked to the same KO participate in the same pathways
        gene_id <- names(organism_genes)[1]
        
        result <- with_retry({
          KEGGREST::keggGet(gene_id)
        })
        
        # If successful, store in cache
        if (!is.null(result) && !inherits(result, "kegg_error")) {
          assign(cache_key, result, envir = kegg_cache)
        } else {
          # If we couldn't get gene info, fall back to generic KO query
          log_message(sprintf("Could not get info for %s gene %s, falling back to generic KO query", 
                             organism, gene_id), "INFO")
          
          result <- with_retry({
            KEGGREST::keggGet(ko_id)
          })
          
          if (!is.null(result) && !inherits(result, "kegg_error")) {
            assign(cache_key, result, envir = kegg_cache)
          }
        }
      } else {
        # If no organism-specific genes found, fall back to generic KO query
        log_message(sprintf("No %s genes found for %s, falling back to generic KO query", 
                           organism, ko_id), "INFO")
        
        result <- with_retry({
          KEGGREST::keggGet(ko_id)
        })
        
        if (!is.null(result) && !inherits(result, "kegg_error")) {
          assign(cache_key, result, envir = kegg_cache)
        }
      }
    } else {
      # For generic KO queries, directly query KEGG
      result <- with_retry({
        KEGGREST::keggGet(ko_id)
      })
      
      # Handle error objects
      if (inherits(result, "kegg_error")) {
        result$ko_id <- ko_id
        return(result)
      }
      
      if (!is.null(result)) {
        assign(cache_key, result, envir = kegg_cache)
      }
    }
  }
  
  # If exists in cache, retrieve it
  if (exists(cache_key, envir = kegg_cache)) {
    return(get(cache_key, envir = kegg_cache, inherits = FALSE))
  }
  
  # If not in cache, return NULL
  NULL
}

#' Retry mechanism for KEGG queries
#' @param expr Expression or function to retry
#' @param max_attempts Maximum number of retry attempts
#' @return Result or error
#' @noRd
with_retry <- function(expr, max_attempts = getOption("ggpicrust2.max_retries", 3)) {
  attempt <- 1
  
  while (attempt <= max_attempts) {
    result <- tryCatch({
      if (is.function(expr)) expr() else if (is.expression(expr)) eval(expr) else expr
    }, error = function(e) {
      # Special handling for HTTP 404 errors
      if (grepl("HTTP 404", e$message)) {
        # Return a special error object instead of NULL
        # This allows the caller to distinguish between "resource not found" and "query failed"
        return(structure(list(
          status = "not_found",
          message = e$message,
          ko_id = NULL  # Caller can populate this
        ), class = "kegg_error"))
      }
      
      # Handle other errors
      if (attempt == max_attempts) {
        # Return a special error object
        return(structure(list(
          status = "failed",
          message = e$message,
          ko_id = NULL
        ), class = "kegg_error"))
      }
      
      log_message(sprintf("Attempt %d failed: %s, retrying...", attempt, e$message), "WARN")
      Sys.sleep(min(2^attempt, 30))  # Limit maximum wait time to 30 seconds
      NULL
    })
    
    # If result is not NULL or is a special error object, return it
    if (!is.null(result) || inherits(result, "kegg_error")) {
      return(result)
    }
    
    attempt <- attempt + 1
  }
  
  # If all attempts fail, return NULL
  NULL
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

#' Validate input parameters
#' @noRd
validate_inputs <- function(file, pathway, daa_results_df, ko_to_kegg) {
  # File validation
  if (!is.null(file)) {
    if (!is.character(file) || length(file) != 1) {
      stop("'file' must be a single character string")
    }
  }
  
  # Pathway validation
  if (!is.null(pathway)) {
    valid_pathways <- c("KO", "EC", "MetaCyc")
    if (!pathway %in% valid_pathways) {
      stop(sprintf("'pathway' must be one of: %s", 
                  paste(valid_pathways, collapse = ", ")))
    }
  }
  
  # Data frame validation
  if (!is.null(daa_results_df)) {
    if (!is.data.frame(daa_results_df)) {
      stop("'daa_results_df' must be a data frame")
    }
    if (ko_to_kegg) {
      required_cols <- c("feature", "p_adjust")
      missing_cols <- setdiff(required_cols, colnames(daa_results_df))
      if (length(missing_cols) > 0) {
        stop(sprintf("Missing required columns: %s", 
                    paste(missing_cols, collapse = ", ")))
      }
    }
  }
}

#' Process KEGG annotations with improved error handling and caching
#' @param df Data frame with features to annotate
#' @param organism KEGG organism code (e.g., 'hsa' for human)
#' @return Annotated data frame
#' @noRd
process_kegg_annotations <- function(df, organism = NULL) {
  if (nrow(df) == 0) {
    stop("Empty data frame provided for KEGG annotation")
  }
  
  filtered_df <- df[df$p_adjust < 0.05, ]
  if (nrow(filtered_df) == 0) {
    stop(
      "No statistically significant biomarkers found (p_adjust < 0.05).\n",
      "Consider using a less stringent threshold or reviewing your data."
    )
  }
  
  # Initialize new columns
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
      filtered_df$pathway_description[i] <- if("DESCRIPTION" %in% names(entry[[1]])) safe_extract(entry[[1]], "DESCRIPTION", 1) else NA_character_
      filtered_df$pathway_class[i] <- if("CLASS" %in% names(entry[[1]])) safe_extract(entry[[1]], "CLASS", 1) else NA_character_
      filtered_df$pathway_map[i] <- if("PATHWAY_MAP" %in% names(entry[[1]])) safe_extract(entry[[1]], "PATHWAY_MAP", 1) else NA_character_
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
    if (not_found_count > 0 && error_count == 0) {
      stop("All KO IDs were not found in the KEGG database (HTTP 404). ",
           "This could be due to invalid KO IDs or KEGG database changes.")
    } else if (error_count > 0) {
      stop("Failed to retrieve any KEGG annotations due to errors. ",
           "This could be due to network issues or KEGG API changes.")
    } else {
      stop("Failed to retrieve any KEGG annotations for unknown reasons.")
    }
  }
  
  # Check if any annotations were found
  has_annotations <- any(!is.na(filtered_df$pathway_name))
  if (!has_annotations) {
    warning("No valid KEGG annotations found for any features")
  }
  
  filtered_df
}

#' Annotate pathways using reference data
#' @param data Data frame to annotate
#' @param pathway_type Type of pathway
#' @param ref_data Reference data
#' @return Annotated data frame
#' @noRd
annotate_pathways <- function(data, pathway_type, ref_data) {
  message("Starting pathway annotation...")
  
  # Check if data is from DAA results
  is_daa_results <- all(c("feature", "p_values") %in% colnames(data))
  
  # Extract features
  if (is_daa_results) {
    message("DAA results data frame is not null. Proceeding...")
    features <- data$feature
  } else {
    features <- colnames(data)[-c(1, 2)]  # Skip first two columns (typically ID and description)
  }
  
  # Match features with reference data
  matches <- match(features, ref_data$id)
  
  # Create description column
  descriptions <- rep(NA_character_, length(features))
  valid_matches <- !is.na(matches)
  if (any(valid_matches)) {
    descriptions[valid_matches] <- ref_data$description[matches[valid_matches]]
  }
  
  # Update data frame
  if (is_daa_results) {
    # For DAA results, add description column
    data$description <- descriptions
  } else {
    # For abundance data, update description column
    for (i in seq_along(features)) {
      feature_idx <- which(colnames(data) == features[i])
      if (length(feature_idx) > 0) {
        data$description[feature_idx - 2] <- descriptions[i]  # Adjust for offset
      }
    }
  }
  
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
#' @param file A character string, the path to the PICRUSt2 output file.
#' @param pathway A character string, the type of pathway to annotate. Options are "KO", "EC", or "MetaCyc".
#' @param daa_results_df A data frame, the output from `pathway_daa` function.
#' @param ko_to_kegg A logical, decide if convert KO abundance to KEGG pathway abundance. Default is FALSE. Set to TRUE when using the function for the second use case.
#' @param organism A character string specifying the KEGG organism code (e.g., 'hsa' for human, 'eco' for E. coli). Default is NULL, which retrieves generic KO information not specific to any organism. Only used when ko_to_kegg is TRUE.
#'
#' @return A data frame with annotated pathway information. 
#' If using the function for the first use case, the output data frame will include the following columns:
#' \itemize{
#'   \item \code{id}: The pathway ID.
#'   \item \code{description}: The description of the pathway.
#'   \item \code{sample1, sample2, ...}: Abundance values for each sample.
#' }
#' 
#' If \code{ko_to_kegg} is set to TRUE, the output data frame will also include the following columns:
#' \itemize{
#'   \item \code{pathway_name}: The name of the KEGG pathway.
#'   \item \code{pathway_description}: The description of the KEGG pathway.
#'   \item \code{pathway_class}: The class of the KEGG pathway.
#'   \item \code{pathway_map}: The KEGG pathway map ID.
#' }
#' 
#' When \code{ko_to_kegg} is TRUE, the function queries the KEGG database for pathway information. By default (organism = NULL), it retrieves generic KO information that is not specific to any organism. If you are interested in organism-specific pathway information, you can specify the KEGG organism code using the \code{organism} parameter.
#'
#' @examples
#' \donttest{
#' # Example 1: Annotate pathways from PICRUSt2 output file
#' pathway_annotation(file = "path/to/picrust2_output.tsv",
#'                               pathway = "KO")
#'
#' # Example 2: Annotate pathways from pathway_daa output
#' # and converting KO abundance to KEGG pathway abundance
#' daa_results <- pathway_daa(abundance, metadata, group = "Group")
#' annotated_results <- pathway_annotation(pathway = "KO",
#'                               daa_results_df = daa_results,
#'                               ko_to_kegg = TRUE)
#'
#' # Example 3: Annotate EC pathways
#' ec_results <- pathway_daa(abundance, metadata, group = "Group")
#' annotated_ec <- pathway_annotation(pathway = "EC",
#'                               daa_results_df = ec_results,
#'                               ko_to_kegg = TRUE)
#'
#' # Example 4: Annotate KO pathways with human-specific information
#' daa_results <- pathway_daa(abundance, metadata, group = "Group")
#' human_results <- pathway_annotation(pathway = "KO",
#'                               daa_results_df = daa_results,
#'                               ko_to_kegg = TRUE,
#'                               organism = "hsa")
#' }
#' @export
pathway_annotation <- function(file = NULL,
                             pathway = NULL,
                             daa_results_df = NULL,
                             ko_to_kegg = FALSE,
                             organism = NULL) {
  
  # Input validation
  if (is.null(file) && is.null(daa_results_df)) {
    stop("Please input the picrust2 output or results of pathway_daa daa_results_df")
  }
  
  if (!is.null(daa_results_df) && nrow(daa_results_df) == 0) {
    stop("Input data frame is empty")
  }
  
  # Process file input
  if (!is.null(file)) {
    abundance <- read_abundance_file(file)
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
      message("We are connecting to the KEGG database to get the latest results, please wait patiently.")
      message("Processing pathways in chunks...")
      return(process_kegg_annotations(daa_results_df, organism))
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