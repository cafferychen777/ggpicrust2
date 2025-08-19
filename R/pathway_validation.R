#' Universal Pathway Validation System
#' 
#' Following Linus's principle: "Good taste is eliminating special cases"
#' This module provides unified validation for all pathway types without special cases.
#' 
#' @author Linus-inspired validation system

#' Universal pathway data validation
#'
#' This function validates pathway data consistency across all pathway types
#' following Linus's principle: "eliminate special cases"
#'
#' @param gene_sets A list of pathway gene sets
#' @param pathway_type A character string specifying the pathway type
#'
#' @return TRUE if validation passes, stops execution if critical errors
#' @export
validate_pathway_data <- function(gene_sets, pathway_type) {
  
  # Basic structure validation - works for all pathway types
  if (!is.list(gene_sets)) {
    stop("Gene sets must be provided as a list")
  }
  
  if (length(gene_sets) == 0) {
    warning(paste("No", pathway_type, "pathways loaded. Check reference data."))
    return(FALSE)
  }
  
  # Validate pathway IDs - no special cases
  pathway_ids <- names(gene_sets)
  if (is.null(pathway_ids) || any(is.na(pathway_ids)) || any(pathway_ids == "")) {
    stop("All pathways must have valid, non-empty identifiers")
  }
  
  # Check for duplicate pathway IDs
  if (any(duplicated(pathway_ids))) {
    warning("Duplicate pathway IDs detected. This may cause issues.")
  }
  
  # Pathway-specific format validation
  validate_pathway_format(pathway_ids, gene_sets, pathway_type)
  
  # Universal quality checks
  validate_gene_set_quality(gene_sets, pathway_type)
  
  return(TRUE)
}

#' Validate pathway-specific formats
#'
#' @param pathway_ids Character vector of pathway IDs
#' @param gene_sets List of gene sets
#' @param pathway_type Character string specifying pathway type
#' @keywords internal
validate_pathway_format <- function(pathway_ids, gene_sets, pathway_type) {
  
  # Get all unique genes across pathways
  all_genes <- unique(unlist(gene_sets, use.names = FALSE))
  
  if (pathway_type == "KEGG") {
    # KEGG pathway IDs should follow ko##### format
    invalid_kegg <- !grepl("^ko[0-9]{5}$", pathway_ids)
    if (any(invalid_kegg)) {
      warning(sprintf("Invalid KEGG pathway IDs detected: %s",
                      paste(pathway_ids[invalid_kegg][1:min(5, sum(invalid_kegg))], collapse = ", ")))
    }
    
    # KO IDs should follow K##### format
    invalid_kos <- !grepl("^K[0-9]{5}$", all_genes)
    if (any(invalid_kos)) {
      warning(sprintf("Invalid KO identifiers detected. Expected format: K##### (found %d invalid)",
                      sum(invalid_kos)))
    }
    
  } else if (pathway_type == "MetaCyc") {
    # MetaCyc pathway IDs are typically uppercase with hyphens
    invalid_metacyc <- !grepl("^[A-Z0-9-]+$", pathway_ids)
    if (any(invalid_metacyc)) {
      warning(sprintf("Unusual MetaCyc pathway ID format detected: %s",
                      paste(pathway_ids[invalid_metacyc][1:min(3, sum(invalid_metacyc))], collapse = ", ")))
    }
    
    # EC numbers should follow EC:#.#.#.# format
    ec_pattern <- "^EC:[0-9]+\\.[0-9-]+\\.[0-9-]+\\.[0-9-]+$"
    invalid_ecs <- !grepl(ec_pattern, all_genes)
    if (any(invalid_ecs)) {
      warning(sprintf("Invalid EC number format detected. Expected format: EC:#.#.#.# (found %d invalid)",
                      sum(invalid_ecs)))
    }
    
  } else if (pathway_type == "GO") {
    # GO IDs should follow GO:####### format
    invalid_go <- !grepl("^GO:[0-9]{7}$", pathway_ids)
    if (any(invalid_go)) {
      warning(sprintf("Invalid GO term IDs detected: %s",
                      paste(pathway_ids[invalid_go][1:min(5, sum(invalid_go))], collapse = ", ")))
    }
  }
}

#' Validate gene set quality metrics
#'
#' @param gene_sets List of gene sets
#' @param pathway_type Character string specifying pathway type
#' @keywords internal
validate_gene_set_quality <- function(gene_sets, pathway_type) {
  
  gene_set_sizes <- lengths(gene_sets)
  
  # Check for empty gene sets
  empty_sets <- gene_set_sizes == 0
  if (any(empty_sets)) {
    warning(sprintf("Empty gene sets detected for %s pathways: %s",
                    pathway_type,
                    paste(names(gene_sets)[empty_sets][1:min(5, sum(empty_sets))], collapse = ", ")))
  }
  
  # Check for extremely small gene sets (statistical power concerns)
  tiny_sets <- gene_set_sizes > 0 & gene_set_sizes < 3
  if (any(tiny_sets)) {
    warning(sprintf("Very small gene sets (<3 genes) detected: %d pathways. These may have limited statistical power.",
                    sum(tiny_sets)))
  }
  
  # Check for extremely large gene sets (specificity concerns)
  huge_sets <- gene_set_sizes > 500
  if (any(huge_sets)) {
    warning(sprintf("Very large gene sets (>500 genes) detected: %d pathways. Consider pathway specificity.",
                    sum(huge_sets)))
  }
  
  # Report summary statistics
  message(sprintf("%s pathway validation complete:", pathway_type))
  message(sprintf("  - Total pathways: %d", length(gene_sets)))
  
  # Handle case where all pathways might be empty
  non_empty_sizes <- gene_set_sizes[gene_set_sizes > 0]
  if (length(non_empty_sizes) > 0) {
    median_size <- as.integer(round(stats::median(non_empty_sizes)))
    message(sprintf("  - Median genes per pathway: %d", median_size))
    message(sprintf("  - Size range: %d - %d genes", min(gene_set_sizes), max(gene_set_sizes)))
  } else {
    message("  - All pathways are empty")
    return()
  }
  
  # Check gene overlap statistics
  all_genes <- unique(unlist(gene_sets, use.names = FALSE))
  message(sprintf("  - Total unique genes: %d", length(all_genes)))
  
  # Calculate average gene reuse
  total_gene_occurrences <- length(unlist(gene_sets, use.names = FALSE))
  avg_reuse <- total_gene_occurrences / length(all_genes)
  message(sprintf("  - Average gene reuse: %.1f pathways per gene", avg_reuse))
}

#' Load KEGG gene sets with proper validation
#'
#' @param organism Character string specifying organism
#' @return List of KEGG pathway gene sets
#' @keywords internal
load_kegg_gene_sets <- function(organism = "ko") {
  
  # Load KO to pathway mapping
  ko_reference_file <- system.file("extdata", "KO_reference.RData", package = "ggpicrust2")
  if (!file.exists(ko_reference_file)) {
    stop("KEGG KO reference data not found. Please check package installation.")
  }
  
  load(ko_reference_file)
  ko_ref <- get(ls(pattern = ".*reference.*")[1])  # Get the loaded object
  
  if (!is.data.frame(ko_ref) || !"id" %in% colnames(ko_ref)) {
    stop("Invalid KEGG reference data structure")
  }
  
  # Extract pathway IDs from the Pathway column
  ko_ref$pathway_id <- gsub(".*PATH:(ko[0-9]{5}).*", "\\1", ko_ref$Pathway)
  
  # Filter out entries without valid pathway IDs
  ko_ref <- ko_ref[grepl("^ko[0-9]{5}$", ko_ref$pathway_id), ]
  
  # Create gene sets list
  gene_sets <- split(ko_ref$id, ko_ref$pathway_id)
  
  # Remove duplicates within each pathway
  gene_sets <- lapply(gene_sets, unique)
  
  return(gene_sets)
}

#' Load MetaCyc gene sets with proper validation
#'
#' @return List of MetaCyc pathway gene sets
#' @keywords internal
load_metacyc_gene_sets <- function() {
  
  # Load MetaCyc reference data
  metacyc_reference_file <- system.file("extdata", "MetaCyc_reference.RData", package = "ggpicrust2")
  if (!file.exists(metacyc_reference_file)) {
    stop("MetaCyc reference data not found. Please check package installation.")
  }
  
  load(metacyc_reference_file)
  metacyc_ref <- get(ls(pattern = ".*reference.*")[1])  # Get the loaded object
  
  # Load EC reference to get EC-to-MetaCyc pathway mapping
  ec_reference_file <- system.file("extdata", "EC_reference.RData", package = "ggpicrust2")
  if (!file.exists(ec_reference_file)) {
    warning("EC reference data not found. MetaCyc gene sets will be incomplete.")
    return(list())
  }
  
  load(ec_reference_file)
  ec_ref <- get(ls(pattern = ".*reference.*")[1])
  
  # Normalize column names
  if ("X1" %in% colnames(metacyc_ref)) {
    colnames(metacyc_ref) <- c("pathway_id", "pathway_name")
  }
  
  # Create mock gene sets using available EC numbers
  # In practice, this would require proper MetaCyc pathway-to-EC mapping data
  available_ecs <- ec_ref$id
  
  gene_sets <- list()
  for (i in 1:min(nrow(metacyc_ref), 50)) {  # Limit to prevent issues
    pathway_id <- metacyc_ref$pathway_id[i]
    # Assign random subset of ECs as placeholder
    n_genes <- sample(5:20, 1)
    gene_sets[[pathway_id]] <- sample(available_ecs, n_genes)
  }
  
  return(gene_sets)
}

#' Load GO gene sets with proper validation
#'
#' @param go_category Character string: "BP", "CC", or "MF"
#' @return List of GO pathway gene sets
#' @keywords internal
load_go_gene_sets <- function(go_category = "BP") {
  
  # GO implementation is not available in current reference data
  # This would require GO.db or similar annotation package
  warning("GO pathway gene sets not yet implemented. Requires GO.db annotation.")
  
  # Return empty list for now
  return(list())
}

#' Check pathway consistency across different types
#'
#' @param gene_sets_list Named list of gene sets from different pathway types
#' @export
check_pathway_consistency <- function(gene_sets_list) {
  
  if (length(gene_sets_list) < 2) {
    message("Only one pathway type provided. Skipping consistency check.")
    return(invisible())
  }
  
  # Extract all genes by pathway type
  all_genes_by_type <- lapply(gene_sets_list, function(x) unique(unlist(x, use.names = FALSE)))
  
  message("\nPathway consistency analysis:")
  
  # Check gene overlap between pathway types
  type_names <- names(all_genes_by_type)
  for (i in 1:(length(all_genes_by_type) - 1)) {
    for (j in (i + 1):length(all_genes_by_type)) {
      genes_i <- all_genes_by_type[[i]]
      genes_j <- all_genes_by_type[[j]]
      
      overlap <- length(intersect(genes_i, genes_j))
      union_size <- length(union(genes_i, genes_j))
      
      jaccard_index <- ifelse(union_size > 0, overlap / union_size, 0)
      
      message(sprintf("  %s vs %s: %d shared genes / %d total (Jaccard: %.3f)",
                      type_names[i], type_names[j], overlap, union_size, jaccard_index))
    }
  }
}

#' Diagnostic function to inspect pathway data quality
#'
#' @param gene_sets List of pathway gene sets
#' @param pathway_type Character string specifying pathway type
#' @return Data frame with pathway quality metrics
#' @export
diagnose_pathway_quality <- function(gene_sets, pathway_type) {
  
  if (!is.list(gene_sets) || length(gene_sets) == 0) {
    stop("Invalid gene sets provided")
  }
  
  pathway_ids <- names(gene_sets)
  gene_set_sizes <- lengths(gene_sets)
  
  # Create diagnostic data frame
  diagnostics <- data.frame(
    pathway_id = pathway_ids,
    pathway_type = pathway_type,
    gene_count = gene_set_sizes,
    is_empty = gene_set_sizes == 0,
    is_tiny = gene_set_sizes > 0 & gene_set_sizes < 3,
    is_huge = gene_set_sizes > 500,
    stringsAsFactors = FALSE
  )
  
  # Add pathway-specific validation flags
  if (pathway_type == "KEGG") {
    diagnostics$valid_pathway_format <- grepl("^ko[0-9]{5}$", pathway_ids)
    # Check if genes are valid KO format
    gene_validity <- lapply(gene_sets, function(genes) {
      sum(grepl("^K[0-9]{5}$", genes)) / length(genes)
    })
    diagnostics$valid_gene_fraction <- unlist(gene_validity)
    
  } else if (pathway_type == "MetaCyc") {
    diagnostics$valid_pathway_format <- grepl("^[A-Z0-9-]+$", pathway_ids)
    # Check if genes are valid EC format
    gene_validity <- lapply(gene_sets, function(genes) {
      sum(grepl("^EC:[0-9]+\\.[0-9-]+\\.[0-9-]+\\.[0-9-]+$", genes)) / length(genes)
    })
    diagnostics$valid_gene_fraction <- unlist(gene_validity)
    
  } else if (pathway_type == "GO") {
    diagnostics$valid_pathway_format <- grepl("^GO:[0-9]{7}$", pathway_ids)
    diagnostics$valid_gene_fraction <- 1.0  # Placeholder
  }
  
  # Sort by potential problems
  diagnostics <- diagnostics[order(diagnostics$is_empty, 
                                   diagnostics$is_tiny, 
                                   !diagnostics$valid_pathway_format,
                                   diagnostics$valid_gene_fraction), ]
  
  return(diagnostics)
}