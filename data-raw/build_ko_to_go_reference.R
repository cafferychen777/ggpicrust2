# Build KO-to-GO reference data for ggpicrust2
#
# Data provenance (first principles):
#   1. KEGG REST API — each KO entry's DBLINKS section lists curated GO IDs
#   2. EBI QuickGO API — authoritative GO term names and categories (BP/MF/CC)
#
# No hardcoded or manually curated KO-GO mappings are used.
# All data comes directly from authoritative biological databases.
#
# The KEGG phase is slow (~35 min for ~5000 KOs at ~4s per batch).
# Intermediate results are cached to data-raw/.cache/ so the build
# can be interrupted and resumed without losing progress.
#
# Usage:
#   source("data-raw/build_ko_to_go_reference.R")
#   build_ko_to_go_reference()               # full build (resumable)
#   build_ko_to_go_reference(max_kos = 200)  # limited build for testing

source("data-raw/build_utils.R")

# =============================================================================
# Configuration
# =============================================================================

KEGG_BATCH_SIZE    <- 10    # KEGG /get accepts up to 10 entries per request
KEGG_DELAY         <- 0.4   # seconds between KEGG requests (rate limit)
QUICKGO_BATCH_SIZE <- 25    # QuickGO comma-separated limit
QUICKGO_DELAY      <- 0.15  # seconds between QuickGO requests
MIN_KO_PER_TERM    <- 3     # minimum KO members to include a GO term

# =============================================================================
# Phase 1: KEGG — Extract GO IDs from KO entry DBLINKS
# =============================================================================

#' Parse GO IDs from a batch KEGG /get response.
#' @return Named list: ko_id -> character vector of GO IDs.
parse_kegg_batch <- function(text) {
  entries <- strsplit(text, "///")[[1]]
  results <- list()

  for (entry in entries) {
    entry <- trimws(entry)
    if (nchar(entry) == 0) next

    ko_match <- regmatches(entry, regexpr("ENTRY\\s+(K\\d{5})", entry, perl = TRUE))
    if (length(ko_match) == 0) next
    ko_id <- sub("ENTRY\\s+", "", ko_match)

    go_line <- regmatches(entry, regexpr("GO:\\s+[0-9 ]+", entry, perl = TRUE))
    if (length(go_line) == 0) next

    go_nums <- regmatches(go_line, gregexpr("[0-9]{7}", go_line))[[1]]
    if (length(go_nums) > 0) {
      results[[ko_id]] <- paste0("GO:", go_nums)
    }
  }
  results
}

#' Batch-query KEGG for KO entries, extract GO cross-references.
#' Uses caching so partial runs can be resumed.
#' @return Named list: ko_id -> character vector of GO IDs.
collect_kegg_go_mappings <- function(ko_list) {
  cache_name <- sprintf("kegg_go_mappings_%d", length(ko_list))
  cached <- cache_load(cache_name, max_age_hours = 48)
  if (!is.null(cached)) return(cached)

  n_batches <- ceiling(length(ko_list) / KEGG_BATCH_SIZE)
  message(sprintf("Querying KEGG for %d KOs in %d batches (~%d min)...",
                  length(ko_list), n_batches,
                  ceiling(n_batches * (4 + KEGG_DELAY) / 60)))

  all_mappings <- list()

  for (i in seq_len(n_batches)) {
    start <- (i - 1) * KEGG_BATCH_SIZE + 1
    end   <- min(i * KEGG_BATCH_SIZE, length(ko_list))
    batch <- ko_list[start:end]

    if (i %% 25 == 0 || i == n_batches) {
      show_progress(i, n_batches,
                    sprintf("%d KOs with GO", length(all_mappings)))
    }

    url <- sprintf("%s/get/%s", BUILD_CONFIG$kegg_base,
                   paste(paste0("ko:", batch), collapse = "+"))
    text <- safe_get(url, as = "text")

    if (!is.null(text)) {
      batch_results <- parse_kegg_batch(text)
      all_mappings <- c(all_mappings, batch_results)
    }

    # Save progress every 100 batches
    if (i %% 100 == 0) {
      cache_save(all_mappings, cache_name)
    }

    Sys.sleep(KEGG_DELAY)
  }

  message(sprintf("KEGG: %d/%d KOs have GO cross-references.",
                  length(all_mappings), length(ko_list)))

  cache_save(all_mappings, cache_name)
  all_mappings
}

# =============================================================================
# Phase 2: QuickGO — Retrieve GO term names and categories
# =============================================================================

#' Map QuickGO aspect string to GO category code.
aspect_to_category <- function(aspect) {
  switch(aspect,
         "biological_process"  = "BP",
         "molecular_function"  = "MF",
         "cellular_component"  = "CC",
         NA_character_)
}

#' Batch-query QuickGO for GO term metadata.
#' @return data.frame: go_id, go_name, category.
fetch_go_metadata <- function(go_ids) {
  cache_name <- sprintf("quickgo_metadata_%d", length(go_ids))
  cached <- cache_load(cache_name, max_age_hours = 48)
  if (!is.null(cached)) return(cached)

  n_batches <- ceiling(length(go_ids) / QUICKGO_BATCH_SIZE)
  message(sprintf("Querying QuickGO for %d GO terms in %d batches...",
                  length(go_ids), n_batches))

  rows <- vector("list", length(go_ids))
  idx <- 0

  for (i in seq_len(n_batches)) {
    start <- (i - 1) * QUICKGO_BATCH_SIZE + 1
    end   <- min(i * QUICKGO_BATCH_SIZE, length(go_ids))
    batch <- go_ids[start:end]

    url <- sprintf("%s/%s", BUILD_CONFIG$quickgo_base,
                   paste(batch, collapse = ","))
    data <- safe_get(url, accept = "application/json", as = "parsed")

    if (!is.null(data) && !is.null(data$results)) {
      for (j in seq_len(nrow(data$results))) {
        if (isTRUE(data$results$isObsolete[j])) next
        category <- aspect_to_category(data$results$aspect[j])
        if (is.na(category)) next

        idx <- idx + 1
        rows[[idx]] <- data.frame(
          go_id    = data$results$id[j],
          go_name  = data$results$name[j],
          category = category,
          stringsAsFactors = FALSE
        )
      }
    }

    Sys.sleep(QUICKGO_DELAY)
  }

  go_meta <- do.call(rbind, rows[seq_len(idx)])
  message(sprintf("QuickGO: metadata for %d/%d terms (excluding obsolete).",
                  nrow(go_meta), length(go_ids)))

  cache_save(go_meta, cache_name)
  go_meta
}

# =============================================================================
# Phase 3: Aggregate, validate, save
# =============================================================================

#' Aggregate KO-GO mappings and GO metadata into reference format.
aggregate_reference <- function(ko_go_map, go_meta) {
  long_df <- do.call(rbind, lapply(names(ko_go_map), function(ko) {
    data.frame(ko_id = ko, go_id = ko_go_map[[ko]], stringsAsFactors = FALSE)
  }))

  message(sprintf("Total KO-GO pairs: %d", nrow(long_df)))
  long_df <- long_df[long_df$go_id %in% go_meta$go_id, ]

  go_groups <- split(long_df$ko_id, long_df$go_id)

  rows <- lapply(names(go_groups), function(go_id) {
    kos <- sort(unique(go_groups[[go_id]]))
    if (length(kos) < MIN_KO_PER_TERM) return(NULL)

    meta_row <- go_meta[go_meta$go_id == go_id, ][1, ]
    data.frame(
      go_id      = go_id,
      go_name    = meta_row$go_name,
      category   = meta_row$category,
      ko_members = paste(kos, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })

  valid_rows <- Filter(Negate(is.null), rows)
  if (length(valid_rows) == 0) {
    stop(sprintf("No GO terms have >= %d KO members. Try querying more KOs.",
                 MIN_KO_PER_TERM), call. = FALSE)
  }

  ref <- do.call(rbind, valid_rows)
  rownames(ref) <- NULL
  ref <- ref[order(as.character(ref$category), as.character(ref$go_id)), ]

  stopifnot(all(grepl("^GO:\\d{7}$", ref$go_id)))
  stopifnot(all(ref$category %in% c("BP", "MF", "CC")))

  message(sprintf("Kept %d GO terms (>= %d KO members each).", nrow(ref), MIN_KO_PER_TERM))
  ref
}

# =============================================================================
# Main entry point
# =============================================================================

#' Build ko_to_go_reference from authoritative sources.
#'
#' @param max_kos Maximum KOs to query (NULL = all). Use 200 for quick testing.
#' @param use_cache Whether to use cached intermediate results (default TRUE).
#' @return data.frame: go_id, go_name, category, ko_members (invisible).
build_ko_to_go_reference <- function(max_kos = NULL, use_cache = TRUE) {
  message("=== Building ko_to_go_reference ===")
  message("Sources: KEGG REST API (DBLINKS) + EBI QuickGO\n")

  if (!use_cache) cache_clear()

  # Get KO list from package data
  data("ko_abundance", package = "ggpicrust2")
  ko_list <- ko_abundance[["#NAME"]]
  ko_list <- ko_list[grepl("^K[0-9]{5}$", ko_list)]
  message(sprintf("Found %d valid KO IDs in ko_abundance.", length(ko_list)))

  if (!is.null(max_kos) && length(ko_list) > max_kos) {
    message(sprintf("Limiting to first %d KOs.", max_kos))
    ko_list <- ko_list[seq_len(max_kos)]
  }

  # Phase 1: KEGG
  message("\n--- Phase 1: KEGG DBLINKS ---")
  ko_go_map <- collect_kegg_go_mappings(ko_list)

  if (length(ko_go_map) == 0) {
    stop("No KO-GO mappings found. Check network connectivity.", call. = FALSE)
  }

  unique_go_ids <- unique(unlist(ko_go_map, use.names = FALSE))
  message(sprintf("Unique GO terms to annotate: %d\n", length(unique_go_ids)))

  # Phase 2: QuickGO
  message("--- Phase 2: QuickGO metadata ---")
  go_meta <- fetch_go_metadata(unique_go_ids)

  if (is.null(go_meta) || nrow(go_meta) == 0) {
    stop("No GO metadata retrieved. Check network connectivity.", call. = FALSE)
  }

  # Phase 3: Aggregate
  message("\n--- Phase 3: Aggregate and validate ---")
  ko_to_go_reference <- aggregate_reference(ko_go_map, go_meta)

  # Metadata
  attr(ko_to_go_reference, "creation_date") <- Sys.Date()
  attr(ko_to_go_reference, "data_sources") <- "KEGG REST API (DBLINKS) + EBI QuickGO"
  attr(ko_to_go_reference, "version") <- "3.0"
  attr(ko_to_go_reference, "total_kos_queried") <- length(ko_list)
  attr(ko_to_go_reference, "min_ko_per_term") <- MIN_KO_PER_TERM

  validate_reference(ko_to_go_reference,
                     c("go_id", "go_name", "category", "ko_members"),
                     "ko_to_go_reference")
  save_reference(ko_to_go_reference, "ko_to_go_reference")

  # Report
  cat_counts <- as.integer(table(factor(ko_to_go_reference$category,
                                        levels = c("BP", "MF", "CC"))))
  ko_counts <- sapply(strsplit(ko_to_go_reference$ko_members, ";"), length)
  message(sprintf("\n=== Results ==="))
  message(sprintf("GO terms: %d (BP: %d, MF: %d, CC: %d)",
                  nrow(ko_to_go_reference),
                  cat_counts[1], cat_counts[2], cat_counts[3]))
  message(sprintf("KO members per term: min=%d, median=%.0f, max=%d",
                  min(ko_counts), median(ko_counts), max(ko_counts)))

  invisible(ko_to_go_reference)
}
