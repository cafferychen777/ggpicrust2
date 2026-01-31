# Shared utilities for data-raw/ build scripts
#
# Usage: source("data-raw/build_utils.R") at the top of each build script.

library(httr)
library(jsonlite)

# =============================================================================
# Configuration
# =============================================================================

BUILD_CONFIG <- list(
  kegg_base     = "https://rest.kegg.jp",
  quickgo_base  = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms",
  output_dir    = "data",
  cache_dir     = "data-raw/.cache",
  http_timeout  = 30,
  max_retries   = 3
)

# =============================================================================
# HTTP
# =============================================================================

#' Robust HTTP GET with exponential backoff.
#'
#' @param url URL to fetch.
#' @param accept Optional Accept header (e.g. "application/json").
#' @param max_retries Maximum number of attempts.
#' @param as How to return: "response" (httr response), "text" (character), "parsed" (JSON).
#' @return Requested content or NULL on failure.
safe_get <- function(url,
                     accept = NULL,
                     max_retries = BUILD_CONFIG$max_retries,
                     as = "response") {
  for (attempt in seq_len(max_retries)) {
    resp <- tryCatch({
      headers <- if (!is.null(accept)) add_headers(Accept = accept)
      GET(url, timeout(BUILD_CONFIG$http_timeout), headers)
    }, error = function(e) NULL)

    if (!is.null(resp) && status_code(resp) == 200) {
      return(switch(as,
        "text"   = content(resp, "text", encoding = "UTF-8"),
        "parsed" = fromJSON(content(resp, "text", encoding = "UTF-8")),
        resp
      ))
    }

    if (attempt < max_retries) Sys.sleep(2^(attempt - 1))
  }
  NULL
}

# =============================================================================
# Progress reporting
# =============================================================================

#' Print a progress line (overwrites previous line on terminal).
show_progress <- function(current, total, detail = "") {
  pct <- sprintf("%.0f%%", current / total * 100)
  msg <- sprintf("  [%s] %d/%d", pct, current, total)
  if (nchar(detail) > 0) msg <- paste0(msg, " — ", detail)
  message(msg)
}

# =============================================================================
# Caching (for long-running API builds)
# =============================================================================

#' Get the path to a cache file.
cache_path <- function(name) {
  dir <- BUILD_CONFIG$cache_dir
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  file.path(dir, paste0(name, ".rds"))
}

#' Save intermediate results to cache.
cache_save <- function(obj, name) {
  saveRDS(obj, cache_path(name))
}

#' Load cached results. Returns NULL if not found or expired.
#' @param name Cache key.
#' @param max_age_hours Maximum age in hours (default 24).
cache_load <- function(name, max_age_hours = 24) {
  path <- cache_path(name)
  if (!file.exists(path)) return(NULL)

  info <- file.info(path)
  age_hours <- as.numeric(difftime(Sys.time(), info$mtime, units = "hours"))
  if (age_hours > max_age_hours) {
    message(sprintf("  Cache '%s' expired (%.1f hours old), re-fetching.", name, age_hours))
    return(NULL)
  }

  message(sprintf("  Using cached '%s' (%.1f hours old).", name, age_hours))
  readRDS(path)
}

#' Clear all cached data.
cache_clear <- function() {
  dir <- BUILD_CONFIG$cache_dir
  if (dir.exists(dir)) unlink(dir, recursive = TRUE)
  message("Cache cleared.")
}

# =============================================================================
# Validation helpers
# =============================================================================

#' Validate a reference data frame has the expected schema.
validate_reference <- function(df, required_cols, name = "reference") {
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing columns: %s", name, paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  if (nrow(df) == 0) {
    stop(sprintf("%s has zero rows.", name), call. = FALSE)
  }
  message(sprintf("  Validated %s: %d rows, columns [%s]",
                  name, nrow(df), paste(colnames(df), collapse = ", ")))
  invisible(TRUE)
}

#' Save a reference data frame to data/*.rda with consistent settings.
save_reference <- function(obj, obj_name) {
  path <- file.path(BUILD_CONFIG$output_dir, paste0(obj_name, ".rda"))
  env <- new.env(parent = emptyenv())
  env[[obj_name]] <- obj
  save(list = obj_name, file = path, envir = env, compress = "xz")
  size <- file.size(path)
  message(sprintf("  Saved %s (%s)", path, format_bytes(size)))
}

#' Human-readable byte size.
format_bytes <- function(bytes) {
  if (bytes < 1024) return(paste(bytes, "B"))
  if (bytes < 1024^2) return(sprintf("%.1f KB", bytes / 1024))
  sprintf("%.1f MB", bytes / 1024^2)
}
