# Build ko_to_kegg_reference from KEGG BRITE hierarchy
#
# Data source: KO_20240812_ff.csv (contributed by @KitHubb, GitHub Discussion #113)
# Original: https://www.genome.jp/kegg-bin/get_htext?htext=KO
#
# Usage:
#   source("data-raw/build_ko_to_kegg_reference.R")
#   build_ko_to_kegg_reference()

source("data-raw/build_utils.R")

library(dplyr)
library(stringr)
library(readr)

# =============================================================================
# Configuration
# =============================================================================

KEGG_CSV_URL  <- "https://github.com/user-attachments/files/16577153/KO_20240812_ff.csv"
KEGG_CSV_PATH <- "data-raw/KO_20240812_ff.csv"

# =============================================================================
# Build
# =============================================================================

#' Build ko_to_kegg_reference from pre-processed BRITE hierarchy CSV.
#'
#' The CSV contains the full KEGG BRITE ko00001 hierarchy, pre-parsed into
#' columns: KEGG, Level1, Level2, Level3, KO_num, KO, EC.
#'
#' @return data.frame: pathway_id, pathway_number, pathway_name, ko_id,
#'   ko_description, ec_number, level1, level2, level3 (invisible).
build_ko_to_kegg_reference <- function() {
  message("=== Building ko_to_kegg_reference ===")
  message("Source: KEGG BRITE hierarchy (KO_20240812_ff.csv)\n")

  # Download if not present
  if (!file.exists(KEGG_CSV_PATH)) {
    message("Downloading from GitHub...")
    download.file(KEGG_CSV_URL, KEGG_CSV_PATH, method = "auto")
  }

  # Load and transform
  raw <- read_csv(KEGG_CSV_PATH, show_col_types = FALSE)
  message(sprintf("Loaded: %d rows x %d columns", nrow(raw), ncol(raw)))

  ko_to_kegg_reference <- raw %>%
    filter(!is.na(KEGG) & KEGG != "", !is.na(KO_num) & KO_num != "") %>%
    transmute(
      pathway_id     = str_extract(KEGG, "^ko\\d+"),
      pathway_number = str_extract(Level3, "^\\d{4}"),
      pathway_name   = str_replace(Level3, "^\\d{4}\\s+(.+?)\\s*\\[PATH:.*$", "\\1") %>%
        trimws(),
      ko_id          = str_trim(KO_num),
      ko_description = KO,
      ec_number      = ifelse(is.na(EC) | EC == "", "", EC),
      level1         = Level1,
      level2         = Level2,
      level3         = Level3
    ) %>%
    # Drop BRITE classification nodes (4-digit IDs like ko9980);
    # keep only real pathway maps (5-digit IDs like ko00010)
    filter(grepl("^ko\\d{5}$", pathway_id)) %>%
    distinct()

  # Structural assertions
  stopifnot(
    all(!is.na(ko_to_kegg_reference$pathway_id)),
    all(!is.na(ko_to_kegg_reference$ko_id))
  )

  n_pathways <- length(unique(ko_to_kegg_reference$pathway_id))
  n_kos      <- length(unique(ko_to_kegg_reference$ko_id))

  message(sprintf("Mappings: %s  |  Pathways: %s  |  KOs: %s",
                  format(nrow(ko_to_kegg_reference), big.mark = ","),
                  format(n_pathways, big.mark = ","),
                  format(n_kos, big.mark = ",")))

  validate_reference(ko_to_kegg_reference,
                     c("pathway_id", "ko_id", "pathway_name"),
                     "ko_to_kegg_reference")
  save_reference(ko_to_kegg_reference, "ko_to_kegg_reference")

  invisible(ko_to_kegg_reference)
}
