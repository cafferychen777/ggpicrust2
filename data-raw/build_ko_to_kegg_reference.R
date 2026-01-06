#!/usr/bin/env Rscript
# Create ko_to_kegg_reference dataset from user-provided KEGG BRITE hierarchy
# Data source: @KitHubb from Discussion #113
# Original data: https://www.genome.jp/kegg-bin/get_htext?htext=KO
# Date: 2024-08-12 (user's data)
# Purpose: Replace old wide-format (306x326, 84% NA) with new long-format

cat("=== Creating ko_to_kegg_reference dataset ===\n\n")

library(dplyr)
library(stringr)
library(readr)

# Download user's processed data if not exists
if (!file.exists("data-raw/KO_20240812_ff.csv")) {
  cat("Downloading user's data from GitHub...\n")
  download.file(
    url = "https://github.com/user-attachments/files/16577153/KO_20240812_ff.csv",
    destfile = "data-raw/KO_20240812_ff.csv",
    method = "auto"
  )
  cat("Download complete.\n\n")
}

# Load user's data
cat("Loading user's data...\n")
raw_data <- read_csv("data-raw/KO_20240812_ff.csv", show_col_types = FALSE)

cat(sprintf("Loaded: %d rows × %d columns\n", nrow(raw_data), ncol(raw_data)))
cat("Columns:", paste(colnames(raw_data), collapse = ", "), "\n\n")

# Transform to package internal format
cat("Transforming to package format...\n")

ko_to_kegg_reference <- raw_data %>%
  # Filter out empty entries
  filter(
    !is.na(KEGG) & KEGG != "",
    !is.na(KO_num) & KO_num != ""
  ) %>%
  # Create clean columns
  transmute(
    # Pathway information
    pathway_id = str_extract(KEGG, "^ko\\d+"),  # Extract only the pathway ID (ko followed by digits)
    pathway_number = str_extract(Level3, "^\\d{4}"),
    pathway_name = str_replace(Level3, "^\\d{4}\\s+(.+?)\\s*\\[PATH:.*$", "\\1") %>%
      trimws(),

    # KO information
    ko_id = str_trim(KO_num),  # Remove any whitespace
    ko_description = KO,

    # EC number
    ec_number = ifelse(is.na(EC) | EC == "", "", EC),

    # Hierarchical classification
    level1 = Level1,
    level2 = Level2,
    level3 = Level3
  ) %>%
  # Remove duplicates
  distinct()

cat(sprintf("After transformation: %d rows × %d columns\n\n",
            nrow(ko_to_kegg_reference), ncol(ko_to_kegg_reference)))

# Data quality checks
cat("=== Data Quality Checks ===\n\n")

check_results <- list()

# Check 1: No NA in key columns
check_results$no_na_pathway <- sum(is.na(ko_to_kegg_reference$pathway_id)) == 0
check_results$no_na_ko <- sum(is.na(ko_to_kegg_reference$ko_id)) == 0

cat(sprintf("✓ pathway_id has no NA: %s\n",
            ifelse(check_results$no_na_pathway, "PASS", "FAIL")))
cat(sprintf("✓ ko_id has no NA: %s\n",
            ifelse(check_results$no_na_ko, "PASS", "FAIL")))

# Check 2: ID format validation
valid_pathway_format <- all(grepl("^ko\\d{5}$", ko_to_kegg_reference$pathway_id))
valid_ko_format <- all(grepl("^K\\d{5}\\s*$", ko_to_kegg_reference$ko_id))

cat(sprintf("✓ pathway_id format (ko#####): %s\n",
            ifelse(valid_pathway_format, "PASS", "FAIL")))
cat(sprintf("✓ ko_id format (K#####): %s\n",
            ifelse(valid_ko_format, "PASS", "FAIL")))

# Check 3: Coverage statistics
unique_pathways <- length(unique(ko_to_kegg_reference$pathway_id))
unique_kos <- length(unique(ko_to_kegg_reference$ko_id))
total_mappings <- nrow(ko_to_kegg_reference)

check_results$sufficient_pathways <- unique_pathways >= 500
check_results$sufficient_kos <- unique_kos >= 25000
check_results$sufficient_mappings <- total_mappings >= 60000

cat(sprintf("\n✓ Unique pathways: %s (expected >= 500)\n",
            format(unique_pathways, big.mark = ",")))
cat(sprintf("✓ Unique KOs: %s (expected >= 25,000)\n",
            format(unique_kos, big.mark = ",")))
cat(sprintf("✓ Total mappings: %s (expected >= 60,000)\n",
            format(total_mappings, big.mark = ",")))

# Check 4: Compare with existing data (if available)
old_data_file <- "data/ko_to_kegg_reference.rda"
if (file.exists(old_data_file)) {
  old_env <- new.env()
  load(old_data_file, envir = old_env)
  old_data <- old_env$ko_to_kegg_reference
  old_mappings <- nrow(old_data)
  change_pct <- (total_mappings - old_mappings) / old_mappings * 100
  cat(sprintf("\n✓ Previous mappings: %s\n", format(old_mappings, big.mark = ",")))
  cat(sprintf("✓ New mappings: %s\n", format(total_mappings, big.mark = ",")))
  cat(sprintf("✓ Change: %+.1f%%\n", change_pct))
} else {
  cat("\n✓ No previous data to compare (new installation)\n")
}

# Check 5: Sample verification
cat("\n=== Sample Data Verification ===\n\n")
sample_data <- head(ko_to_kegg_reference, 3)
print(sample_data[, c("pathway_id", "ko_id", "pathway_name")])

# Check 6: Pathway size distribution
pathway_sizes <- ko_to_kegg_reference %>%
  group_by(pathway_id) %>%
  summarise(ko_count = n(), .groups = "drop")

cat("\n=== Pathway Size Distribution ===\n\n")
cat(sprintf("Min KOs per pathway: %d\n", min(pathway_sizes$ko_count)))
cat(sprintf("Max KOs per pathway: %d\n", max(pathway_sizes$ko_count)))
cat(sprintf("Mean KOs per pathway: %.1f\n", mean(pathway_sizes$ko_count)))
cat(sprintf("Median KOs per pathway: %.0f\n", median(pathway_sizes$ko_count)))

# Overall check
all_checks_pass <- all(unlist(check_results))

if (all_checks_pass) {
  cat("\n=== ✓ ALL QUALITY CHECKS PASSED ===\n\n")
} else {
  cat("\n=== ✗ SOME QUALITY CHECKS FAILED ===\n\n")
  stop("Data quality checks failed. Please review the data.")
}

# Save to package data
cat("Saving to data/ko_to_kegg_reference.rda...\n")

# Use usethis if available, otherwise use base save
if (requireNamespace("usethis", quietly = TRUE)) {
  usethis::use_data(ko_to_kegg_reference, overwrite = TRUE, compress = "xz")
} else {
  save(ko_to_kegg_reference,
       file = "data/ko_to_kegg_reference.rda",
       compress = "xz")
}

# File size
data_size <- file.size("data/ko_to_kegg_reference.rda")

cat("\n=== ✓ Dataset Creation Complete ===\n\n")
cat("Summary:\n")
cat(sprintf("  • Total mappings: %s\n", format(total_mappings, big.mark = ",")))
cat(sprintf("  • Unique pathways: %s\n", format(unique_pathways, big.mark = ",")))
cat(sprintf("  • Unique KOs: %s\n", format(unique_kos, big.mark = ",")))
cat(sprintf("  • File size: %.1f KB\n", data_size / 1024))
cat("\nFile updated: data/ko_to_kegg_reference.rda\n")
