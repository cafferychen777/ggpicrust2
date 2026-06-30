# =============================================================================
# Tests for taxa contribution analysis
# =============================================================================

# Shared test data generator
create_contrib_test_data <- function() {
  set.seed(123)
  samples <- c("S1", "S2")
  functions <- c("K00001", "K00002", "K00003")
  taxa <- c("ASV1", "ASV2", "ASV3", "ASV4")

  contrib <- expand.grid(
    sample = samples,
    function_id = functions,
    taxon = taxa,
    stringsAsFactors = FALSE
  )
  # Rename function_id -> function for raw input format
  contrib_raw <- contrib
  colnames(contrib_raw)[colnames(contrib_raw) == "function_id"] <- "function"
  contrib_raw$taxon_function_abun <- runif(nrow(contrib_raw))
  contrib_raw$norm_taxon_function_contrib <- runif(nrow(contrib_raw))
  contrib_raw$taxon_rel_function_abun <- contrib_raw$taxon_function_abun / 10

  path_contrib_raw <- contrib_raw
  path_contrib_raw$`function` <- rep(c("ko00010", "ko00020", "PWY-7219"), 8)
  path_contrib_raw$norm_taxon_function_contrib <- NULL

  real_path_contrib_raw <- data.frame(
    sample = c("100CHE6KO", "100CHE6KO", "101CHE6WT"),
    `function` = c("1CMET2-PWY", "ANAEROFRUCAT-PWY", "ANAGLYCOLYSIS-PWY"),
    taxon = c(
      "20e568023c10eaac834f1c110aacea18",
      "23fe12a325dfefcdb23447f43b6b896e",
      "288c8176059111c4c7fdfb0cd5afce64"
    ),
    taxon_abun = c(26, 108, 25.5),
    taxon_rel_abun = c(1.960784, 8.144796, 1.923077),
    genome_function_count = c(1.007415, 0.755562, 1.007416),
    taxon_function_abun = c(26.1928, 81.6007, 25.6891),
    taxon_rel_function_abun = c(1.9753, 6.1539, 1.9373),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Stratified format (wide)
  strat_long <- expand.grid(
    function_id = functions,
    taxon = taxa,
    stringsAsFactors = FALSE
  )
  strat_wide <- data.frame(
    `function` = strat_long$function_id,
    sequence = strat_long$taxon,
    S1 = runif(nrow(strat_long)),
    S2 = runif(nrow(strat_long)),
    check.names = FALSE
  )

  # Taxonomy (QIIME2 format)
  taxonomy_qiime2 <- data.frame(
    `Feature.ID` = taxa,
    Taxon = c(
      "k__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__acidophilus",
      "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__coli",
      "k__Bacteria;p__Firmicutes;c__Clostridia;o__Eubacteriales;f__Lachnospiraceae;g__Roseburia;s__intestinalis",
      "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__fragilis"
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Taxonomy (DADA2 format)
  taxonomy_dada2 <- data.frame(
    ASV = taxa,
    Kingdom = rep("Bacteria", 4),
    Phylum = c("Firmicutes", "Proteobacteria", "Firmicutes", "Bacteroidetes"),
    Class = c("Bacilli", "Gammaproteobacteria", "Clostridia", "Bacteroidia"),
    Order = c("Lactobacillales", "Enterobacterales", "Eubacteriales", "Bacteroidales"),
    Family = c("Lactobacillaceae", "Enterobacteriaceae", "Lachnospiraceae", "Bacteroidaceae"),
    Genus = c("Lactobacillus", "Escherichia", "Roseburia", "Bacteroides"),
    Species = c("acidophilus", "coli", "intestinalis", "fragilis"),
    stringsAsFactors = FALSE
  )

  # Metadata
  metadata <- data.frame(
    sample = c("S1", "S2"),
    group = c("Control", "Treatment"),
    stringsAsFactors = FALSE
  )

  list(
    contrib_raw = contrib_raw,
    path_contrib_raw = path_contrib_raw,
    real_path_contrib_raw = real_path_contrib_raw,
    strat_wide = strat_wide,
    taxonomy_qiime2 = taxonomy_qiime2,
    taxonomy_dada2 = taxonomy_dada2,
    metadata = metadata
  )
}


# ---- read_contrib_file ----

test_that("read_contrib_file parses contrib data correctly", {
  td <- create_contrib_test_data()
  result <- read_contrib_file(data = td$contrib_raw)

  expect_s3_class(result, "data.frame")
  expect_true("function_id" %in% colnames(result))
  expect_false("function" %in% colnames(result))
  expect_true(all(c("sample", "taxon", "taxon_function_abun",
                     "norm_taxon_function_contrib") %in% colnames(result)))
  expect_equal(nrow(result), nrow(td$contrib_raw))
})

test_that("read_contrib_file strips ko: prefix", {
  td <- create_contrib_test_data()
  td$contrib_raw$`function` <- paste0("ko:", td$contrib_raw$`function`)
  result <- read_contrib_file(data = td$contrib_raw)
  expect_false(any(grepl("^ko:", result$function_id)))
})

test_that("read_contrib_file accepts already standardized function_id", {
  td <- create_contrib_test_data()
  contrib_standardized <- td$contrib_raw
  colnames(contrib_standardized)[colnames(contrib_standardized) == "function"] <- "function_id"

  result <- read_contrib_file(data = contrib_standardized)

  expect_true("function_id" %in% colnames(result))
  expect_false("function" %in% colnames(result))
  expect_equal(nrow(result), nrow(contrib_standardized))
})

test_that("read_contrib_file accepts official-style contrib data without norm column", {
  td <- create_contrib_test_data()
  contrib_no_norm <- td$contrib_raw
  contrib_no_norm$norm_taxon_function_contrib <- NULL

  result <- read_contrib_file(data = contrib_no_norm)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("sample", "function_id", "taxon",
                    "taxon_function_abun", "taxon_rel_function_abun") %in%
                    colnames(result)))
  expect_false("norm_taxon_function_contrib" %in% colnames(result))
})

test_that("read_contrib_file rejects mixed gene-family and pathway identifiers", {
  td <- create_contrib_test_data()
  mixed_contrib <- td$contrib_raw
  mixed_contrib$`function`[1] <- "ko00010"

  expect_error(
    read_contrib_file(data = mixed_contrib),
    "mixes pathway-level and gene-family-level"
  )

  expect_error(
    read_contrib_file(data = td$contrib_raw, type = "pathway"),
    "declared as 'pathway'.*gene-family-level"
  )

  expect_error(
    read_contrib_file(data = td$path_contrib_raw, type = "gene_family"),
    "declared as 'gene_family'.*pathway-level"
  )
})

test_that("read_pathway_contrib_file parses path_abun_contrib-style data", {
  td <- create_contrib_test_data()

  result <- read_pathway_contrib_file(data = td$path_contrib_raw)

  expect_s3_class(result, "data.frame")
  expect_true("function_id" %in% colnames(result))
  expect_true(all(result$function_id %in% c("ko00010", "ko00020", "PWY-7219")))
  expect_true(all(result$feature_level == "pathway"))
})

test_that("read_pathway_contrib_file parses real PICRUSt2 pathway schema", {
  td <- create_contrib_test_data()

  result <- read_pathway_contrib_file(data = td$real_path_contrib_raw)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("sample", "function_id", "taxon", "taxon_abun",
                    "taxon_rel_abun", "genome_function_count",
                    "taxon_function_abun", "taxon_rel_function_abun") %in%
                    colnames(result)))
  expect_equal(result$function_id,
               c("1CMET2-PWY", "ANAEROFRUCAT-PWY", "ANAGLYCOLYSIS-PWY"))
  expect_true(all(result$feature_level == "pathway"))
})

test_that("read_pathway_contrib_file reads gzipped pathway contribution files", {
  td <- create_contrib_test_data()
  gz_file <- tempfile(fileext = ".tsv.gz")
  on.exit(unlink(gz_file), add = TRUE)

  con <- gzfile(gz_file, open = "wt")
  utils::write.table(
    td$real_path_contrib_raw,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  close(con)

  result <- read_pathway_contrib_file(gz_file)

  expect_equal(nrow(result), nrow(td$real_path_contrib_raw))
  expect_equal(result$function_id,
               c("1CMET2-PWY", "ANAEROFRUCAT-PWY", "ANAGLYCOLYSIS-PWY"))
  expect_true(all(result$feature_level == "pathway"))
})

test_that("real pathway contribution schema aggregates and annotates MetaCyc IDs", {
  td <- create_contrib_test_data()
  contrib <- read_pathway_contrib_file(data = td$real_path_contrib_raw)

  agg <- aggregate_taxa_contributions(contrib, top_n = 10)
  annotation <- suppressMessages(pathway_annotation(
    data = data.frame(function_id = unique(agg$function_id)),
    pathway = "MetaCyc"
  ))

  expect_gt(nrow(agg), 0)
  expect_true(all(c("sample", "function_id", "taxon_label", "contribution") %in%
                    colnames(agg)))
  expect_true(all(!is.na(agg$contribution)))
  expect_true(any(!is.na(annotation$description)))
})

test_that("read_contrib_file errors on missing columns", {
  bad_df <- data.frame(sample = "S1", x = 1)
  expect_error(read_contrib_file(data = bad_df), "Missing required columns")
})

test_that("read_contrib_file rejects missing contribution keys", {
  td <- create_contrib_test_data()

  bad_sample <- td$contrib_raw
  bad_sample$sample[1] <- NA_character_
  expect_error(
    read_contrib_file(data = bad_sample),
    "sample.*non-empty"
  )

  bad_function <- td$contrib_raw
  bad_function$`function`[1] <- ""
  expect_error(
    read_contrib_file(data = bad_function),
    "function_id.*non-empty"
  )

  bad_taxon <- td$contrib_raw
  bad_taxon$taxon[1] <- "  "
  expect_error(
    read_contrib_file(data = bad_taxon),
    "taxon.*non-empty"
  )
})

test_that("read_contrib_file requires file or data", {
  expect_error(read_contrib_file(), "Please provide either")
})


# ---- read_strat_file ----

test_that("read_strat_file converts wide to long format", {
  td <- create_contrib_test_data()
  result <- read_strat_file(data = td$strat_wide)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("function_id", "taxon", "sample", "abundance")
                  %in% colnames(result)))
  # 3 functions x 4 taxa x 2 samples = 24 rows
  expect_equal(nrow(result), 24)
})

test_that("read_strat_file strips ko: prefix", {
  td <- create_contrib_test_data()
  td$strat_wide$`function` <- paste0("ko:", td$strat_wide$`function`)
  result <- read_strat_file(data = td$strat_wide)
  expect_false(any(grepl("^ko:", result$function_id)))
})

test_that("read_strat_file accepts already standardized column names", {
  td <- create_contrib_test_data()
  strat_standardized <- td$strat_wide
  colnames(strat_standardized)[colnames(strat_standardized) == "function"] <- "function_id"
  colnames(strat_standardized)[colnames(strat_standardized) == "sequence"] <- "taxon"

  result <- read_strat_file(data = strat_standardized)

  expect_true(all(c("function_id", "taxon", "sample", "abundance")
                  %in% colnames(result)))
  expect_equal(nrow(result), 24)
})

test_that("read_strat_file errors on too few columns", {
  bad_df <- data.frame(`function` = "K00001", check.names = FALSE)
  expect_error(read_strat_file(data = bad_df), "at least 3 columns")
})

test_that("read_strat_file rejects missing keys and invalid sample columns", {
  td <- create_contrib_test_data()

  bad_function <- td$strat_wide
  bad_function$`function`[1] <- NA_character_
  expect_error(
    read_strat_file(data = bad_function),
    "function_id.*non-empty"
  )

  bad_taxon <- td$strat_wide
  bad_taxon$sequence[1] <- ""
  expect_error(
    read_strat_file(data = bad_taxon),
    "taxon.*non-empty"
  )

  duplicate_samples <- data.frame(
    `function` = "K00001",
    sequence = "ASV1",
    S1 = 1,
    S1 = 2,
    check.names = FALSE
  )
  expect_error(
    read_strat_file(data = duplicate_samples),
    "sample column names are duplicated"
  )

  empty_sample <- data.frame(
    `function` = "K00001",
    sequence = "ASV1",
    S1 = 1,
    check.names = FALSE
  )
  colnames(empty_sample)[3] <- ""
  expect_error(
    read_strat_file(data = empty_sample),
    "sample column names must be non-empty"
  )
})


# ---- aggregate_taxa_contributions ----

test_that("aggregate_taxa_contributions returns expected shape", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)
  agg <- aggregate_taxa_contributions(contrib, top_n = 4)

  expect_s3_class(agg, "data.frame")
  expect_true(all(c("sample", "function_id", "taxon_label", "contribution")
                  %in% colnames(agg)))
  expect_gt(nrow(agg), 0)
})

test_that("aggregate_taxa_contributions with QIIME2 taxonomy uses genus labels", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)
  agg <- aggregate_taxa_contributions(
    contrib,
    taxonomy = td$taxonomy_qiime2,
    tax_level = "Genus",
    top_n = 10
  )

  expect_true(any(agg$taxon_label %in%
                    c("Lactobacillus", "Escherichia", "Roseburia", "Bacteroides")))
})

test_that("aggregate_taxa_contributions with DADA2 taxonomy uses genus labels", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)
  agg <- aggregate_taxa_contributions(
    contrib,
    taxonomy = td$taxonomy_dada2,
    tax_level = "Genus",
    top_n = 10
  )

  expect_true(any(agg$taxon_label %in%
                    c("Lactobacillus", "Escherichia", "Roseburia", "Bacteroides")))
})

test_that("aggregate_taxa_contributions does not duplicate rows for repeated taxonomy IDs", {
  contrib <- data.frame(
    sample = "S1",
    function_id = "K00001",
    taxon = "ASV1",
    taxon_function_abun = 1,
    norm_taxon_function_contrib = 1,
    stringsAsFactors = FALSE
  )
  taxonomy <- data.frame(
    Feature.ID = c("ASV1", "ASV1"),
    Taxon = c("k__Bacteria;g__Alpha", "k__Bacteria;g__Beta"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  expect_warning(
    agg <- aggregate_taxa_contributions(contrib, taxonomy = taxonomy, top_n = 10),
    "Duplicate taxon IDs found"
  )
  expect_equal(nrow(agg), 1)
  expect_equal(sum(agg$contribution), 1)
})

test_that("aggregate_taxa_contributions creates 'Other' category with top_n", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)
  agg <- aggregate_taxa_contributions(contrib, top_n = 2)

  expect_true("Other" %in% agg$taxon_label)
  # Should have at most 3 unique labels (top 2 + Other)
  expect_lte(length(unique(agg$taxon_label)), 3)
})

test_that("aggregate_taxa_contributions ranks top taxa by total contribution", {
  contrib <- data.frame(
    sample = c("S1", "S1", "S2", "S2", "S1"),
    function_id = c("K00001", "K00002", "K00001", "K00002", "K00001"),
    taxon = c("Common", "Common", "Common", "Common", "Rare"),
    taxon_function_abun = c(4, 4, 4, 4, 10),
    stringsAsFactors = FALSE
  )

  agg <- aggregate_taxa_contributions(
    contrib,
    top_n = 1,
    contribution_col = "taxon_function_abun"
  )

  expect_true("Common" %in% agg$taxon_label)
  expect_false("Rare" %in% agg$taxon_label)
  expect_true("Other" %in% agg$taxon_label)
  expect_equal(sum(agg$contribution), sum(contrib$taxon_function_abun))
})

test_that("aggregate_taxa_contributions works with strat file input", {
  td <- create_contrib_test_data()
  strat <- read_strat_file(data = td$strat_wide)
  agg <- aggregate_taxa_contributions(strat, top_n = 4)

  expect_s3_class(agg, "data.frame")
  expect_true(all(c("sample", "function_id", "taxon_label", "contribution")
                  %in% colnames(agg)))
})

test_that("aggregate_taxa_contributions auto-selects absolute contribution when norm column is absent", {
  td <- create_contrib_test_data()
  contrib_no_norm <- read_contrib_file(data = {
    x <- td$contrib_raw
    x$norm_taxon_function_contrib <- NULL
    x
  })

  agg <- aggregate_taxa_contributions(contrib_no_norm, top_n = 4)

  expected_total <- sum(contrib_no_norm$taxon_function_abun)
  expect_equal(sum(agg$contribution), expected_total)
})

test_that("aggregate_taxa_contributions supports explicit contribution column", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)

  agg <- aggregate_taxa_contributions(
    contrib,
    top_n = 4,
    contribution_col = "taxon_rel_function_abun"
  )

  expect_equal(sum(agg$contribution), sum(contrib$taxon_rel_function_abun))
})

test_that("aggregate_taxa_contributions supports non-syntactic contribution columns", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)
  contrib[["taxon function contribution"]] <- contrib$taxon_function_abun

  agg <- aggregate_taxa_contributions(
    contrib,
    top_n = 4,
    contribution_col = "taxon function contribution"
  )

  expect_equal(sum(agg$contribution), sum(contrib[["taxon function contribution"]]))
})

test_that("aggregate_taxa_contributions validates p-values and contribution values", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)

  expect_error(
    aggregate_taxa_contributions(contrib, top_n = 0),
    "top_n.*positive"
  )

  daa_missing <- data.frame(feature = "K00001", stringsAsFactors = FALSE)
  expect_error(
    aggregate_taxa_contributions(contrib, daa_results_df = daa_missing),
    "Missing required columns"
  )

  daa_bad <- data.frame(
    feature = "K00001",
    p_adjust = 1.2,
    stringsAsFactors = FALSE
  )
  expect_error(
    aggregate_taxa_contributions(contrib, daa_results_df = daa_bad),
    "between 0 and 1"
  )
  expect_error(
    aggregate_taxa_contributions(contrib, daa_results_df = daa_bad,
                                 p_threshold = 0),
    "range"
  )

  contrib_na <- contrib
  contrib_na$norm_taxon_function_contrib[1] <- NA_real_
  expect_error(
    aggregate_taxa_contributions(contrib_na),
    "finite non-negative"
  )

  contrib_neg <- contrib
  contrib_neg$norm_taxon_function_contrib[1] <- -0.1
  expect_error(
    aggregate_taxa_contributions(contrib_neg),
    "finite non-negative"
  )

  contrib_bad_key <- contrib
  contrib_bad_key$sample[1] <- NA_character_
  expect_error(
    aggregate_taxa_contributions(contrib_bad_key),
    "sample.*non-empty"
  )

  daa_bad_feature <- data.frame(
    feature = NA_character_,
    p_adjust = 0.01,
    stringsAsFactors = FALSE
  )
  expect_error(
    aggregate_taxa_contributions(contrib, daa_results_df = daa_bad_feature),
    "feature.*non-empty"
  )

  expect_error(
    aggregate_taxa_contributions(contrib, pathway_ids = ""),
    "features.*non-empty"
  )
})

test_that("aggregate_taxa_contributions directly filters pathway-level contribution IDs", {
  td <- create_contrib_test_data()
  contrib <- read_pathway_contrib_file(data = td$path_contrib_raw)

  agg <- aggregate_taxa_contributions(
    contrib,
    pathway_ids = "ko00010",
    top_n = 4
  )

  expect_true(all(agg$function_id == "ko00010"))
  expect_gt(nrow(agg), 0)
})

test_that("aggregate_taxa_contributions maps KEGG pathways to KOs only for KO-level contribution data", {
  ko_to_kegg_reference <- ggpicrust2:::load_reference_data("ko_to_kegg")
  ko00010_kos <- head(unique(ko_to_kegg_reference$ko_id[
    ko_to_kegg_reference$pathway_id == "ko00010"
  ]), 2)
  skip_if(length(ko00010_kos) < 1)

  contrib <- data.frame(
    sample = "S1",
    `function` = ko00010_kos,
    taxon = paste0("ASV", seq_along(ko00010_kos)),
    taxon_function_abun = seq_along(ko00010_kos),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  contrib <- read_contrib_file(data = contrib)

  agg <- aggregate_taxa_contributions(contrib, pathway_ids = "ko00010", top_n = 10)

  expect_setequal(agg$function_id, ko00010_kos)
})

test_that("aggregate_taxa_contributions rejects mixed contribution feature levels", {
  contrib <- data.frame(
    sample = c("S1", "S1"),
    function_id = c("K00001", "ko00010"),
    taxon = c("ASV1", "ASV2"),
    taxon_function_abun = c(1, 2),
    stringsAsFactors = FALSE
  )

  expect_error(
    aggregate_taxa_contributions(contrib),
    "mixes pathway-level and gene-family-level"
  )

  contrib$feature_level <- c("gene_family", "pathway")
  expect_error(
    aggregate_taxa_contributions(contrib),
    "mixes contribution feature levels"
  )
})

test_that("aggregate_taxa_contributions with daa_results_df filters pathways", {
  skip_if_not_installed("ggpicrust2")

  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)

  # Create mock DAA results with a pathway that maps to our KOs
  daa_results <- data.frame(
    feature = "ko00010",
    p_adjust = 0.01,
    stringsAsFactors = FALSE
  )

  # This may error if ko00010 KOs don't match our test KOs, which is expected
  # The important thing is the pathway→KO mapping logic works
  result <- tryCatch(
    aggregate_taxa_contributions(contrib, daa_results_df = daa_results),
    error = function(e) e
  )

  # Should either succeed or fail with informative message about no matching data
  if (inherits(result, "error")) {
    expect_match(result$message, "No data remaining|No significant")
  } else {
    expect_s3_class(result, "data.frame")
  }
})


# ---- taxa_contribution_bar ----

test_that("taxa_contribution_bar returns ggplot", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)
  agg <- aggregate_taxa_contributions(contrib, top_n = 4)

  p <- taxa_contribution_bar(agg, td$metadata, group = "group")
  expect_s3_class(p, "ggplot")
})

test_that("taxa contribution visualizations validate contribution values and counts", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)
  agg <- aggregate_taxa_contributions(contrib, top_n = 4)

  expect_error(
    taxa_contribution_bar(agg, td$metadata, group = "group",
                          n_functions = 0),
    "n_functions.*positive"
  )
  expect_error(
    taxa_contribution_heatmap(agg, n_functions = 0),
    "n_functions.*positive"
  )

  agg_bad <- agg
  agg_bad$contribution[1] <- NA_real_
  expect_error(
    taxa_contribution_bar(agg_bad, td$metadata, group = "group"),
    "finite non-negative"
  )
  expect_error(
    taxa_contribution_heatmap(agg_bad),
    "finite non-negative"
  )

  agg_bad_key <- agg
  agg_bad_key$taxon_label[1] <- NA_character_
  expect_error(
    taxa_contribution_bar(agg_bad_key, td$metadata, group = "group"),
    "taxon_label.*non-empty"
  )
  expect_error(
    taxa_contribution_heatmap(agg_bad_key),
    "taxon_label.*non-empty"
  )

  expect_error(
    taxa_contribution_bar(agg, td$metadata, group = "group",
                          function_ids = ""),
    "function_ids.*non-empty"
  )
  expect_error(
    taxa_contribution_bar(agg, td$metadata, group = "group",
                          function_ids = "not-present"),
    "No contribution rows match"
  )
})

test_that("taxa_contribution_bar accepts facet_by = 'group'", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)
  agg <- aggregate_taxa_contributions(contrib, top_n = 4)

  p <- taxa_contribution_bar(agg, td$metadata, group = "group",
                             facet_by = "group")
  expect_s3_class(p, "ggplot")
})

test_that("taxa_contribution_bar supports metadata with a single group", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)
  agg <- aggregate_taxa_contributions(contrib, top_n = 4)
  td$metadata$group <- "Control"

  p <- taxa_contribution_bar(agg, td$metadata, group = "group")
  expect_s3_class(p, "ggplot")
})

test_that("taxa_contribution_bar ranks default functions by sample-level variance", {
  agg <- data.frame(
    sample = c("S1", "S1", "S2", "S2", "S1", "S2", "S2"),
    function_id = c(
      "F_variable", "F_variable", "F_variable", "F_variable",
      "F_constant", "F_constant", "F_constant"
    ),
    taxon_label = c("Taxon1", "Taxon2", "Taxon1", "Taxon2",
                    "Taxon1", "Taxon1", "Taxon2"),
    contribution = c(5, 5, 6, 6, 20, 10, 10),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(
    sample = c("S1", "S2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )

  p <- taxa_contribution_bar(
    agg,
    metadata,
    group = "group",
    n_functions = 1,
    show_percentage = FALSE
  )

  expect_setequal(unique(p$data$function_id), "F_variable")
})

test_that("taxa_contribution_bar rejects zero-total percentage denominators", {
  agg <- data.frame(
    sample = c("S1", "S1"),
    function_id = c("K00001", "K00001"),
    taxon_label = c("Taxon1", "Taxon2"),
    contribution = c(0, 0),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(
    sample = "S1",
    group = "A",
    stringsAsFactors = FALSE
  )

  expect_error(
    taxa_contribution_bar(agg, metadata, group = "group",
                          show_percentage = TRUE),
    "total contribution <= 0"
  )

  p <- taxa_contribution_bar(agg, metadata, group = "group",
                             show_percentage = FALSE)
  expect_s3_class(p, "ggplot")
  expect_equal(sum(p$data$contribution), 0)
})


# ---- taxa_contribution_heatmap ----

test_that("taxa_contribution_heatmap returns ggplot", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)
  agg <- aggregate_taxa_contributions(contrib, top_n = 4)

  p <- taxa_contribution_heatmap(agg, cluster_rows = FALSE, cluster_cols = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("taxa_contribution_heatmap treats absent sparse contributions as zero", {
  agg <- data.frame(
    sample = c("S1", "S1", "S2"),
    function_id = c("K00001", "K00001", "K00001"),
    taxon_label = c("Rare", "Common", "Common"),
    contribution = c(10, 4, 4),
    stringsAsFactors = FALSE
  )

  p <- taxa_contribution_heatmap(
    agg,
    n_functions = 1,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )

  plot_df <- p$data
  expect_equal(
    plot_df$value[as.character(plot_df$taxon) == "Rare"],
    5
  )
  expect_equal(
    plot_df$value[as.character(plot_df$taxon) == "Common"],
    4
  )
})

test_that("taxa_contribution_heatmap relabels columns using pathway_annotation output", {
  # Regression: the heatmap looked for an `annotation_data$pathway` ID
  # column, but pathway_annotation() emits the ID in `feature` (with
  # `description` as the human-readable label). Under the old code, the
  # left-join returned all-NA names and the axis silently kept the raw
  # function IDs. Accept both `feature`/`description` (non-ko_to_kegg)
  # and `pathway`/`pathway_name` (ko_to_kegg) column naming.
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)
  agg <- aggregate_taxa_contributions(contrib, top_n = 4)

  ann <- data.frame(
    feature = c("K00001", "K00002", "K00003"),
    description = c("alpha desc", "beta desc", "gamma desc"),
    stringsAsFactors = FALSE
  )
  p <- taxa_contribution_heatmap(
    agg,
    annotation_data = ann,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
  expect_s3_class(p, "ggplot")
  func_levels <- levels(p$data$func)
  # At least one of the labels should be the pathway_annotation description.
  expect_true(any(c("alpha desc", "beta desc", "gamma desc") %in% func_levels))
})

test_that("taxa_contribution_heatmap rejects conflicting annotation labels", {
  agg <- data.frame(
    sample = c("S1", "S1", "S2", "S2"),
    function_id = c("K00001", "K00002", "K00001", "K00002"),
    taxon_label = c("Taxon1", "Taxon1", "Taxon2", "Taxon2"),
    contribution = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )

  same_label_annotation <- data.frame(
    feature = c("K00001", "K00001"),
    description = c("shared label", "shared label"),
    stringsAsFactors = FALSE
  )
  p <- taxa_contribution_heatmap(
    agg,
    annotation_data = same_label_annotation,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
  expect_s3_class(p, "ggplot")
  expect_true("shared label" %in% levels(p$data$func))

  conflicting_annotation <- data.frame(
    feature = c("K00001", "K00001"),
    description = c("first label", "second label"),
    stringsAsFactors = FALSE
  )
  expect_error(
    taxa_contribution_heatmap(
      agg,
      annotation_data = conflicting_annotation,
      cluster_rows = FALSE,
      cluster_cols = FALSE
    ),
    "multiple labels"
  )
})

test_that("taxa_contribution_heatmap with clustering returns patchwork", {
  skip_if_not_installed("ggdendro")

  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)
  agg <- aggregate_taxa_contributions(contrib, top_n = 4)

  p <- taxa_contribution_heatmap(agg, cluster_rows = TRUE, cluster_cols = TRUE)
  # With dendrograms, patchwork wraps the result
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})


# ---- Edge cases ----

test_that("aggregate_taxa_contributions errors on empty data after filtering", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)

  # Use a syntactically valid but non-existent pathway ID so the filter
  # reaches the "no matching KOs in contrib_data" branch rather than the
  # strict ID-shape check.
  expect_error(
    aggregate_taxa_contributions(contrib, pathway_ids = "ko99999"),
    "No data remaining"
  )
})

test_that("aggregate_taxa_contributions errors with invalid tax_level", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)

  expect_error(
    aggregate_taxa_contributions(contrib, taxonomy = td$taxonomy_qiime2,
                                 tax_level = "Superclass"),
    "Invalid tax_level"
  )
})
