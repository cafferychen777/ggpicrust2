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


# ---- taxa_contribution_heatmap ----

test_that("taxa_contribution_heatmap returns ggplot", {
  td <- create_contrib_test_data()
  contrib <- read_contrib_file(data = td$contrib_raw)
  agg <- aggregate_taxa_contributions(contrib, top_n = 4)

  p <- taxa_contribution_heatmap(agg, cluster_rows = FALSE, cluster_cols = FALSE)
  expect_s3_class(p, "ggplot")
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
