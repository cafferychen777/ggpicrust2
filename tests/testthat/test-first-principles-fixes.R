# Regression tests for 8 issues fixed in the first-principles audit pass.
# Each `test_that()` block corresponds to one issue in that audit and
# reproduces the original broken behavior to prove the fix sticks.

# -----------------------------------------------------------------------------
# Issue 1: align_samples() must reject duplicate sample IDs when the caller
# passes `sample_col` explicitly. Before the fix, match() silently collapsed
# duplicates to the first hit and dropped other rows without warning.
# -----------------------------------------------------------------------------
test_that("align_samples() rejects duplicate IDs under explicit sample_col", {
  align <- getFromNamespace("align_samples", "ggpicrust2")

  abundance <- matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    nrow = 3,
    dimnames = list(paste0("f", 1:3), c("S1", "S2", "S3"))
  )
  metadata <- data.frame(
    sample_id = c("S1", "S1", "S2"),   # duplicated "S1"
    group = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )

  expect_error(
    align(abundance, metadata, sample_col = "sample_id", verbose = FALSE),
    regexp = "S1"
  )
  expect_error(
    align(abundance, metadata, sample_col = "sample_id", verbose = FALSE),
    regexp = "duplicated"
  )
})

test_that("align_samples() still succeeds on unique explicit sample_col", {
  align <- getFromNamespace("align_samples", "ggpicrust2")

  abundance <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 2,
    dimnames = list(paste0("f", 1:2), c("S1", "S2", "S3"))
  )
  metadata <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    group = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )

  res <- align(abundance, metadata, sample_col = "sample_id", verbose = FALSE)
  expect_equal(res$n_samples, 3)
  expect_equal(colnames(res$abundance), c("S1", "S2", "S3"))
})


# -----------------------------------------------------------------------------
# Issue 2: aggregate_taxa_contributions() used to treat every DAA feature as
# a pathway ID, then remap to KOs. That silently zeroed-out KO-level DAA
# results. The fix dispatches on ID shape: KO -> direct match, pathway ->
# expand via ko_to_kegg.
# -----------------------------------------------------------------------------

make_contrib <- function() {
  data.frame(
    sample = rep(c("S1", "S2"), each = 4),
    function_id = rep(c("K00001", "K00002", "K00003", "K00004"), 2),
    taxon = rep(c("ASV1", "ASV2"), each = 2, times = 2),
    taxon_function_abun = seq_len(8),
    norm_taxon_function_contrib = seq(0.1, 0.8, length.out = 8),
    stringsAsFactors = FALSE
  )
}

test_that("aggregate_taxa_contributions() accepts KO-level DAA results directly", {
  contrib <- make_contrib()
  daa <- data.frame(
    feature = c("K00001", "K00002"),
    p_adjust = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )
  # Before the fix this call errored with "No data remaining after filtering".
  agg <- aggregate_taxa_contributions(contrib, daa_results_df = daa, top_n = 2)
  expect_true(nrow(agg) > 0)
  expect_true(all(agg$function_id %in% c("K00001", "K00002")))
})

test_that("aggregate_taxa_contributions() errors clearly on unrecognized feature IDs", {
  contrib <- make_contrib()
  daa <- data.frame(
    feature = c("GENE_X", "GENE_Y"),
    p_adjust = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )
  expect_error(
    aggregate_taxa_contributions(contrib, daa_results_df = daa),
    regexp = "KO|pathway"
  )
})

test_that("aggregate_taxa_contributions() still honors pathway-level features", {
  # Use a pathway ID that ko_to_kegg maps to at least one KO present in contrib.
  ref <- getFromNamespace("load_reference_data", "ggpicrust2")("ko_to_kegg")
  kos <- unique(ref$ko_id[ref$pathway_id == ref$pathway_id[1]])
  skip_if(length(kos) == 0, "ko_to_kegg reference unavailable")

  contrib <- data.frame(
    sample = rep(c("S1", "S2"), each = length(kos)),
    function_id = rep(kos, 2),
    taxon = rep(c("ASV1", "ASV2"), length.out = 2 * length(kos)),
    taxon_function_abun = seq_len(2 * length(kos)),
    norm_taxon_function_contrib = seq(0.01, 1, length.out = 2 * length(kos)),
    stringsAsFactors = FALSE
  )
  daa <- data.frame(
    feature = ref$pathway_id[1],
    p_adjust = 0.01,
    stringsAsFactors = FALSE
  )
  agg <- aggregate_taxa_contributions(contrib, daa_results_df = daa, top_n = 2)
  expect_true(nrow(agg) > 0)
})


# -----------------------------------------------------------------------------
# Issue 3: pathway_heatmap() with cluster_rows = TRUE crashed on constant
# rows because scale() turned them into NA, and dist()/hclust() then
# received NAs. The fix coerces post-scale NA back to 0.
# -----------------------------------------------------------------------------
test_that("pathway_heatmap() survives constant pathway rows with clustering", {
  abundance <- matrix(
    c(
      10, 20, 30, 40,   # varying
       5,  5,  5,  5,   # constant - would have NA'd after scale()
       8, 18, 28, 38    # varying
    ),
    nrow = 3, byrow = TRUE,
    dimnames = list(paste0("Pathway", 1:3), paste0("S", 1:4))
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  expect_silent_run <- tryCatch({
    p <- pathway_heatmap(
      abundance = abundance,
      metadata = metadata,
      group = "group",
      cluster_rows = TRUE
    )
    TRUE
  }, error = function(e) e)

  expect_true(isTRUE(expect_silent_run) || inherits(expect_silent_run, "ggplot"))
  # Check the data went through with zero (not NA) for the constant row.
  # Extract the long data directly for assertion.
  df <- p$data
  const_vals <- df$Value[df$rowname == "Pathway2"]
  expect_true(all(const_vals == 0))
  expect_false(any(is.na(df$Value)))
})


# -----------------------------------------------------------------------------
# Issue 4: compare_gsea_daa(plot_type = "scatter") crashed deep inside merge()
# when DAA results lacked log2_fold_change. The fix validates required
# columns up front with an actionable error.
# -----------------------------------------------------------------------------
test_that("compare_gsea_daa('scatter') errors when effect-size columns missing", {
  gsea <- data.frame(
    pathway_id = c("ko00010", "ko00020"),
    p.adjust = c(0.01, 0.2),
    stringsAsFactors = FALSE
  )
  daa <- data.frame(
    feature = c("ko00010", "ko00020"),
    p_adjust = c(0.02, 0.5),
    stringsAsFactors = FALSE
  )
  expect_error(
    compare_gsea_daa(gsea, daa, plot_type = "scatter"),
    regexp = "log2_fold_change|NES"
  )
})

test_that("compare_gsea_daa('venn') still works without effect-size columns", {
  gsea <- data.frame(
    pathway_id = c("ko00010", "ko00020"),
    p.adjust = c(0.01, 0.2),
    stringsAsFactors = FALSE
  )
  daa <- data.frame(
    feature = c("ko00010", "ko00030"),
    p_adjust = c(0.02, 0.01),
    stringsAsFactors = FALSE
  )
  res <- compare_gsea_daa(gsea, daa, plot_type = "venn")
  expect_equal(res$results$n_overlap, 1)
  expect_equal(res$results$overlap, "ko00010")
})

test_that("compare_gsea_daa() counts significant pathways with set semantics", {
  gsea <- data.frame(
    pathway_id = c("ko00010", "ko00020"),
    p.adjust = c(0.01, 0.2),
    stringsAsFactors = FALSE
  )
  daa <- data.frame(
    feature = c("ko00010", "ko00010", "ko00030"),
    p_adjust = c(0.02, 0.03, 0.01),
    stringsAsFactors = FALSE
  )

  res <- compare_gsea_daa(gsea, daa, plot_type = "venn")
  expect_equal(res$results$n_gsea_total, 1)
  expect_equal(res$results$n_daa_total, 2)
  expect_equal(res$results$n_overlap, 1)
  expect_equal(res$results$daa_only, "ko00030")
})

test_that("compare_gsea_daa() validates p-value inputs and keeps zero p-values plottable", {
  gsea <- data.frame(
    pathway_id = c("ko00010", "ko00020"),
    NES = c(1.5, -1.2),
    p.adjust = c(0, 0.2),
    group1 = c("Treatment", "Treatment"),
    group2 = c("Control", "Control"),
    stringsAsFactors = FALSE
  )
  daa <- data.frame(
    feature = c("ko00010", "ko00020"),
    log2_fold_change = c(2, -2),
    p_adjust = c(0.02, 0.5),
    group1 = c("Control", "Control"),
    group2 = c("Treatment", "Treatment"),
    stringsAsFactors = FALSE
  )

  res <- compare_gsea_daa(gsea, daa, plot_type = "scatter")
  expect_true(all(is.finite(res$plot$data$gsea_neg_log10_p_adjust)))

  expect_error(
    compare_gsea_daa(gsea, daa, p_threshold = 0),
    "range"
  )

  gsea_bad <- gsea
  gsea_bad$p.adjust[1] <- 1.2
  expect_error(
    compare_gsea_daa(gsea_bad, daa),
    "between 0 and 1"
  )
})

test_that("compare_gsea_daa('scatter') rejects duplicate effect-size rows", {
  gsea_dup <- data.frame(
    pathway_id = c("ko00010", "ko00010"),
    NES = c(1.5, 1.8),
    p.adjust = c(0.01, 0.02),
    group1 = c("Treatment", "Treatment"),
    group2 = c("Control", "Control"),
    stringsAsFactors = FALSE
  )
  daa <- data.frame(
    feature = "ko00010",
    log2_fold_change = 2,
    p_adjust = 0.02,
    group1 = "Control",
    group2 = "Treatment",
    stringsAsFactors = FALSE
  )

  expect_error(
    compare_gsea_daa(gsea_dup, daa, plot_type = "scatter"),
    "one effect-size row per pathway|Duplicate pathway_id"
  )

  gsea <- gsea_dup[1, , drop = FALSE]
  daa_dup <- data.frame(
    feature = c("ko00010", "ko00010"),
    log2_fold_change = c(2, 3),
    p_adjust = c(0.02, 0.03),
    group1 = c("Control", "Control"),
    group2 = c("Treatment", "Treatment"),
    stringsAsFactors = FALSE
  )
  expect_error(
    compare_gsea_daa(gsea, daa_dup, plot_type = "scatter"),
    "one effect-size row per pathway|Duplicate feature"
  )
})

test_that("compare_gsea_daa() rejects missing pathway identifiers", {
  gsea <- data.frame(
    pathway_id = c("ko00010", NA_character_),
    p.adjust = c(0.01, 0.2),
    stringsAsFactors = FALSE
  )
  daa <- data.frame(
    feature = c("ko00010", "ko00020"),
    p_adjust = c(0.02, 0.5),
    stringsAsFactors = FALSE
  )

  expect_error(
    compare_gsea_daa(gsea, daa, plot_type = "venn"),
    "pathway_id.*non-empty"
  )

  daa$feature[2] <- ""
  gsea$pathway_id[2] <- "ko00020"
  expect_error(
    compare_gsea_daa(gsea, daa, plot_type = "venn"),
    "feature.*non-empty"
  )
})


# -----------------------------------------------------------------------------
# Issue 5: compare_daa_results() used to double-filter diff_features so
# that a feature present in 2 of 3 methods was silently collapsed to 0
# disagreement. The fix defines diff = method_features \ intersection.
# -----------------------------------------------------------------------------
test_that("compare_daa_results() counts partial-agreement features as diffs", {
  # Feature A is in all three methods -> common to all.
  # Feature B is in methods 1 and 2 only -> disagreement from method 3's
  # perspective, and should appear in method 1 & 2 diff lists.
  m1 <- data.frame(feature = c("A", "B"),
                   group1 = "g1", group2 = "g2",
                   p_adjust = c(0.01, 0.02),
                   stringsAsFactors = FALSE)
  m2 <- data.frame(feature = c("A", "B"),
                   group1 = "g1", group2 = "g2",
                   p_adjust = c(0.01, 0.03),
                   stringsAsFactors = FALSE)
  m3 <- data.frame(feature = c("A"),
                   group1 = "g1", group2 = "g2",
                   p_adjust = c(0.04),
                   stringsAsFactors = FALSE)

  suppressMessages(
    out <- compare_daa_results(
      list(m1, m2, m3),
      method_names = c("m1", "m2", "m3"),
      p_values_threshold = 0.05
    )
  )

  # A is the single unanimous feature.
  expect_equal(unique(out$num_common_features), 1)

  # m1 and m2 each report B as a diff (not endorsed by m3).
  # m3 has no diffs because A is unanimous and B is absent.
  diffs <- setNames(out$num_diff_features, out$method)
  expect_equal(unname(diffs["m1"]), 1)
  expect_equal(unname(diffs["m2"]), 1)
  expect_equal(unname(diffs["m3"]), 0)

  # String column reflects the same thing.
  expect_equal(out$diff_features[out$method == "m1"], "B")
  expect_equal(out$diff_features[out$method == "m2"], "B")
})

test_that("compare_daa_results() validates adjusted p-values and threshold", {
  valid <- data.frame(
    feature = c("A", "B"),
    group1 = "g1",
    group2 = "g2",
    p_adjust = c(0.01, 0.2),
    stringsAsFactors = FALSE
  )

  missing_p <- valid[, setdiff(names(valid), "p_adjust"), drop = FALSE]
  expect_error(
    compare_daa_results(list(missing_p), "m1"),
    "Missing required columns"
  )

  invalid_p <- valid
  invalid_p$p_adjust[1] <- 1.2
  expect_error(
    compare_daa_results(list(invalid_p), "m1"),
    "between 0 and 1"
  )

  expect_error(
    compare_daa_results(list(valid), "m1", p_values_threshold = NA_real_),
    "single finite numeric"
  )
})

test_that("compare_daa_results() canonicalizes the legacy ALDEx2 KW label", {
  df <- data.frame(feature = "A",
                   group1 = "g1", group2 = "g2",
                   p_adjust = 0.01,
                   stringsAsFactors = FALSE)

  suppressMessages(
    out <- compare_daa_results(
      list(df),
      method_names = "ALDEx2_Kruskal-Wallace test",
      p_values_threshold = 0.05
    )
  )

  expect_equal(out$method, "ALDEx2_Kruskal-Wallis test")
})

test_that("compare_daa_results() compares feature/group-pair units", {
  m1 <- data.frame(
    feature = "A",
    group1 = "g1",
    group2 = "g2",
    p_adjust = 0.01,
    stringsAsFactors = FALSE
  )
  m2 <- data.frame(
    feature = "A",
    group1 = "g2",
    group2 = "g3",
    p_adjust = 0.01,
    stringsAsFactors = FALSE
  )

  suppressMessages(
    out <- compare_daa_results(
      list(m1, m2),
      method_names = c("m1", "m2"),
      p_values_threshold = 0.05
    )
  )

  expect_equal(unique(out$num_common_features), 0)
  expect_equal(out$num_diff_features, c(1L, 1L))
  expect_equal(out$diff_features[out$method == "m1"], "A [g1 vs g2]")
  expect_equal(out$diff_features[out$method == "m2"], "A [g2 vs g3]")
})

test_that("compare_daa_results() treats reversed two-group pairs as the same comparison", {
  m1 <- data.frame(
    feature = "A",
    group1 = "control",
    group2 = "treatment",
    p_adjust = 0.01,
    stringsAsFactors = FALSE
  )
  m2 <- data.frame(
    feature = "A",
    group1 = "treatment",
    group2 = "control",
    p_adjust = 0.02,
    stringsAsFactors = FALSE
  )

  suppressMessages(
    out <- compare_daa_results(
      list(m1, m2),
      method_names = c("m1", "m2"),
      p_values_threshold = 0.05
    )
  )

  expect_equal(unique(out$num_common_features), 1)
  expect_equal(out$num_diff_features, c(0L, 0L))
  expect_equal(out$common_features, c("A", "A"))
})

test_that("compare_daa_results() rejects ambiguous method names and missing IDs", {
  df <- data.frame(
    feature = "A",
    group1 = "g1",
    group2 = "g2",
    p_adjust = 0.01,
    stringsAsFactors = FALSE
  )

  expect_error(
    compare_daa_results(list(df, df), method_names = c("m1", "m1")),
    "duplicated"
  )

  df$feature <- NA_character_
  expect_error(
    compare_daa_results(list(df), method_names = "m1"),
    "feature.*non-empty"
  )
})


# -----------------------------------------------------------------------------
# Issue 6: visualize_gsea() heatmap branch used to rely solely on
# rownames(metadata). When metadata had a sample column + default
# integer rownames, column annotation came out as all-NA. The fix
# routes through align_samples() so it works the same as every other
# entry point in the package.
# -----------------------------------------------------------------------------
test_that("visualize_gsea('heatmap') aligns via sample column, not just rownames", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")

  gsea_results <- create_gsea_test_results(n_pathways = 5)
  gsea_results$leading_edge <- vapply(seq_len(nrow(gsea_results)), function(i) {
    paste(paste0("K", sprintf("%05d", ((i - 1) * 5 + 1):(i * 5))), collapse = ";")
  }, character(1))

  set.seed(42)
  abundance <- as.data.frame(matrix(
    runif(50 * 8, 1, 100), nrow = 50, ncol = 8,
    dimnames = list(paste0("K", sprintf("%05d", 1:50)), paste0("Sample", 1:8))
  ))

  # Metadata uses sample_name column; rownames are default integers.
  metadata <- data.frame(
    sample_name = paste0("Sample", 1:8),
    group = rep(c("G1", "G2"), each = 4),
    stringsAsFactors = FALSE
  )
  stopifnot(identical(rownames(metadata), as.character(1:8)))

  hm <- visualize_gsea(
    gsea_results = gsea_results,
    plot_type = "heatmap",
    abundance = abundance,
    metadata = metadata,
    group = "group"
  )
  expect_s4_class(hm, "Heatmap")

  # Extract the HeatmapAnnotation and confirm Group values are not all-NA.
  # Access through top_annotation slot for ComplexHeatmap objects.
  top_anno <- hm@top_annotation
  anno_df <- top_anno@anno_list[[1]]@color_mapping@levels
  expect_true(all(c("G1", "G2") %in% anno_df))
})


# -----------------------------------------------------------------------------
# Issue 7: visualize_gsea() unnecessarily required enrichplot for
# enrichment_plot/dotplot/barplot. Those branches are pure ggplot2.
# -----------------------------------------------------------------------------
test_that("visualize_gsea() pure-ggplot branches do not require enrichplot", {
  # We cannot actually uninstall enrichplot in the test harness, so
  # instead confirm the visualize_gsea() body no longer names enrichplot
  # for these plot types. This guards against accidental re-introduction
  # of the requireNamespace() check.
  body_src <- paste(deparse(body(ggpicrust2::visualize_gsea)), collapse = "\n")
  # The helpful inline comment still mentions enrichplot, but the
  # active check -- requireNamespace("enrichplot", ...) -- should be gone
  # for enrichment_plot/dotplot/barplot.
  expect_false(grepl("requireNamespace\\(\\s*\"enrichplot\"", body_src))
})


# -----------------------------------------------------------------------------
# Issue 8: run_limma_gsea() must not replace a failed voom fit with a log2
# transform and unit weights because that changes the statistical model.
# -----------------------------------------------------------------------------
test_that("run_limma_gsea() has no unit-weight voom fallback", {
  run_limma_gsea <- getFromNamespace("run_limma_gsea", "ggpicrust2")
  body_src <- paste(deparse(body(run_limma_gsea)), collapse = "\n")

  expect_false(grepl("v\\s*<<-", body_src))
  expect_false(grepl("class\\(v\\)\\s*<<-", body_src))
  expect_false(grepl("log2\\(abundance_mat \\+ 0\\.5\\)", body_src))
  expect_false(grepl("weights\\s*=\\s*matrix\\(1", body_src))
})

test_that("run_limma_gsea() passes raw counts to voom", {
  run_limma_gsea <- getFromNamespace("run_limma_gsea", "ggpicrust2")
  body_src <- paste(deparse(body(run_limma_gsea)), collapse = "\n")

  # voom() expects raw counts and applies its own 0.5 offset when computing
  # logCPM values. Adding a pseudocount before voom changes library sizes and
  # the estimated mean-variance relationship.
  expect_false(grepl("abundance_mat\\s*<-\\s*abundance_mat\\s*\\+\\s*0\\.5", body_src))
  expect_match(body_src, "limma::voom\\(abundance_mat, design")
})
