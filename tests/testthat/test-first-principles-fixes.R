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
# Issue 8: run_limma_gsea() voom fallback used <<- to inject v into the
# caller's frame. Verify the source no longer does that.
# -----------------------------------------------------------------------------
test_that("run_limma_gsea() no longer uses <<- in its voom fallback", {
  run_limma_gsea <- getFromNamespace("run_limma_gsea", "ggpicrust2")
  body_src <- paste(deparse(body(run_limma_gsea)), collapse = "\n")
  expect_false(grepl("v\\s*<<-", body_src))
  expect_false(grepl("class\\(v\\)\\s*<<-", body_src))
})
