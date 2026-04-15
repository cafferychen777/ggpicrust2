# Helper: create standard DAA test data
# 3 pathways, configurable samples/groups. Values are deterministic (no randomness).
create_daa_test_data <- function(n_samples = 4, n_groups = 2) {
  # Pool of abundance values per pathway (cycled if n_samples > 6)
  pool <- list(
    c(10, 20, 15, 30, 35, 25),
    c(20, 30, 25, 40, 45, 35),
    c(30, 40, 35, 50, 55, 45)
  )
  abundance <- data.frame(
    lapply(seq_len(n_samples), function(i) {
      sapply(pool, function(p) p[((i - 1) %% length(p)) + 1])
    })
  )
  colnames(abundance) <- paste0("sample", seq_len(n_samples))
  rownames(abundance) <- paste0("pathway", 1:3)

  spg <- n_samples %/% n_groups
  groups <- rep(c("control", "treatment", "other")[seq_len(n_groups)], each = spg)
  if (length(groups) < n_samples) groups <- c(groups, rep("other", n_samples - length(groups)))
  metadata <- data.frame(
    sample = colnames(abundance),
    group = groups,
    stringsAsFactors = FALSE
  )
  list(abundance = abundance, metadata = metadata)
}

test_that("pathway_daa works with basic inputs", {
  td <- create_daa_test_data(n_samples = 6)

  result <- pathway_daa(td$abundance, td$metadata, "group", daa_method = "ALDEx2")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("feature", "method", "p_values", "adj_method", "p_adjust") %in% colnames(result)))
  expect_gte(nrow(result), nrow(td$abundance))
  expect_true(all(grepl("ALDEx2", result$method)))
  expect_true(all(!is.na(result$p_values)))
  expect_true(all(result$p_values >= 0 & result$p_values <= 1))
})

test_that("pathway_daa validates inputs correctly", {
  abundance <- data.frame(
    sample1 = c(10, 20, 30), sample2 = c(20, 30, 40),
    sample3 = c(30, 40, 50), sample4 = c(40, 50, 60),
    row.names = c("pathway1", "pathway2", "pathway3")
  )

  # Mismatched sample names
  expect_error(
    pathway_daa(abundance, data.frame(wrong = paste0("x", 1:4), group = rep(c("a","b"), each = 2)), "group"),
    "Cannot find matching sample identifiers between abundance and metadata"
  )

  # Single group
  expect_error(
    pathway_daa(abundance, data.frame(sample = paste0("sample", 1:4), group = rep("control", 4)), "group"),
    "At least 2 groups are required"
  )

  # Insufficient sample size
  small_abundance <- abundance[, 1:3, drop = FALSE]
  expect_error(
    pathway_daa(small_abundance, data.frame(sample = paste0("sample", 1:3), group = c("a","a","b")), "group"),
    "At least 4 samples are required"
  )
})

test_that("pathway_daa core methods produce expected results", {
  n_samples <- 10
  n_features <- 3

  set.seed(123)
  abundance <- as.data.frame(matrix(
    rpois(n_samples * n_features, lambda = 20),
    nrow = n_features, ncol = n_samples,
    dimnames = list(paste0("pathway", 1:n_features), paste0("sample", 1:n_samples))
  ))

  metadata <- data.frame(
    sample = paste0("sample", 1:n_samples),
    group = rep(c("control", "treatment"), each = n_samples / 2)
  )

  core_methods <- c("ALDEx2", "limma voom", "edgeR")
  for (method in core_methods) {
    result <- suppressWarnings(pathway_daa(abundance, metadata, "group", daa_method = method))

    expect_true(is.data.frame(result))
    expect_true(all(c("feature", "method", "p_values") %in% colnames(result)))

    if (method == "ALDEx2") {
      expect_gte(nrow(result), n_features)
    } else {
      expect_equal(nrow(result), n_features)
    }

    expect_true(all(result$p_values >= 0 & result$p_values <= 1))
  }
})

test_that("pathway_daa extended methods run when explicitly enabled", {
  skip_if(
    Sys.getenv("GGPICRUST2_RUN_EXTENDED_DAA_TESTS", "false") != "true",
    "Set GGPICRUST2_RUN_EXTENDED_DAA_TESTS=true to run extended DAA method tests."
  )

  n_samples <- 10
  n_features <- 3
  set.seed(123)
  abundance <- as.data.frame(matrix(
    rpois(n_samples * n_features, lambda = 20),
    nrow = n_features, ncol = n_samples,
    dimnames = list(paste0("pathway", 1:n_features), paste0("sample", 1:n_samples))
  ))
  metadata <- data.frame(
    sample = paste0("sample", 1:n_samples),
    group = rep(c("control", "treatment"), each = n_samples / 2)
  )

  method_pkg <- c(
    "DESeq2" = "DESeq2",
    "metagenomeSeq" = "metagenomeSeq",
    "Maaslin2" = "Maaslin2"
  )
  extended_methods <- c("DESeq2", "metagenomeSeq", "LinDA", "Maaslin2")

  for (method in extended_methods) {
    if (method %in% names(method_pkg)) {
      skip_if_not_installed(method_pkg[[method]])
    }

    # Capture noisy method output to keep default test logs readable.
    captured <- capture.output({
      result <- suppressWarnings(pathway_daa(abundance, metadata, "group", daa_method = method))
    }, type = "output")
    ignore <- captured

    expect_true(is.data.frame(result))
    expect_true(all(c("feature", "method", "p_values") %in% colnames(result)))
  }
})

test_that("pathway_daa handles sample selection correctly", {
  # Use 6 samples so selecting 4 still meets the minimum requirement
  td <- create_daa_test_data(n_samples = 6)

  # Select a true subset (4 of 6): 2 control + 2 treatment
  selected <- c("sample1", "sample2", "sample4", "sample5")
  result <- pathway_daa(td$abundance, td$metadata, "group",
                       daa_method = "ALDEx2", select = selected)
  expect_s3_class(result, "data.frame")

  # Invalid sample selection
  expect_error(
    pathway_daa(td$abundance, td$metadata, "group",
                daa_method = "ALDEx2",
                select = c("sample1", "invalid_sample")),
    "Some selected samples not in abundance data"
  )
})

test_that("pathway_daa select= keeps metadata rows aligned with abundance columns", {
  # Regression: previously the select= branch reordered abundance columns to
  # match `select`, but only filtered metadata rows without reordering them.
  # Group labels then drifted relative to abundance, producing wrong p-values
  # (and sometimes significance flips) when `select` was not in natural order.
  #
  # Build a strong-signal dataset and compare results against a pre-subset
  # baseline. Internal select= must equal the pre-subset path.
  set.seed(7)
  n_features <- 3
  ctl <- matrix(rpois(n_features * 3, lambda = 5),  nrow = n_features)
  trt <- matrix(rpois(n_features * 3, lambda = 50), nrow = n_features)
  abundance <- as.data.frame(cbind(ctl, trt))
  colnames(abundance) <- paste0("S", 1:6)
  rownames(abundance) <- paste0("p", seq_len(n_features))

  metadata <- data.frame(
    sample = colnames(abundance),
    group  = rep(c("ctl", "trt"), each = 3),
    stringsAsFactors = FALSE
  )

  # Pick samples in a deliberately non-natural order that mixes the two groups.
  selected <- c("S4", "S2", "S6", "S1")

  via_select <- pathway_daa(abundance, metadata, "group",
                            daa_method = "ALDEx2", select = selected)
  via_preset <- pathway_daa(abundance[, selected, drop = FALSE],
                            metadata[metadata$sample %in% selected, ],
                            "group",
                            daa_method = "ALDEx2")

  # Same features, same method rows, same (group1, group2) assignment.
  expect_equal(via_select$feature, via_preset$feature)
  expect_equal(via_select$method,  via_preset$method)
  expect_equal(via_select$group1,  via_preset$group1)
  expect_equal(via_select$group2,  via_preset$group2)

  # p-values should agree closely (ALDEx2 is Monte Carlo-based but
  # deterministic enough here for a tight tolerance).
  expect_true(cor(via_select$p_values, via_preset$p_values) > 0.99)
})

test_that("limma voom multi-group labels align with p-values and coefficients", {
  # Regression: group2 was assigned a length-(k-1) vector and R's recycling
  # produced interleaved B,C,B,C,... labels for a 3-group design. The correct
  # labeling is B repeated N_features times, then C repeated N_features times,
  # matching as.vector()'s column-major flattening of fit$p.value[,-1].
  n_features <- 4
  n_per_group <- 4
  set.seed(42)
  abundance <- as.data.frame(matrix(
    rpois(n_features * n_per_group * 3, lambda = 20),
    nrow = n_features,
    dimnames = list(paste0("p", seq_len(n_features)),
                    paste0("S", seq_len(n_per_group * 3)))
  ))
  metadata <- data.frame(
    sample = colnames(abundance),
    group  = rep(c("A", "B", "C"), each = n_per_group),
    stringsAsFactors = FALSE
  )

  res <- suppressWarnings(
    pathway_daa(abundance, metadata, "group",
                daa_method = "limma voom", reference = "A")
  )

  # 4 features * 2 non-reference contrasts = 8 rows.
  expect_equal(nrow(res), n_features * 2)

  # Within each contrast block, every feature must appear exactly once
  # (not interleaved across contrasts).
  expect_equal(sort(res$feature[res$group2 == "B"]), paste0("p", 1:n_features))
  expect_equal(sort(res$feature[res$group2 == "C"]), paste0("p", 1:n_features))

  # group1 must be the reference for every row.
  expect_true(all(res$group1 == "A"))
})

test_that("DESeq2 respects the user-supplied reference level", {
  skip_if_not_installed("DESeq2")
  td <- create_daa_test_data(n_samples = 6)

  captured <- capture.output({
    res_default   <- suppressWarnings(
      pathway_daa(td$abundance, td$metadata, "group", daa_method = "DESeq2")
    )
    res_reference <- suppressWarnings(
      pathway_daa(td$abundance, td$metadata, "group",
                  daa_method = "DESeq2", reference = "treatment")
    )
  }, type = "output")
  ignore <- captured

  # Default behavior unchanged: alphabetical first level is the reference.
  expect_true(all(res_default$group1 == "control"))
  expect_true(all(res_default$group2 == "treatment"))

  # Explicit reference inverts the contrast direction.
  expect_true(all(res_reference$group1 == "treatment"))
  expect_true(all(res_reference$group2 == "control"))

  # Flipping the reference flips the sign of log2 fold changes. DESeq2
  # re-fits dispersion per call so the magnitudes are not bit-identical,
  # but should match to several decimals.
  merged <- merge(
    res_default[, c("feature", "log2_fold_change")],
    res_reference[, c("feature", "log2_fold_change")],
    by = "feature", suffixes = c("_def", "_ref")
  )
  finite <- is.finite(merged$log2_fold_change_def) &
            is.finite(merged$log2_fold_change_ref)
  expect_true(all(abs(
    merged$log2_fold_change_def[finite] + merged$log2_fold_change_ref[finite]
  ) < 1e-4))
})

test_that("DESeq2 handles multi-group input", {
  skip_if_not_installed("DESeq2")

  td <- create_daa_test_data(n_samples = 9, n_groups = 3)

  captured <- capture.output({
    res <- suppressWarnings(
      pathway_daa(td$abundance, td$metadata, "group",
                  daa_method = "DESeq2", reference = "control")
    )
  }, type = "output")
  ignore <- captured

  # 3 features * 2 non-reference contrasts = 6 rows.
  expect_equal(nrow(res), nrow(td$abundance) * 2)
  expect_true(all(res$group1 == "control"))
  expect_setequal(unique(res$group2), c("treatment", "other"))
})

test_that("Maaslin2 multi-group produces one row per (feature, contrast)", {
  skip_if_not_installed("Maaslin2")

  n_features <- 4
  n_per_group <- 4
  set.seed(123)
  abundance <- as.data.frame(matrix(
    rpois(n_features * n_per_group * 3, lambda = 30),
    nrow = n_features,
    dimnames = list(paste0("p", seq_len(n_features)),
                    paste0("S", seq_len(n_per_group * 3)))
  ))
  metadata <- data.frame(
    sample = colnames(abundance),
    group  = rep(c("A", "B", "C"), each = n_per_group),
    stringsAsFactors = FALSE
  )

  captured <- capture.output({
    res <- suppressWarnings(
      pathway_daa(abundance, metadata, "group",
                  daa_method = "Maaslin2", reference = "A")
    )
  }, type = "output")
  ignore <- captured

  # Was flattened to n_features rows pre-fix; expect n_features * (k-1).
  expect_equal(nrow(res), n_features * 2)
  expect_true(all(res$group1 == "A"))
  expect_setequal(unique(res$group2), c("B", "C"))

  # Features preserved after Maaslin2's hyphen-to-dot renaming round-trip.
  expect_true(all(res$feature %in% rownames(abundance)))
})

test_that("metagenomeSeq works with non-default sample column name", {
  skip_if_not_installed("metagenomeSeq")

  # cumNormStatFast() needs enough inter-sample variance to estimate a
  # scaling quantile, so use Poisson-generated counts with n=10 samples
  # (matching the existing extended-method test pattern) rather than the
  # tiny deterministic pool from create_daa_test_data().
  n_samples  <- 10
  n_features <- 3
  set.seed(321)
  abundance <- as.data.frame(matrix(
    rpois(n_samples * n_features, lambda = 20),
    nrow = n_features,
    dimnames = list(paste0("p", seq_len(n_features)),
                    paste0("S", seq_len(n_samples)))
  ))
  # Exercise the branch that previously hardcoded metadata$sample by
  # naming the sample column something other than "sample".
  metadata <- data.frame(
    SampleID = colnames(abundance),
    group    = rep(c("control", "treatment"), each = n_samples / 2),
    stringsAsFactors = FALSE
  )

  captured <- capture.output({
    res <- suppressWarnings(
      pathway_daa(abundance, metadata, "group", daa_method = "metagenomeSeq")
    )
  }, type = "output")
  ignore <- captured

  expect_s3_class(res, "data.frame")
  expect_true(all(c("feature", "method", "p_values") %in% colnames(res)))
  expect_equal(nrow(res), n_features)
})

test_that("metagenomeSeq survives degenerate cumNormStat input", {
  skip_if_not_installed("metagenomeSeq")

  # Minimal-but-legal input: 4 samples, 3 features, monotonic counts.
  # metagenomeSeq::cumNormStatFast() returns NaN for this shape because the
  # per-sample quantile search has nothing to stabilize on, and the package
  # then aborts with the cryptic "missing value where TRUE/FALSE needed".
  # We pre-compute the normalization factor with a fallback to p = 0.5,
  # which is metagenomeSeq's own documented default.
  abundance <- data.frame(
    S1 = c(10, 20, 30), S2 = c(15, 25, 35),
    S3 = c(30, 40, 50), S4 = c(35, 45, 55),
    row.names = paste0("p", 1:3)
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group  = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  captured <- capture.output({
    res <- suppressWarnings(
      pathway_daa(abundance, metadata, "group", daa_method = "metagenomeSeq")
    )
  }, type = "output")
  ignore <- captured

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(abundance))
  expect_true(all(c("feature", "method", "p_values") %in% colnames(res)))
})

test_that("pathway_daa rejects negative abundance values", {
  # Regression: validate_daa_input() previously skipped numeric-matrix quality
  # checks, so negative entries silently flowed into downstream methods and
  # produced cryptic failures or nonsense results.
  abundance <- data.frame(
    S1 = c(10, 20, 30), S2 = c(20, -5, 40),
    S3 = c(30, 40, 50), S4 = c(40, 50, 60),
    row.names = paste0("p", 1:3)
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group  = c("a", "a", "b", "b"),
    stringsAsFactors = FALSE
  )
  expect_error(
    pathway_daa(abundance, metadata, "group", daa_method = "ALDEx2"),
    "Negative values found"
  )
})

test_that("pathway_daa handles factor levels correctly with subset", {
  # GitHub issue #158: 3 groups, select only 2
  abundance <- data.frame(
    sample1 = c(10, 20, 30), sample2 = c(20, 30, 40),
    sample3 = c(15, 25, 35), sample4 = c(30, 40, 50),
    sample5 = c(25, 35, 45),
    row.names = paste0("pathway", 1:3)
  )

  metadata <- data.frame(
    sample = paste0("sample", 1:5),
    group = c("A", "A", "B", "B", "C")
  )

  selected_samples <- c("sample1", "sample2", "sample3", "sample4")

  # Internal subsetting via select=
  set.seed(42)
  result1 <- pathway_daa(abundance, metadata, "group",
                        daa_method = "ALDEx2", select = selected_samples)

  # Pre-subsetting
  set.seed(42)
  result2 <- pathway_daa(abundance[, selected_samples, drop = FALSE],
                        metadata[metadata$sample %in% selected_samples, ],
                        "group", daa_method = "ALDEx2")

  expect_equal(nrow(result1), nrow(result2))
  expect_equal(result1$feature, result2$feature)
  expect_equal(result1$method, result2$method)

  # Groups should only include A and B, not C
  expect_true(all(result1$group1 %in% c("A", "B")))
  expect_true(all(result1$group2 %in% c("A", "B")))
  expect_false(any(c(result1$group1, result1$group2) == "C"))

  p_correlation <- cor(result1$p_values, result2$p_values, use = "complete.obs")
  expect_true(p_correlation > 0.9)
})

test_that("pathway_daa handles multiple groups correctly", {
  td <- create_daa_test_data(n_samples = 6, n_groups = 3)

  suppressWarnings({
    result <- pathway_daa(td$abundance, td$metadata, "group",
                         daa_method = "limma voom", reference = "control")
  })

  expect_s3_class(result, "data.frame")
  expect_true(all(!is.na(result$p_values)))

  # A correctly wired multi-group result should have one row per
  # (feature, non-reference level) pair -- `N_features * (k - 1)` total --
  # with `group1` pinned to the reference and `group2` covering every
  # other level. Asserting only `!is.na(p_values)` misses label/row-count
  # regressions like the one fixed in 2.5.15 where group2 was recycled
  # into interleaved B,C,B,C,... labels.
  expected_contrasts <- setdiff(unique(td$metadata$group), "control")
  expect_equal(nrow(result), nrow(td$abundance) * length(expected_contrasts))
  expect_true(all(result$group1 == "control"))
  expect_setequal(unique(result$group2), expected_contrasts)
  # Every feature must appear exactly once within each contrast block,
  # confirming `group2` labels align with per-feature p-values instead
  # of being recycled across rows.
  for (lvl in expected_contrasts) {
    expect_equal(sort(result$feature[result$group2 == lvl]),
                 sort(rownames(td$abundance)))
  }
})

test_that("pathway_daa handles p-value adjustment correctly", {
  td <- create_daa_test_data(n_samples = 4)

  # ALDEx2 uses its own pre-computed BH correction (Monte Carlo-based),
  # which is more accurate than simple p.adjust. The user's p_adjust_method
  # choice is not applied for ALDEx2. Other methods still honor the user's choice.
  result <- pathway_daa(td$abundance, td$metadata, "group",
                       daa_method = "ALDEx2", p_adjust_method = "BH")

  expect_true(all(result$adj_method == "BH (method-specific)"))
  expect_true(all(!is.na(result$p_adjust)))
  expect_true(all(result$p_adjust >= 0 & result$p_adjust <= 1))
})

test_that("pathway_daa include_abundance_stats parameter works correctly", {
  td <- create_daa_test_data(n_samples = 4)

  # Columns contributed by abundance stats (relative-abundance means/SDs).
  # log2_fold_change is intentionally excluded here: when the DAA method
  # already provides a log2_fold_change column (ALDEx2 with effect size on),
  # the relative-abundance ratio is suppressed to avoid conflating two
  # different effect-size definitions. That behavior is covered by a
  # dedicated test below.
  abundance_cols <- c("mean_rel_abundance_group1", "sd_rel_abundance_group1",
                     "mean_rel_abundance_group2", "sd_rel_abundance_group2")

  # Without abundance stats (default). Disable effect size too so this test
  # isolates the abundance-stats toggle from the ALDEx2 effect-size toggle.
  result_basic <- pathway_daa(td$abundance, td$metadata, "group",
                             daa_method = "ALDEx2",
                             include_abundance_stats = FALSE,
                             include_effect_size = FALSE)

  expect_false(any(abundance_cols %in% colnames(result_basic)))
  expect_false("log2_fold_change" %in% colnames(result_basic))

  # With abundance stats, still isolating from effect size.
  result_enhanced <- pathway_daa(td$abundance, td$metadata, "group",
                                daa_method = "ALDEx2",
                                include_abundance_stats = TRUE,
                                include_effect_size = FALSE)

  expect_true(all(abundance_cols %in% colnames(result_enhanced)))
  # With no method-native log2FC in the result, abundance stats should
  # contribute its relative-abundance-ratio log2_fold_change column.
  expect_true("log2_fold_change" %in% colnames(result_enhanced))

  for (col in c(abundance_cols, "log2_fold_change")) {
    expect_true(is.numeric(result_enhanced[[col]]))
  }
  expect_true(all(is.finite(result_enhanced$log2_fold_change)))
  expect_true(all(result_enhanced$sd_rel_abundance_group1 >= 0, na.rm = TRUE))
  expect_true(all(result_enhanced$sd_rel_abundance_group2 >= 0, na.rm = TRUE))
})

test_that("ALDEx2 returns effect size columns by default", {
  td <- create_daa_test_data(n_samples = 4)

  # Default call: no explicit include_effect_size.
  result_default <- pathway_daa(td$abundance, td$metadata, "group",
                                daa_method = "ALDEx2")

  effect_cols <- c("effect_size", "diff_btw", "log2_fold_change",
                   "rab_all", "overlap")
  expect_true(all(effect_cols %in% colnames(result_default)))
  for (col in effect_cols) {
    expect_true(is.numeric(result_default[[col]]))
  }

  # Opt-out recovers the minimal p-value-only output.
  result_opt_out <- pathway_daa(td$abundance, td$metadata, "group",
                                daa_method = "ALDEx2",
                                include_effect_size = FALSE)
  expect_false(any(effect_cols %in% colnames(result_opt_out)))
})

test_that("include_abundance_stats does not collide with method-native log2FC", {
  td <- create_daa_test_data(n_samples = 4)

  # ALDEx2 with both flags on: method-native log2_fold_change (CLR-space
  # diff.btw) must be preserved, and merging abundance stats must not
  # introduce .x / .y suffixes.
  result <- pathway_daa(td$abundance, td$metadata, "group",
                        daa_method = "ALDEx2",
                        include_abundance_stats = TRUE,
                        include_effect_size = TRUE)

  expect_true("log2_fold_change" %in% colnames(result))
  expect_false(any(c("log2_fold_change.x", "log2_fold_change.y") %in%
                     colnames(result)))
  # The retained log2_fold_change should equal ALDEx2's diff_btw (CLR space),
  # not the relative-abundance ratio.
  expect_equal(result$log2_fold_change, result$diff_btw)
})

test_that("format_linda_output handles malformed LinDA outputs robustly", {
  malformed_output <- list(
    groupB = data.frame(
      pvalue = I(matrix(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06), nrow = 3)),
      log2FoldChange = c(1.2, 0.8, 0.5),
      row.names = c("path1", "path2", "path3")
    ),
    groupC = data.frame(
      some_other_col = c(1, 2),
      row.names = c("path4", "path5")
    )
  )

  result <- suppressWarnings(
    ggpicrust2:::format_linda_output(
      linda_output = malformed_output,
      group = "group",
      reference = "A",
      Level = c("A", "B", "C")
    )
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
  expect_true(all(c("feature", "method", "group1", "group2", "p_values", "log2_fold_change") %in% colnames(result)))

  group_b_rows <- result$group2 == "B"
  expect_true(any(group_b_rows))
  expect_true(all(is.na(result$p_values[group_b_rows])))
  expect_true(all(!is.na(result$log2_fold_change[group_b_rows])))

  group_c_rows <- result$group2 == "C"
  expect_true(any(group_c_rows))
  expect_true(all(is.na(result$p_values[group_c_rows])))
  expect_true(all(is.na(result$log2_fold_change[group_c_rows])))
})

test_that("pathway_daa Lefser fails fast for multi-group input", {
  skip_if_not_installed("lefser")

  td <- create_daa_test_data(n_samples = 6, n_groups = 3)
  expect_error(
    pathway_daa(
      abundance = td$abundance,
      metadata = td$metadata,
      group = "group",
      daa_method = "Lefser"
    ),
    "requires exactly 2 groups"
  )
})

test_that("pathway_daa LinDA honors user-specified reference", {
  # Regression: LinDA used `~ group` without releveling, so the factor's
  # natural first level was always used as reference while the result's
  # group1 column was labeled with the user-supplied `reference`. That
  # produced rows with group1 == group2 and log2FC that did not reflect
  # the requested contrast direction.
  skip_if_not_installed("MicrobiomeStat")

  set.seed(42)
  n_feat <- 12; n_samp <- 12
  abund <- matrix(rpois(n_feat * n_samp, 40),
                  nrow = n_feat, ncol = n_samp,
                  dimnames = list(paste0("f", seq_len(n_feat)),
                                  paste0("S", seq_len(n_samp))))
  meta <- data.frame(
    sample = paste0("S", seq_len(n_samp)),
    Env = c(rep("control", 6), rep("treatment", 6)),
    stringsAsFactors = FALSE
  )

  r_ctrl <- suppressMessages(suppressWarnings(
    pathway_daa(abundance = abund, metadata = meta, group = "Env",
                daa_method = "LinDA", reference = "control")
  ))
  r_trt <- suppressMessages(suppressWarnings(
    pathway_daa(abundance = abund, metadata = meta, group = "Env",
                daa_method = "LinDA", reference = "treatment")
  ))

  expect_true(all(r_ctrl$group1 == "control"))
  expect_true(all(r_ctrl$group2 == "treatment"))
  expect_true(all(r_trt$group1 == "treatment"))
  expect_true(all(r_trt$group2 == "control"))
  # Flipping the reference must flip the sign of log2 fold change.
  expect_equal(r_ctrl$log2_fold_change, -r_trt$log2_fold_change, tolerance = 1e-6)
})

test_that("pathway_daa Maaslin2 honors user-specified reference in 2-group case", {
  # Regression: the 2-group branch passed `reference = NULL` to Maaslin2,
  # so Maaslin2 used its alphabetical default regardless of the user's
  # `reference` argument. Combined with the result-labeling that used
  # the user's `reference` for group1, this produced group1 == group2
  # and a coefficient whose sign did not flip when the user flipped
  # the reference.
  skip_if_not_installed("Maaslin2")

  set.seed(42)
  n_feat <- 12; n_samp <- 12
  abund <- matrix(rpois(n_feat * n_samp, 40),
                  nrow = n_feat, ncol = n_samp,
                  dimnames = list(paste0("f", seq_len(n_feat)),
                                  paste0("S", seq_len(n_samp))))
  meta <- data.frame(
    sample = paste0("S", seq_len(n_samp)),
    Env = c(rep("control", 6), rep("treatment", 6)),
    stringsAsFactors = FALSE
  )

  m_ctrl <- suppressMessages(suppressWarnings(
    pathway_daa(abundance = abund, metadata = meta, group = "Env",
                daa_method = "Maaslin2", reference = "control")
  ))
  m_trt <- suppressMessages(suppressWarnings(
    pathway_daa(abundance = abund, metadata = meta, group = "Env",
                daa_method = "Maaslin2", reference = "treatment")
  ))

  expect_true(all(m_ctrl$group1 == "control"))
  expect_true(all(m_ctrl$group2 == "treatment"))
  expect_true(all(m_trt$group1 == "treatment"))
  expect_true(all(m_trt$group2 == "control"))
  # Align rows by feature before comparing — Maaslin2 may reorder.
  ord_ctrl <- order(m_ctrl$feature)
  ord_trt <- order(m_trt$feature)
  expect_equal(m_ctrl$log2_fold_change[ord_ctrl],
               -m_trt$log2_fold_change[ord_trt],
               tolerance = 1e-6)
})

test_that("pathway_daa edgeR honors user-specified reference", {
  # Regression: edgeR's exactTest() used raw factor order (pair = c(1, 2)),
  # so the `reference` argument was silently ignored and result labels
  # were always Level[1]/Level[2]. Relevel the grouping factor so edgeR
  # and the labels agree with the documented semantics.
  skip_if_not_installed("edgeR")

  set.seed(7)
  n_feat <- 10; n_samp <- 12
  abund <- matrix(rpois(n_feat * n_samp, 40),
                  nrow = n_feat, ncol = n_samp,
                  dimnames = list(paste0("f", seq_len(n_feat)),
                                  paste0("S", seq_len(n_samp))))
  meta <- data.frame(
    sample = paste0("S", seq_len(n_samp)),
    Env = c(rep("control", 6), rep("treatment", 6)),
    stringsAsFactors = FALSE
  )

  e_ctrl <- suppressWarnings(suppressMessages(
    pathway_daa(abundance = abund, metadata = meta, group = "Env",
                daa_method = "edgeR", reference = "control")
  ))
  e_trt <- suppressWarnings(suppressMessages(
    pathway_daa(abundance = abund, metadata = meta, group = "Env",
                daa_method = "edgeR", reference = "treatment")
  ))

  expect_true(all(e_ctrl$group1 == "control"))
  expect_true(all(e_ctrl$group2 == "treatment"))
  expect_true(all(e_trt$group1 == "treatment"))
  expect_true(all(e_trt$group2 == "control"))
  expect_equal(e_ctrl$log2_fold_change, -e_trt$log2_fold_change, tolerance = 1e-6)
})

test_that("pathway_daa metagenomeSeq honors user-specified reference", {
  # Regression: metagenomeSeq labels were fixed to Level[1]/Level[2] and
  # the model matrix used raw factor order, so `reference = "treatment"`
  # left both labels and p-values unchanged.
  skip_if_not_installed("metagenomeSeq")

  set.seed(7)
  n_feat <- 10; n_samp <- 12
  abund <- matrix(rpois(n_feat * n_samp, 40),
                  nrow = n_feat, ncol = n_samp,
                  dimnames = list(paste0("f", seq_len(n_feat)),
                                  paste0("S", seq_len(n_samp))))
  meta <- data.frame(
    sample = paste0("S", seq_len(n_samp)),
    Env = c(rep("control", 6), rep("treatment", 6)),
    stringsAsFactors = FALSE
  )

  m_ctrl <- suppressWarnings(suppressMessages(
    pathway_daa(abundance = abund, metadata = meta, group = "Env",
                daa_method = "metagenomeSeq", reference = "control")
  ))
  m_trt <- suppressWarnings(suppressMessages(
    pathway_daa(abundance = abund, metadata = meta, group = "Env",
                daa_method = "metagenomeSeq", reference = "treatment")
  ))

  expect_true(all(m_ctrl$group1 == "control"))
  expect_true(all(m_ctrl$group2 == "treatment"))
  expect_true(all(m_trt$group1 == "treatment"))
  expect_true(all(m_trt$group2 == "control"))
})

test_that("pathway_daa re-validates group count after align/select", {
  # Regression: validate_group() only checks the raw metadata. If
  # align_samples() or a narrow `select =` filter removes every sample
  # of a level, the single-group leftover could propagate into backends
  # with a confusing downstream error.
  td <- create_daa_test_data(n_samples = 6, n_groups = 2)

  expect_error(
    pathway_daa(
      abundance = td$abundance,
      metadata = td$metadata,
      group = "group",
      daa_method = "ALDEx2",
      select = td$metadata$sample[td$metadata$group == "control"]
    ),
    "at least 2 groups with samples"
  )
})

test_that("pathway_daa rejects unsupported daa_method with a typo suggestion", {
  # Regression: pathway_daa() used to dispatch via switch() with no
  # default branch, so a misspelled method like "linDA" fell through and
  # returned NULL silently. Now we whitelist-validate up front and hint
  # the canonical spelling for common typos.
  td <- create_daa_test_data(n_samples = 6, n_groups = 2)

  expect_error(
    pathway_daa(
      abundance = td$abundance,
      metadata = td$metadata,
      group = "group",
      daa_method = "linDA"
    ),
    "Did you mean 'LinDA'"
  )

  expect_error(
    pathway_daa(
      abundance = td$abundance,
      metadata = td$metadata,
      group = "group",
      daa_method = "Lefse"
    ),
    "Did you mean 'Lefser'"
  )

  expect_error(
    pathway_daa(
      abundance = td$abundance,
      metadata = td$metadata,
      group = "group",
      daa_method = "totally_unknown"
    ),
    "Unsupported daa_method"
  )
})

test_that("pathway_daa rejects samples with zero total abundance", {
  # Regression: `x / sum(x)` inside calculate_abundance_stats() produced
  # NaN for zero-sum sample columns; the surrounding `mean(..., na.rm =
  # TRUE)` then silently dropped those NaN values, so the reported group
  # means were computed from fewer samples than the user had supplied.
  # The zero-sum column must now surface as an actionable error that
  # names the offending sample.
  skip_if_not_installed("ALDEx2")

  td <- create_daa_test_data(n_samples = 6, n_groups = 2)
  bad_abundance <- td$abundance
  zero_sample <- colnames(bad_abundance)[3]
  bad_abundance[, zero_sample] <- 0

  expect_error(
    pathway_daa(
      abundance = bad_abundance,
      metadata = td$metadata,
      group = "group",
      daa_method = "ALDEx2",
      include_abundance_stats = TRUE
    ),
    regexp = zero_sample
  )
})
