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

test_that("pathway_daa rejects invalid reference levels instead of silently falling back", {
  td <- create_daa_test_data(n_samples = 6)

  expect_error(
    pathway_daa(td$abundance, td$metadata, "group",
                daa_method = "ALDEx2", p_adjust_method = "bogus"),
    "'p_adjust_method' must be one of"
  )

  expect_error(
    pathway_daa(td$abundance, td$metadata, "group",
                daa_method = "ALDEx2", reference = "missing"),
    "Reference level 'missing' was not found"
  )
  expect_error(
    pathway_daa(td$abundance, td$metadata, "group",
                daa_method = "ALDEx2", reference = ""),
    "'reference' must be NULL or a single non-empty character string"
  )

  td_three <- create_daa_test_data(n_samples = 9, n_groups = 3)
  selected_samples <- td_three$metadata$sample[td_three$metadata$group != "control"]
  expect_error(
    pathway_daa(td_three$abundance, td_three$metadata, "group",
                daa_method = "ALDEx2", reference = "control",
                select = selected_samples),
    "Reference level 'control' was not found"
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

test_that("formula-based DAA methods handle non-syntactic group columns", {
  method_pkg <- c(
    "DESeq2" = "DESeq2",
    "metagenomeSeq" = "metagenomeSeq",
    "LinDA" = "MicrobiomeStat",
    "Maaslin2" = "Maaslin2"
  )
  installed_methods <- names(method_pkg)[vapply(
    method_pkg,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )]
  skip_if(length(installed_methods) == 0,
          "No formula-based DAA backend packages are installed.")

  set.seed(123)
  n_features <- 12
  n_samples <- 8
  abundance <- matrix(
    rpois(n_features * n_samples, lambda = 40) + 1,
    nrow = n_features,
    dimnames = list(paste0("p", seq_len(n_features)),
                    paste0("S", seq_len(n_samples)))
  )
  abundance[seq_len(3), 5:8] <- abundance[seq_len(3), 5:8] + 30
  abundance <- as.data.frame(abundance)

  metadata <- data.frame(
    sample = colnames(abundance),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  metadata[["treatment group"]] <- rep(
    c("baseline group", "treated group"),
    each = n_samples / 2
  )
  metadata[[".ggpicrust2_group"]] <- "preexisting column"

  for (method in installed_methods) {
    captured <- capture.output({
      result <- suppressWarnings(suppressMessages(pathway_daa(
        abundance,
        metadata,
        "treatment group",
        daa_method = method,
        reference = "baseline group"
      )))
    }, type = "output")
    ignore <- captured

    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), n_features)
    expect_true(all(result$group1 == "baseline group"))
    expect_true(all(result$group2 == "treated group"))
    expect_true(all(is.na(result$p_values) |
                      (result$p_values >= 0 & result$p_values <= 1)))
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

test_that("DESeq2 preserves method-native adjusted p-values", {
  skip_if_not_installed("DESeq2")

  set.seed(202)
  n_features <- 8
  n_samples <- 12
  abundance <- as.data.frame(matrix(
    rpois(n_features * n_samples, lambda = 40),
    nrow = n_features,
    dimnames = list(paste0("p", seq_len(n_features)),
                    paste0("S", seq_len(n_samples)))
  ))
  metadata <- data.frame(
    sample = colnames(abundance),
    group = rep(c("control", "treatment"), each = n_samples / 2),
    stringsAsFactors = FALSE
  )

  captured <- capture.output({
    res <- suppressWarnings(
      pathway_daa(
        abundance,
        metadata,
        "group",
        daa_method = "DESeq2",
        p_adjust_method = "bonferroni"
      )
    )
  }, type = "output")
  ignore <- captured

  expect_true("p_adjust" %in% colnames(res))
  expect_true(all(res$adj_method == "bonferroni (method-specific)"))
  expect_true(all(is.na(res$p_adjust) | (res$p_adjust >= 0 & res$p_adjust <= 1)))
})

test_that("DESeq2 scientific diagnostic warnings reach the caller", {
  skip_if_not_installed("DESeq2")

  abundance <- matrix(
    rep(c(10, 20, 30, 40, 50), 6),
    nrow = 5,
    dimnames = list(paste0("f", 1:5), paste0("S", 1:6))
  )
  metadata <- data.frame(
    sample = colnames(abundance),
    group = rep(c("A", "B"), each = 3),
    stringsAsFactors = FALSE
  )

  expect_warning(
    result <- suppressMessages(pathway_daa(
      abundance,
      metadata,
      "group",
      daa_method = "DESeq2"
    )),
    "all genes have equal values for all samples"
  )
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), nrow(abundance))
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
  expect_true("p_adjust" %in% colnames(res))
  expect_true(all(res$adj_method == "BH (method-specific)"))

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
  expect_true("log2_fold_change" %in% colnames(res))
  expect_true(all(is.finite(res$log2_fold_change)))
})

test_that("metagenomeSeq extracts feature-aligned log fold changes", {
  skip_if_not_installed("metagenomeSeq")

  set.seed(99)
  n_features <- 20
  n_samples <- 12
  abundance <- matrix(
    rpois(n_features * n_samples, lambda = 40),
    nrow = n_features,
    dimnames = list(
      c("f_up", "f_down", paste0("f", sprintf("%02d", 3:n_features))),
      paste0("S", seq_len(n_samples))
    )
  )
  abundance["f_up", 1:6] <- c(10, 11, 12, 13, 14, 15)
  abundance["f_up", 7:12] <- c(120, 125, 130, 135, 140, 145)
  abundance["f_down", 1:6] <- c(140, 135, 130, 125, 120, 115)
  abundance["f_down", 7:12] <- c(15, 14, 13, 12, 11, 10)
  metadata <- data.frame(
    sample = colnames(abundance),
    group = rep(c("A", "B"), each = n_samples / 2),
    stringsAsFactors = FALSE
  )

  res <- suppressWarnings(suppressMessages(
    pathway_daa(
      abundance,
      metadata,
      "group",
      daa_method = "metagenomeSeq",
      reference = "A"
    )
  ))

  expect_equal(nrow(res), n_features)
  expect_true(all(is.finite(res$log2_fold_change)))
  expect_gt(res$log2_fold_change[res$feature == "f_up"], 0)
  expect_lt(res$log2_fold_change[res$feature == "f_down"], 0)
})

test_that("metagenomeSeq survives degenerate cumNormStat input", {
  skip_if_not_installed("metagenomeSeq")

  # Minimal-but-legal input: 4 samples, 3 features, monotonic counts.
  # metagenomeSeq::cumNormStatFast() returns NaN for this shape because the
  # per-sample quantile search has nothing to stabilize on, and the package
  # then aborts with the cryptic "missing value where TRUE/FALSE needed".
  # We keep the existing p = 0.5 fallback for this degenerate search path,
  # but it must be visible to callers as a warning rather than a silent
  # normalization change.
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

  res <- NULL
  expect_warning(
    res <- pathway_daa(
      abundance,
      metadata,
      "group",
      daa_method = "metagenomeSeq"
    ),
    "could not estimate a stable metagenomeSeq CSS normalization quantile"
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), nrow(abundance))
  expect_true(all(c("feature", "method", "p_values") %in% colnames(res)))
})

test_that("metagenomeSeq fails when CSS normalization is undefined for sparse samples", {
  skip_if_not_installed("metagenomeSeq")

  abundance <- data.frame(
    S1 = c(10, 0, 0),
    S2 = c(12, 0, 0),
    S3 = c(0, 10, 0),
    S4 = c(0, 11, 0),
    row.names = paste0("p", 1:3)
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    pathway_daa(abundance, metadata, "group", daa_method = "metagenomeSeq"),
    "requires at least two positive features in every sample"
  )
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
    "non-negative values"
  )
})

test_that("validate_daa_input rejects missing and non-finite abundance values", {
  validate_daa <- getFromNamespace("validate_daa_input", "ggpicrust2")
  mat <- matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8),
    nrow = 2,
    dimnames = list(c("p1", "p2"), paste0("S", 1:4))
  )

  mat_missing <- mat
  mat_missing[1, 2] <- NA_real_
  expect_error(
    validate_daa(mat_missing),
    "must not contain missing values"
  )

  mat_infinite <- mat
  mat_infinite[2, 3] <- Inf
  expect_error(
    validate_daa(mat_infinite),
    "must contain only finite values"
  )
})

test_that("pathway_daa warns when count-based methods round non-integer abundance", {
  skip_if_not_installed("edgeR")

  abundance <- data.frame(
    S1 = c(10.2, 20.7, 30.1),
    S2 = c(11.4, 21.2, 31.8),
    S3 = c(40.5, 15.3, 25.1),
    S4 = c(41.6, 16.9, 26.4),
    row.names = paste0("p", 1:3)
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    suppressMessages(
      pathway_daa(abundance, metadata, "group", daa_method = "edgeR")
    ),
    "non-integer abundance value\\(s\\) were rounded"
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

test_that("wrapper-computed DAA p-values are adjusted within each contrast", {
  result <- data.frame(
    feature = rep(c("p1", "p2"), 2),
    method = "edgeR",
    group1 = "A",
    group2 = rep(c("B", "C"), each = 2),
    p_values = c(0.01, 0.02, 0.01, 0.90),
    stringsAsFactors = FALSE
  )

  adjusted <- ggpicrust2:::adjust_daa_p_values(result, "bonferroni")

  expect_equal(
    adjusted$p_adjust[adjusted$group2 == "B"],
    stats::p.adjust(c(0.01, 0.02), method = "bonferroni")
  )
  expect_equal(
    adjusted$p_adjust[adjusted$group2 == "C"],
    stats::p.adjust(c(0.01, 0.90), method = "bonferroni")
  )
  expect_equal(adjusted$adj_method, rep("bonferroni", nrow(adjusted)))
})

test_that("DAA p-value adjustment validates raw p-values before adjustment", {
  result <- data.frame(
    feature = "p1",
    method = "edgeR",
    group1 = "A",
    group2 = "B",
    p_values = 1.2,
    stringsAsFactors = FALSE
  )

  expect_error(
    ggpicrust2:::adjust_daa_p_values(result, "BH"),
    "Column 'p_values' in DAA result must contain finite values between 0 and 1 or NA"
  )
})

test_that("DAA p-value adjustment validates method-specific adjusted p-values", {
  result <- data.frame(
    feature = "p1",
    method = "DESeq2",
    group1 = "A",
    group2 = "B",
    p_values = 0.2,
    p_adjust = 1.5,
    stringsAsFactors = FALSE
  )
  attr(result, "adj_method") <- "BH (method-specific)"

  expect_error(
    ggpicrust2:::adjust_daa_p_values(result, "BH"),
    "Column 'p_adjust' in DAA result must contain finite values between 0 and 1 or NA"
  )
})

test_that("DAA p-value adjustment preserves valid missing p-values", {
  result <- data.frame(
    feature = c("p1", "p2", "p3"),
    method = "limma voom",
    group1 = "A",
    group2 = "B",
    p_values = c(NA_real_, 0.01, 0.02),
    stringsAsFactors = FALSE
  )

  adjusted <- ggpicrust2:::adjust_daa_p_values(result, "BH")

  expect_equal(adjusted$p_adjust, c(NA_real_, 0.02, 0.02))
  expect_equal(adjusted$adj_method, rep("BH", nrow(adjusted)))
})

test_that("backend result values are aligned and validated by feature IDs", {
  align_values <- getFromNamespace(
    "validate_and_align_backend_result_values",
    "ggpicrust2"
  )

  aligned <- align_values(
    values = c(0.3, 0.2, 0.1),
    output_ids = c("f3", "f2", "f1"),
    feature_ids = c("f1", "f2", "f3"),
    value_name = "pvalue",
    context = "synthetic backend",
    value_type = "probability",
    allow_na = FALSE
  )
  expect_equal(aligned, c(0.1, 0.2, 0.3))

  expect_error(
    align_values(
      values = c(0.3, 0.2, 0.1),
      output_ids = c("f3", "fX", "f1"),
      feature_ids = c("f1", "f2", "f3"),
      value_name = "pvalue",
      context = "synthetic backend",
      value_type = "probability",
      allow_na = FALSE
    ),
    "backend feature identifiers do not match"
  )
  expect_error(
    align_values(
      values = c(0.3, 1.2, 0.1),
      output_ids = c("f3", "f2", "f1"),
      feature_ids = c("f1", "f2", "f3"),
      value_name = "pvalue",
      context = "synthetic backend",
      value_type = "probability",
      allow_na = FALSE
    ),
    "between 0 and 1"
  )
  expect_error(
    align_values(
      values = c(3, Inf, 1),
      output_ids = c("f3", "f2", "f1"),
      feature_ids = c("f1", "f2", "f3"),
      value_name = "logFC",
      context = "synthetic backend",
      value_type = "finite_numeric",
      allow_na = FALSE
    ),
    "finite numeric"
  )
})

test_that("limma voom backend output is aligned by feature row names", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")

  td <- create_daa_test_data(n_samples = 6)
  feature_ids <- rownames(td$abundance)
  backend_ids <- rev(feature_ids)
  backend_p <- setNames(seq_along(backend_ids) / 10, backend_ids)
  backend_lfc <- setNames(seq_along(backend_ids) * 10, backend_ids)

  testthat::local_mocked_bindings(
    eBayes = function(fit, ...) {
      coef_names <- colnames(fit$coefficients)
      p_mat <- matrix(
        0.9,
        nrow = length(backend_ids),
        ncol = length(coef_names),
        dimnames = list(backend_ids, coef_names)
      )
      coef_mat <- matrix(
        0,
        nrow = length(backend_ids),
        ncol = length(coef_names),
        dimnames = list(backend_ids, coef_names)
      )
      p_mat[, 2] <- unname(backend_p[backend_ids])
      coef_mat[, 2] <- unname(backend_lfc[backend_ids])
      fit$p.value <- p_mat
      fit$coefficients <- coef_mat
      fit
    },
    .package = "limma"
  )

  res <- suppressWarnings(suppressMessages(
    pathway_daa(
      td$abundance,
      td$metadata,
      "group",
      daa_method = "limma voom"
    )
  ))

  expect_equal(res$feature, feature_ids)
  expect_equal(res$p_values, unname(backend_p[feature_ids]))
  expect_equal(res$log2_fold_change, unname(backend_lfc[feature_ids]))
})

test_that("edgeR backend output is aligned by exactTest row names", {
  skip_if_not_installed("edgeR")

  td <- create_daa_test_data(n_samples = 6)
  feature_ids <- rownames(td$abundance)
  backend_ids <- rev(feature_ids)
  backend_p <- setNames(seq_along(backend_ids) / 10, backend_ids)
  backend_lfc <- setNames(seq_along(backend_ids) * 10, backend_ids)
  fake_table <- data.frame(
    logFC = unname(backend_lfc[backend_ids]),
    logCPM = seq_along(backend_ids),
    PValue = unname(backend_p[backend_ids]),
    row.names = backend_ids
  )

  testthat::local_mocked_bindings(
    exactTest = function(...) {
      list(table = fake_table)
    },
    .package = "edgeR"
  )

  res <- suppressWarnings(suppressMessages(
    pathway_daa(
      td$abundance,
      td$metadata,
      "group",
      daa_method = "edgeR"
    )
  ))

  expect_equal(res$feature, feature_ids)
  expect_equal(res$p_values, unname(backend_p[feature_ids]))
  expect_equal(res$log2_fold_change, unname(backend_lfc[feature_ids]))
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

test_that("ALDEx2 effect-size output is aligned and validated by feature", {
  validate_effects <- getFromNamespace(
    "validate_and_align_aldex2_effect_results",
    "ggpicrust2"
  )
  effect_results <- data.frame(
    effect = c(2, 1),
    diff.btw = c(20, 10),
    rab.all = c(200, 100),
    overlap = c(0.2, 0.1),
    row.names = c("feature2", "feature1")
  )

  aligned <- validate_effects(
    effect_results,
    c("feature1", "feature2")
  )
  expect_equal(rownames(aligned), c("feature1", "feature2"))
  expect_equal(aligned$diff.btw, c(10, 20))

  expect_error(
    validate_effects(
      effect_results[, setdiff(colnames(effect_results), "effect")],
      c("feature1", "feature2")
    ),
    "missing required column.*effect"
  )
  invalid_overlap <- effect_results
  invalid_overlap$overlap[1] <- 1.5
  expect_error(
    validate_effects(invalid_overlap, c("feature1", "feature2")),
    "overlap.*between 0 and 1"
  )
})

test_that("ALDEx2 requested effect-size failures stop the analysis", {
  skip_if_not_installed("ALDEx2")

  td <- create_daa_test_data(n_samples = 4)
  local_mocked_bindings(
    aldex.effect = function(...) {
      stop("synthetic effect failure")
    },
    .package = "ALDEx2"
  )

  expect_error(
    pathway_daa(
      td$abundance,
      td$metadata,
      "group",
      daa_method = "ALDEx2",
      include_effect_size = TRUE
    ),
    paste0(
      "effect-size calculation failed while .*include_effect_size = TRUE.*",
      "synthetic effect failure"
    )
  )
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

test_that("include_abundance_stats fails when backend features cannot be summarized", {
  skip_if_not_installed("limma")

  td <- create_daa_test_data(n_samples = 4)
  testthat::local_mocked_bindings(
    perform_limma_voom_analysis = function(...) {
      data.frame(
        feature = "feature_not_in_abundance",
        method = "limma voom",
        group1 = "control",
        group2 = "treatment",
        p_values = 0.5,
        stringsAsFactors = FALSE
      )
    }
  )

  expect_error(
    pathway_daa(
      td$abundance,
      td$metadata,
      "group",
      daa_method = "limma voom",
      include_abundance_stats = TRUE
    ),
    "Failed to calculate abundance statistics.*missing from abundance row names"
  )
})

test_that("pathway_daa ALDEx2 honors user-specified reference", {
  # Regression: ALDEx2 received `as.numeric(Group)` from the raw factor order,
  # so `reference` changed neither group labels nor the direction of the
  # CLR-space `diff.btw`/`log2_fold_change` effect size.
  skip_if_not_installed("ALDEx2")

  counts <- matrix(c(
    100, 100, 100, 100, 300, 300, 300, 300,
    300, 300, 300, 300, 100, 100, 100, 100,
    50, 50, 50, 50, 50, 50, 50, 50,
    80, 80, 80, 80, 80, 80, 80, 80
  ), nrow = 4, byrow = TRUE)
  rownames(counts) <- paste0("f", seq_len(nrow(counts)))
  colnames(counts) <- paste0("S", seq_len(ncol(counts)))
  meta <- data.frame(
    sample = colnames(counts),
    Env = c(rep("control", 4), rep("treatment", 4)),
    stringsAsFactors = FALSE
  )

  set.seed(24)
  a_ctrl <- suppressMessages(suppressWarnings(
    pathway_daa(abundance = counts, metadata = meta, group = "Env",
                daa_method = "ALDEx2", reference = "control")
  ))
  set.seed(24)
  a_trt <- suppressMessages(suppressWarnings(
    pathway_daa(abundance = counts, metadata = meta, group = "Env",
                daa_method = "ALDEx2", reference = "treatment")
  ))

  expect_true(all(a_ctrl$group1 == "control"))
  expect_true(all(a_ctrl$group2 == "treatment"))
  expect_true(all(a_trt$group1 == "treatment"))
  expect_true(all(a_trt$group2 == "control"))

  ctrl_welch <- a_ctrl[a_ctrl$method == "ALDEx2_Welch's t test", ]
  trt_welch <- a_trt[a_trt$method == "ALDEx2_Welch's t test", ]
  ord_ctrl <- order(ctrl_welch$feature)
  ord_trt <- order(trt_welch$feature)
  expect_equal(
    ctrl_welch$log2_fold_change[ord_ctrl],
    -trt_welch$log2_fold_change[ord_trt],
    tolerance = 0.05
  )
  expect_true(ctrl_welch$log2_fold_change[ctrl_welch$feature == "f1"] > 0)
  expect_true(trt_welch$log2_fold_change[trt_welch$feature == "f1"] < 0)
  expect_true(ctrl_welch$log2_fold_change[ctrl_welch$feature == "f2"] < 0)
  expect_true(trt_welch$log2_fold_change[trt_welch$feature == "f2"] > 0)
  expect_equal(ctrl_welch$log2_fold_change, ctrl_welch$diff_btw)
})

test_that("format_linda_output rejects malformed LinDA outputs", {
  malformed_shape <- list(
    groupB = data.frame(
      pvalue = I(matrix(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06), nrow = 3)),
      padj = c(0.02, 0.03, 0.04),
      log2FoldChange = c(1.2, 0.8, 0.5),
      row.names = c("path1", "path2", "path3")
    )
  )
  missing_column <- list(
    groupC = data.frame(
      some_other_col = c(1, 2),
      row.names = c("path4", "path5")
    )
  )
  invalid_probability <- list(
    groupB = data.frame(
      pvalue = c(0.01, 1.2),
      padj = c(0.02, 0.5),
      log2FoldChange = c(1.2, 0.8),
      row.names = c("path1", "path2")
    )
  )

  expect_error(
    ggpicrust2:::format_linda_output(
      linda_output = malformed_shape,
      group = "group",
      reference = "A",
      Level = c("A", "B", "C")
    ),
    "pvalue.*invalid length or shape"
  )
  expect_error(
    ggpicrust2:::format_linda_output(
      linda_output = missing_column,
      group = "group",
      reference = "A",
      Level = c("A", "B", "C")
    ),
    "missing required column 'pvalue'"
  )
  expect_error(
    ggpicrust2:::format_linda_output(
      linda_output = invalid_probability,
      group = "group",
      reference = "A",
      Level = c("A", "B")
    ),
    "between 0 and 1"
  )
})

test_that("format_linda_output preserves LinDA adjusted p-values", {
  linda_output <- list(
    groupB = data.frame(
      pvalue = c(0.01, 0.2),
      padj = c(0.02, 0.2),
      log2FoldChange = c(1.5, -0.3),
      row.names = c("path1", "path2")
    )
  )

  result <- ggpicrust2:::format_linda_output(
    linda_output = linda_output,
    group = "group",
    reference = "A",
    Level = c("A", "B")
  )

  expect_true("p_adjust" %in% colnames(result))
  expect_equal(result$p_adjust, c(0.02, 0.2))
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

test_that("pathway_daa Lefser honors reference and tests relative abundance", {
  skip_if_not_installed("lefser")

  abundance <- data.frame(
    S1 = c(900, 100),
    S2 = c(800, 200),
    S3 = c(90000, 10000),
    S4 = c(80000, 20000),
    row.names = c("f1", "f2")
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    Env = c("control", "control", "treatment", "treatment"),
    stringsAsFactors = FALSE
  )

  res <- suppressMessages(suppressWarnings(
    pathway_daa(
      abundance = abundance,
      metadata = metadata,
      group = "Env",
      daa_method = "Lefser",
      reference = "treatment"
    )
  ))

  rel_abundance <- sweep(as.matrix(abundance), 2, colSums(abundance), "/") * 1e6
  group_vector <- factor(metadata$Env, levels = c("treatment", "control"))
  expected_p <- apply(rel_abundance, 1, function(feature_values) {
    stats::kruskal.test(feature_values ~ group_vector)$p.value
  })

  expect_true(all(res$group1 == "treatment"))
  expect_true(all(res$group2 == "control"))
  expect_equal(
    res$p_values[match(rownames(rel_abundance), res$feature)],
    unname(expected_p),
    tolerance = 1e-12
  )
  expect_equal(unname(expected_p), c(1, 1), tolerance = 1e-12)
})

test_that("pathway_daa preserves a leading feature ID column before sample alignment", {
  skip_if_not_installed("limma")

  abundance <- data.frame(
    feature = c("F_A", "F_B", "F_C", "F_D", "F_E"),
    S1 = c(10, 1, 1, 2, 3),
    S2 = c(11, 1, 1, 2, 3),
    S3 = c(1, 10, 1, 2, 3),
    S4 = c(1, 11, 1, 2, 3),
    check.names = FALSE
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    Env = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  res <- suppressWarnings(pathway_daa(
    abundance = abundance,
    metadata = metadata,
    group = "Env",
    daa_method = "limma voom"
  ))

  expect_setequal(res$feature, abundance$feature)
  expect_false(any(res$feature %in% as.character(seq_len(nrow(abundance)))))
})

test_that("pathway_daa requires explicit unique feature identifiers", {
  skip_if_not_installed("limma")

  metadata <- data.frame(
    sample = paste0("S", 1:4),
    Env = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  no_feature_ids <- data.frame(
    S1 = c(10, 20, 30),
    S2 = c(11, 21, 31),
    S3 = c(30, 10, 20),
    S4 = c(31, 11, 21),
    check.names = FALSE
  )
  expect_error(
    pathway_daa(
      abundance = no_feature_ids,
      metadata = metadata,
      group = "Env",
      daa_method = "limma voom"
    ),
    "feature identifiers|row names"
  )

  duplicated_feature_ids <- matrix(
    c(
      10, 11, 30, 31,
      20, 21, 10, 11,
      30, 31, 20, 21
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("F1", "F1", "F3"), paste0("S", 1:4))
  )
  expect_error(
    pathway_daa(
      abundance = duplicated_feature_ids,
      metadata = metadata,
      group = "Env",
      daa_method = "limma voom"
    ),
    "duplicated feature identifiers"
  )
})

test_that("perform_lefser_analysis surfaces backend errors instead of dropping LDA scores", {
  skip_if_not_installed("lefser")

  abundance <- matrix(
    c(1, 2, 3, 4,
      4, 3, 2, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("f1", "f2"), paste0("S", 1:4))
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    Env = rep("A", 4),
    stringsAsFactors = FALSE
  )

  expect_error(
    ggpicrust2:::perform_lefser_analysis(
      abundance_mat = abundance,
      metadata = metadata,
      group = "Env",
      reference = "A",
      Level = c("A", "B")
    ),
    "Lefser analysis failed"
  )
})

test_that("calculate_lefser_kruskal_p_value fails instead of returning p=1 for invalid tests", {
  calc <- ggpicrust2:::calculate_lefser_kruskal_p_value

  expect_equal(
    calc(
      feature_values = c(5, 5, 5, 5),
      group_vector = factor(c("A", "A", "B", "B")),
      feature_id = "constant_feature"
    ),
    1
  )

  expect_error(
    calc(
      feature_values = c(1, 2, 3),
      group_vector = factor(c("A", "A", "B", "B")),
      feature_id = "bad_length"
    ),
    "bad_length.*length"
  )
  expect_error(
    calc(
      feature_values = c(1, 2, Inf, 4),
      group_vector = factor(c("A", "A", "B", "B")),
      feature_id = "bad_value"
    ),
    "bad_value.*finite"
  )
  expect_error(
    calc(
      feature_values = c(1, 2, 3, 4),
      group_vector = factor(c("A", "A", NA, NA)),
      feature_id = "bad_group"
    ),
    "bad_group"
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

test_that("perform_linda_analysis surfaces backend errors instead of returning empty results", {
  skip_if_not_installed("MicrobiomeStat")

  abundance <- matrix(
    rpois(12 * 4, lambda = 20),
    nrow = 12,
    dimnames = list(paste0("f", 1:12), paste0("S", 1:4))
  )
  metadata <- data.frame(
    sample = paste0("S", 1:4),
    Env = rep("A", 4),
    stringsAsFactors = FALSE
  )

  expect_error(
    ggpicrust2:::perform_linda_analysis(
      abundance = abundance,
      metadata = metadata,
      group = "Env",
      reference = "A",
      Level = c("A", "B"),
      length_Level = 2,
      p_adjust_method = "BH"
    ),
    "LinDA analysis failed"
  )
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

test_that("Maaslin2 feature IDs are mapped through make.names without ambiguity", {
  prepare_map <- getFromNamespace("prepare_maaslin2_feature_mapping", "ggpicrust2")
  validate_results <- getFromNamespace("validate_maaslin2_results", "ggpicrust2")

  expect_error(
    prepare_map(c("path-a", "path.a", "path b")),
    "make.names\\(\\).*ambiguous"
  )

  feature_map <- prepare_map(c("path-a", "path b"))
  expect_equal(names(feature_map), c("path.a", "path.b"))

  maaslin2_results <- data.frame(
    feature = c("path.a", "path.b"),
    metadata = c("Env", "Env"),
    value = c("B", "B"),
    coef = c(1.2, -0.5),
    pval = c(0.01, 0.20),
    qval = c(0.02, 0.20),
    stringsAsFactors = FALSE
  )
  validated <- validate_results(maaslin2_results, "Env", feature_map)
  expect_equal(validated$feature, c("path-a", "path b"))

  bad_feature <- maaslin2_results
  bad_feature$feature[1] <- "unexpected"
  expect_error(
    validate_results(bad_feature, "Env", feature_map),
    "do not match the sanitized"
  )

  bad_p <- maaslin2_results
  bad_p$pval[1] <- 1.2
  expect_error(
    validate_results(bad_p, "Env", feature_map),
    "between 0 and 1"
  )
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

test_that("pathway_daa metagenomeSeq emits one block per non-reference level for >=3 groups", {
  # Regression: metagenomeSeq previously built a full k-column model matrix,
  # called fit_feature_model() once, read `coef = 2`, and hard-coded the
  # output as `group1 = Level[1] / group2 = Level[2]` regardless of how
  # many groups existed. Contrasts against any level beyond the first
  # non-reference one were silently dropped. The result should now have
  # `(k - 1) * n_features` rows, one (ref, non-ref) block per contrast,
  # with both the coefficient and the labels redone per pair (since
  # fitFeatureModel is documented as a two-group entry point).
  skip_if_not_installed("metagenomeSeq")

  set.seed(11)
  n_feat <- 8
  per_group <- 6
  groups <- rep(c("A", "B", "C"), each = per_group)
  n_samp <- length(groups)
  abund <- matrix(rpois(n_feat * n_samp, 50),
                  nrow = n_feat, ncol = n_samp,
                  dimnames = list(paste0("f", seq_len(n_feat)),
                                  paste0("S", seq_len(n_samp))))
  meta <- data.frame(
    sample = paste0("S", seq_len(n_samp)),
    Env = groups,
    stringsAsFactors = FALSE
  )

  res <- suppressWarnings(suppressMessages(
    pathway_daa(abundance = abund, metadata = meta, group = "Env",
                daa_method = "metagenomeSeq")
  ))

  # Shape: k = 3 groups -> 2 non-reference contrasts -> 2 * n_feat rows.
  expect_equal(nrow(res), 2 * n_feat)

  # Reference defaults to the first level ("A"). group1 must always be A,
  # group2 must cover exactly {B, C} once each per feature.
  expect_true(all(res$group1 == "A"))
  expect_setequal(unique(res$group2), c("B", "C"))
  expect_equal(sum(res$group2 == "B"), n_feat)
  expect_equal(sum(res$group2 == "C"), n_feat)

  # reference = "B" must flip group1 to B and restrict group2 to {A, C}.
  res_refB <- suppressWarnings(suppressMessages(
    pathway_daa(abundance = abund, metadata = meta, group = "Env",
                daa_method = "metagenomeSeq", reference = "B")
  ))
  expect_equal(nrow(res_refB), 2 * n_feat)
  expect_true(all(res_refB$group1 == "B"))
  expect_setequal(unique(res_refB$group2), c("A", "C"))

  # p_values and feature-aligned log2_fold_change columns must both be
  # present across the full output.
  expect_true("p_values" %in% colnames(res))
  expect_true("log2_fold_change" %in% colnames(res))
  expect_true(all(is.finite(res$log2_fold_change)))
})

test_that("pathway_daa re-validates group count after align/select", {
  # Regression: validate_group() only checks the raw metadata. If
  # align_samples() or a narrow `select =` filter removes every sample
  # of a level, the single-group leftover could propagate into backends
  # with a confusing downstream error.
  td <- create_daa_test_data(n_samples = 8, n_groups = 2)

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

test_that("pathway_daa re-validates sample count and uniqueness after select", {
  # Regression: select= was applied after the aligned sample-count check.
  # A caller could start with a valid dataset, select fewer than four
  # samples, and still reach backend fitting.
  td <- create_daa_test_data(n_samples = 6, n_groups = 2)

  expect_error(
    pathway_daa(
      abundance = td$abundance,
      metadata = td$metadata,
      group = "group",
      daa_method = "ALDEx2",
      select = c("sample1", "sample2", "sample4")
    ),
    "At least 4 samples required for DAA after `select` filtering"
  )

  expect_error(
    pathway_daa(
      abundance = td$abundance,
      metadata = td$metadata,
      group = "group",
      daa_method = "ALDEx2",
      select = c("sample1", "sample1", "sample4", "sample5")
    ),
    "'select' must contain unique sample names"
  )
})

test_that("pathway_daa rejects missing groups and singleton groups after alignment/select", {
  td <- create_daa_test_data(n_samples = 6, n_groups = 2)

  metadata_missing <- td$metadata
  metadata_missing$group[2] <- NA
  expect_error(
    pathway_daa(
      abundance = td$abundance,
      metadata = metadata_missing,
      group = "group",
      daa_method = "ALDEx2"
    ),
    "contains missing or empty values"
  )

  td_unbalanced <- create_daa_test_data(n_samples = 6, n_groups = 2)
  td_unbalanced$metadata$group <- c("control", "control", "control",
                                    "treatment", "treatment", "other")
  expect_error(
    pathway_daa(
      abundance = td_unbalanced$abundance,
      metadata = td_unbalanced$metadata,
      group = "group",
      daa_method = "ALDEx2"
    ),
    "at least 2 samples per group"
  )

  td_three <- create_daa_test_data(n_samples = 9, n_groups = 3)
  selected <- c(
    td_three$metadata$sample[td_three$metadata$group == "control"],
    td_three$metadata$sample[td_three$metadata$group == "treatment"][1]
  )
  expect_error(
    pathway_daa(
      abundance = td_three$abundance,
      metadata = td_three$metadata,
      group = "group",
      daa_method = "ALDEx2",
      select = selected
    ),
    "at least 2 samples per group"
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

test_that("pathway_daa rejects backend arguments passed through dots", {
  # Regression: pathway_daa() documented `...` as backend passthrough, but
  # never forwarded those arguments. That is especially dangerous for
  # formula/fixed-effect/covariate arguments because users may believe the
  # DAA model was adjusted when it was not.
  td <- create_daa_test_data(n_samples = 6, n_groups = 2)

  expect_error(
    pathway_daa(
      abundance = td$abundance,
      metadata = td$metadata,
      group = "group",
      daa_method = "ALDEx2",
      formula = ~ group + batch
    ),
    "does not currently pass arguments in `...`"
  )

  expect_error(
    pathway_daa(
      abundance = td$abundance,
      metadata = td$metadata,
      group = "group",
      daa_method = "ALDEx2",
      select = NULL,
      p_adjust_method = "BH",
      reference = NULL,
      include_abundance_stats = FALSE,
      include_effect_size = TRUE,
      p.adjust = NULL,
      .pre_aligned = FALSE,
      .sample_col = NULL,
      "unused"
    ),
    "Unsupported argument\\(s\\): ..1"
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
