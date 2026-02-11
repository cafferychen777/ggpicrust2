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

  # Without abundance stats (default)
  result_basic <- pathway_daa(td$abundance, td$metadata, "group",
                             daa_method = "ALDEx2", include_abundance_stats = FALSE)

  abundance_cols <- c("mean_rel_abundance_group1", "sd_rel_abundance_group1",
                     "mean_rel_abundance_group2", "sd_rel_abundance_group2",
                     "log2_fold_change")

  expect_false(any(abundance_cols %in% colnames(result_basic)))

  # With abundance stats
  result_enhanced <- pathway_daa(td$abundance, td$metadata, "group",
                                daa_method = "ALDEx2", include_abundance_stats = TRUE)

  expect_true(all(abundance_cols %in% colnames(result_enhanced)))

  for (col in abundance_cols) {
    expect_true(is.numeric(result_enhanced[[col]]))
  }
  expect_true(all(is.finite(result_enhanced$log2_fold_change)))
  expect_true(all(result_enhanced$sd_rel_abundance_group1 >= 0, na.rm = TRUE))
  expect_true(all(result_enhanced$sd_rel_abundance_group2 >= 0, na.rm = TRUE))
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
