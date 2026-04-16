# Tests for pathway_annotation function

# Regression: file mode silently ignored ko_to_kegg and organism. A caller
# who passed file + ko_to_kegg = TRUE expected KEGG API annotations but got
# local reference annotations with no indication. Now we warn.
test_that("pathway_annotation warns when file mode receives KEGG-specific params", {
  temp_file <- tempfile(fileext = ".tsv")
  test_data <- data.frame(
    feature = c("K00001", "K00002"),
    sample1 = c(1, 2),
    sample2 = c(3, 4)
  )
  suppressMessages(write.table(test_data, temp_file, sep = "\t", row.names = FALSE))
  on.exit(unlink(temp_file))

  expect_warning(
    suppressMessages(pathway_annotation(file = temp_file, pathway = "KO", ko_to_kegg = TRUE)),
    "ko_to_kegg is ignored in file mode"
  )
  expect_warning(
    suppressMessages(pathway_annotation(file = temp_file, pathway = "KO", organism = "hsa")),
    "organism is ignored in file mode"
  )
})

test_that("pathway_annotation basic functionality works", {
  temp_file <- tempfile(fileext = ".tsv")
  test_data <- data.frame(
    feature = c("K00001", "K00002"),
    sample1 = c(1, 2),
    sample2 = c(3, 4)
  )
  suppressMessages(write.table(test_data, temp_file, sep = "\t", row.names = FALSE))

  result <- suppressMessages(pathway_annotation(file = temp_file, pathway = "KO", ko_to_kegg = FALSE))

  expect_s3_class(result, "data.frame")
  expect_true("description" %in% colnames(result))
  expect_equal(nrow(result), 2)

  # Regression: the file-mode branch previously took sample column names as
  # features, so `description` came back all-NA even for well-known KOs.
  # After the fix, rows whose feature IDs live in the reference data should
  # carry non-NA descriptions.
  expect_true(any(!is.na(result$description)))

  unlink(temp_file)
})

test_that("pathway_annotation works with daa_results_df input", {
  test_daa_df <- data.frame(
    feature = c("K00001", "K00002"),
    p_adjust = c(0.04, 0.03),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(pathway_annotation(pathway = "KO", daa_results_df = test_daa_df, ko_to_kegg = FALSE))

  expect_s3_class(result, "data.frame")
  expect_true("description" %in% colnames(result))
  expect_equal(nrow(result), 2)
})

test_that("pathway_annotation works with ko_to_kegg", {
  skip_if(
    Sys.getenv("GGPICRUST2_RUN_NETWORK_TESTS", "false") != "true",
    "Set GGPICRUST2_RUN_NETWORK_TESTS=true to run network-dependent KEGG tests."
  )
  skip_if_offline()

  test_daa_df <- data.frame(
    feature = c("K00001", "K00002"),
    p_adjust = c(0.04, 0.03),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(pathway_annotation(pathway = "KO", daa_results_df = test_daa_df, ko_to_kegg = TRUE))

  expect_s3_class(result, "data.frame")
  expect_true(all(c("pathway_name", "pathway_class") %in% colnames(result)))
})

test_that("pathway_annotation fills description/class/map for pathway IDs with ko_to_kegg", {
  skip_if(
    Sys.getenv("GGPICRUST2_RUN_NETWORK_TESTS", "false") != "true",
    "Set GGPICRUST2_RUN_NETWORK_TESTS=true to run network-dependent KEGG tests."
  )
  skip_if_offline()

  # KEGG pathway IDs from typical ko2kegg_abundance output
  test_daa_df <- data.frame(
    feature = c("ko05340", "ko00564"),
    p_adjust = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(pathway_annotation(pathway = "KO", daa_results_df = test_daa_df, ko_to_kegg = TRUE))

  expect_s3_class(result, "data.frame")
  expect_true(all(c("pathway_name", "pathway_description", "pathway_class", "pathway_map") %in% colnames(result)))
  # pathway_description is optional for some pathway entries in KEGG;
  # regression check is that it is no longer uniformly NA.
  expect_true(any(!is.na(result$pathway_description)))
  expect_true(any(!is.na(result$pathway_class)))
  expect_true(any(!is.na(result$pathway_map)))
})

test_that("pathway_annotation works with all pathway types", {
  # Each annotator gets its own correctly-typed features
  ko_df <- data.frame(feature = c("K00001", "K00002"), p_adjust = c(.04, .03), stringsAsFactors = FALSE)
  ec_df <- data.frame(feature = c("EC:1.1.1.1", "EC:2.7.1.1"), p_adjust = c(.04, .03), stringsAsFactors = FALSE)
  metacyc_df <- data.frame(feature = c("PWY-7219", "GLYCOLYSIS"), p_adjust = c(.04, .03), stringsAsFactors = FALSE)

  ko_result <- suppressMessages(pathway_annotation(pathway = "KO", daa_results_df = ko_df, ko_to_kegg = FALSE))
  expect_s3_class(ko_result, "data.frame")

  ec_result <- suppressMessages(pathway_annotation(pathway = "EC", daa_results_df = ec_df, ko_to_kegg = FALSE))
  expect_s3_class(ec_result, "data.frame")

  metacyc_result <- suppressMessages(pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_df, ko_to_kegg = FALSE))
  expect_s3_class(metacyc_result, "data.frame")
})

test_that("pathway_annotation ko_to_kegg returns every input row (merges back)", {
  # Regression: process_kegg_annotations() used to return only the
  # p_adjust < threshold subset, silently dropping non-significant
  # features that downstream code (e.g. ggpicrust2()'s
  # plot_result_list$daa_results_df) expected to still be present.
  # Annotations should be populated for significant rows and NA for
  # the rest, but every row must survive the call.
  cache <- getFromNamespace("kegg_cache", "ggpicrust2")

  # Pre-populate the cache with fake KEGG entries so the function never
  # touches the network.
  fake_entry <- function(name, path_id, path_desc) {
    pathway_vec <- stats::setNames(path_desc, path_id)
    map_vec <- stats::setNames("map_title", path_id)
    list(list(
      NAME = name,
      PATHWAY = pathway_vec,
      PATHWAY_MAP = map_vec,
      CLASS = "Metabolism; Carbohydrate metabolism"
    ))
  }
  assign("K00001", fake_entry("alcohol dehydrogenase", "ko00010", "Glycolysis"), envir = cache)
  assign("K00002", fake_entry("alcohol dehydrogenase NADP+", "ko00010", "Glycolysis"), envir = cache)
  on.exit({
    if (exists("K00001", envir = cache)) rm("K00001", envir = cache)
    if (exists("K00002", envir = cache)) rm("K00002", envir = cache)
    if (exists("K99999", envir = cache)) rm("K99999", envir = cache)
  }, add = TRUE)
  assign("K99999", fake_entry("non-significant placeholder", "ko99999", "Other"),
         envir = cache)

  daa_df <- data.frame(
    feature  = c("K00001", "K00002", "K99999"),
    p_values = c(0.001, 0.002, 0.4),
    p_adjust = c(0.01, 0.02, 0.9),
    method   = "ALDEx2",
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(suppressWarnings(
    pathway_annotation(pathway = "KO", daa_results_df = daa_df,
                       ko_to_kegg = TRUE, p_adjust_threshold = 0.05)
  ))

  # All input rows preserved.
  expect_equal(nrow(result), 3)
  expect_setequal(result$feature, c("K00001", "K00002", "K99999"))

  # Significant rows carry annotations.
  sig_rows <- result[result$feature %in% c("K00001", "K00002"), , drop = FALSE]
  expect_true(all(!is.na(sig_rows$pathway_name)))

  # Non-significant row has NA annotation columns (filtered out of the
  # annotation step, but still present as a row).
  ns_row <- result[result$feature == "K99999", , drop = FALSE]
  expect_true(is.na(ns_row$pathway_name))
})

test_that("pathway_annotation validates inputs", {
  expect_error(pathway_annotation(), "Please input")

  # Test invalid pathway type through daa_results_df path (avoids file-existence check)
  test_df <- data.frame(feature = "K00001", p_adjust = 0.04, stringsAsFactors = FALSE)
  expect_error(
    suppressMessages(pathway_annotation(pathway = "INVALID", daa_results_df = test_df, ko_to_kegg = FALSE)),
    "Unknown reference type"
  )

  empty_df <- data.frame(feature = character(0), p_adjust = numeric(0))
  expect_error(suppressMessages(pathway_annotation(pathway = "KO", daa_results_df = empty_df)), "empty")
})

# Regression: ko_to_kegg = TRUE used to silently ignore `pathway`, so a caller
# that passed pathway = "EC" or "MetaCyc" together with ko_to_kegg = TRUE was
# quietly routed into the KEGG (KO-only) code path. That surfaced either as a
# generic HTTP 404 or as mis-shaped enzyme records. The parameter combination
# is semantically contradictory and should fail at the function boundary.
test_that("pathway_annotation rejects ko_to_kegg = TRUE with non-KO pathway", {
  test_df <- data.frame(
    feature = c("EC:1.1.1.1", "EC:2.7.1.1"),
    p_adjust = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  expect_error(
    pathway_annotation(pathway = "EC", daa_results_df = test_df, ko_to_kegg = TRUE),
    "ko_to_kegg = TRUE requires pathway = 'KO'"
  )
  expect_error(
    pathway_annotation(pathway = "MetaCyc", daa_results_df = test_df, ko_to_kegg = TRUE),
    "ko_to_kegg = TRUE requires pathway = 'KO'"
  )
})

test_that("pathway_annotation allows ko_to_kegg = TRUE with pathway = 'KO' or NULL", {
  # Both forms are the legitimate KO-to-KEGG use case. We only exercise input
  # validation here (network calls are gated by GGPICRUST2_RUN_NETWORK_TESTS),
  # so wrap the network hop in tryCatch and assert we get past the parameter
  # check.
  test_df <- data.frame(
    feature = c("K00001"),
    p_adjust = c(0.01),
    stringsAsFactors = FALSE
  )

  past_check <- function(expr) {
    result <- tryCatch(expr, error = function(e) e)
    # If validation rejected the combination it would raise the KO-only
    # message; any other error means we made it past the boundary check.
    if (inherits(result, "error")) {
      expect_false(grepl("ko_to_kegg = TRUE requires pathway = 'KO'",
                         conditionMessage(result), fixed = TRUE))
    }
  }

  past_check(suppressMessages(
    pathway_annotation(pathway = "KO", daa_results_df = test_df, ko_to_kegg = TRUE)
  ))
  past_check(suppressMessages(
    pathway_annotation(pathway = NULL, daa_results_df = test_df, ko_to_kegg = TRUE)
  ))
})
