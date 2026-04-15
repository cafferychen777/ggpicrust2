# Tests for ggpicrust2() return structure

test_that("ggpicrust2 returns all expected fields", {
  skip_if(
    Sys.getenv("GGPICRUST2_RUN_E2E_TESTS", "false") != "true",
    "Set GGPICRUST2_RUN_E2E_TESTS=true to run full ggpicrust2 end-to-end tests."
  )
  skip_on_cran()
  skip_if_not_installed("ALDEx2")

  data(ko_abundance, package = "ggpicrust2")
  data(metadata, package = "ggpicrust2")

  result <- suppressMessages(ggpicrust2(
    data = ko_abundance,
    metadata = metadata,
    group = "Environment",
    pathway = "KO",
    daa_method = "ALDEx2",
    ko_to_kegg = TRUE
  ))

  # Check structure
  expect_type(result, "list")
  expect_true(all(c("abundance", "metadata", "group", "daa_results_df", "ko_to_kegg") %in% names(result)))

  # Check field types
  expect_true(is.data.frame(result$abundance) || is.matrix(result$abundance))
  expect_true(is.data.frame(result$metadata))
  expect_true(is.data.frame(result$daa_results_df))
  expect_equal(result$group, "Environment")
  expect_true(result$ko_to_kegg)

  # Backward compatibility: old API returned unnamed list elements accessible
  # by position. result[[1]] is the first DA method's plot+results.
  expect_true(!is.null(result[[1]]))
  expect_true(all(c("plot", "results") %in% names(result[[1]])))
})

test_that("ggpicrust2 aligns Group vector to abundance sample order before plotting", {
  captured_group <- NULL

  mock_abundance <- data.frame(
    S_B = c(10, 20),
    S_A = c(30, 40),
    row.names = c("K00001", "K00002"),
    check.names = FALSE
  )

  metadata <- data.frame(
    sample = c("S_A", "S_B"),
    Environment = c("A", "B"),
    stringsAsFactors = FALSE
  )

  mock_pathway_daa <- function(abundance, metadata, group, daa_method, select, p_adjust_method, reference) {
    data.frame(
      feature = rownames(abundance)[1],
      method = "mock_method",
      p_values = 0.01,
      adj_method = "BH",
      p_adjust = 0.01,
      group1 = "A",
      group2 = "B",
      stringsAsFactors = FALSE
    )
  }

  mock_pathway_annotation <- function(pathway, ko_to_kegg, daa_results_df, p_adjust_threshold) {
    daa_results_df$pathway_name <- daa_results_df$feature
    daa_results_df$pathway_class <- "Class1"
    daa_results_df
  }

  mock_pathway_errorbar <- function(abundance, daa_results_df, Group, ko_to_kegg, p_value_bar, order, colors, select, x_lab, p_values_threshold) {
    captured_group <<- Group
    NULL
  }

  result <- testthat::with_mocked_bindings(
    ggpicrust2(
      data = data.frame(function. = rownames(mock_abundance), mock_abundance, check.names = FALSE),
      metadata = metadata,
      group = "Environment",
      pathway = "KO",
      daa_method = "ALDEx2",
      ko_to_kegg = FALSE
    ),
    pathway_daa = mock_pathway_daa,
    pathway_annotation = mock_pathway_annotation,
    pathway_errorbar = mock_pathway_errorbar
  )

  expect_type(result, "list")
  expect_equal(names(captured_group), c("S_B", "S_A"))
  expect_equal(unname(captured_group), c("B", "A"))
})

test_that("ggpicrust2 rejects inconsistent ko_to_kegg / pathway combinations", {
  # Catch the user-intent conflict up front: ko_to_kegg = TRUE aggregates KO
  # counts into KEGG pathways, so it is only meaningful when pathway = "KO".
  # Previously this combination produced a cryptic "No features in abundance
  # data" error after the matrix collapsed to zero rows downstream.
  minimal_data <- data.frame(
    function. = c("1.1.1.1", "1.1.1.2"),
    S1 = c(10, 20),
    S2 = c(15, 25),
    stringsAsFactors = FALSE
  )
  minimal_metadata <- data.frame(
    sample = c("S1", "S2"),
    Condition = c("A", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    ggpicrust2(
      data = minimal_data,
      metadata = minimal_metadata,
      group = "Condition",
      pathway = "EC",
      daa_method = "LinDA",
      ko_to_kegg = TRUE
    ),
    "requires .pathway = .KO."
  )

  expect_error(
    ggpicrust2(
      data = minimal_data,
      metadata = minimal_metadata,
      group = "Condition",
      pathway = "MetaCyc",
      daa_method = "LinDA",
      ko_to_kegg = TRUE
    ),
    "requires .pathway = .KO."
  )
})

test_that("ggpicrust2 does not emit p.adjust deprecation warning on a normal call", {
  # Regression: ggpicrust2() had been renamed to accept `p_adjust_method`,
  # but its internal pathway_daa() call still passed `p.adjust = p.adjust`.
  # That made every normal ggpicrust2() call trigger the deprecation
  # warning that is supposed to fire only when a user explicitly supplies
  # the legacy argument.
  skip_if_not_installed("MicrobiomeStat")

  set.seed(42)
  n_ko <- 15; n_samp <- 8
  abund <- matrix(rpois(n_ko * n_samp, lambda = 50),
                  nrow = n_ko, ncol = n_samp)
  rownames(abund) <- paste0("K", sprintf("%05d", seq_len(n_ko)))
  colnames(abund) <- paste0("S", seq_len(n_samp))
  abund_df <- cbind(
    `function` = rownames(abund),
    as.data.frame(abund)
  )
  meta <- data.frame(
    sample_name = colnames(abund),
    Env = rep(c("A", "B"), each = n_samp / 2),
    stringsAsFactors = FALSE
  )

  dep_msgs <- character(0)
  tryCatch(
    withCallingHandlers(
      suppressMessages(ggpicrust2(
        data = abund_df,
        metadata = meta,
        group = "Env",
        pathway = "KO",
        daa_method = "LinDA",
        ko_to_kegg = FALSE
      )),
      warning = function(w) {
        if (grepl("p.adjust.*deprecated", conditionMessage(w))) {
          dep_msgs <<- c(dep_msgs, conditionMessage(w))
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) invisible(NULL)   # unrelated downstream errors OK
  )
  expect_length(dep_msgs, 0)
})

test_that("ggpicrust2 still warns when caller passes legacy p.adjust argument", {
  skip_if_not_installed("MicrobiomeStat")

  set.seed(42)
  n_ko <- 10; n_samp <- 6
  abund <- matrix(rpois(n_ko * n_samp, lambda = 50),
                  nrow = n_ko, ncol = n_samp)
  rownames(abund) <- paste0("K", sprintf("%05d", seq_len(n_ko)))
  colnames(abund) <- paste0("S", seq_len(n_samp))
  abund_df <- cbind(
    `function` = rownames(abund),
    as.data.frame(abund)
  )
  meta <- data.frame(
    sample_name = colnames(abund),
    Env = rep(c("A", "B"), each = n_samp / 2),
    stringsAsFactors = FALSE
  )

  got_dep <- FALSE
  tryCatch(
    withCallingHandlers(
      suppressMessages(ggpicrust2(
        data = abund_df,
        metadata = meta,
        group = "Env",
        pathway = "KO",
        daa_method = "LinDA",
        ko_to_kegg = FALSE,
        p.adjust = "BH"
      )),
      warning = function(w) {
        if (grepl("p.adjust.*deprecated", conditionMessage(w))) {
          got_dep <<- TRUE
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) invisible(NULL)
  )
  expect_true(got_dep)
})
