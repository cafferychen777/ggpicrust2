depth_fixture <- function() {
  set.seed(2026)
  x <- matrix(rexp(200 * 16), nrow = 200,
              dimnames = list(sprintf("K%05d", 1:200), paste0("s", 1:16)))
  depth <- setNames(round(exp(rnorm(16, log(2e4), 0.5))), colnames(x))
  metadata <- data.frame(sample = colnames(x),
                         group = rep(c("A", "B"), each = 8),
                         z = rnorm(16),
                         reads = unname(depth))
  sets <- split(rownames(x), rep(paste0("set", 1:10), each = 20))
  list(x = x, depth = depth, metadata = metadata, sets = sets)
}

test_that("sequencing_depth adds the centered log depth to the design", {
  skip_if_not_installed("limma")
  d <- depth_fixture()
  manual <- d$metadata
  manual$log_sequencing_depth <- log(d$depth) - mean(log(d$depth))
  for (method in c("camera", "fry")) {
    expected <- suppressMessages(pathway_gsea(d$x, manual, "group", method = method,
      covariates = c("z", "log_sequencing_depth"), gene_sets = d$sets))
    by_vector <- suppressMessages(pathway_gsea(d$x, d$metadata, "group", method = method,
      covariates = "z", sequencing_depth = d$depth, gene_sets = d$sets))
    by_column <- suppressMessages(pathway_gsea(d$x, d$metadata, "group", method = method,
      covariates = "z", sequencing_depth = "reads", gene_sets = d$sets))
    expect_equal(by_vector, expected)
    expect_equal(by_column, expected)
    # Vector input is matched by sample ID, not by position
    shuffled <- suppressMessages(pathway_gsea(d$x, d$metadata, "group", method = method,
      covariates = "z", sequencing_depth = rev(d$depth), gene_sets = d$sets))
    expect_equal(shuffled, expected)
    unadjusted <- suppressMessages(pathway_gsea(d$x, d$metadata, "group", method = method,
      covariates = "z", gene_sets = d$sets))
    expect_false(isTRUE(all.equal(unadjusted$pvalue, expected$pvalue)))
  }
})

test_that("sequencing_depth is validated", {
  skip_if_not_installed("limma")
  d <- depth_fixture()
  run <- function(...) suppressMessages(pathway_gsea(d$x, d$metadata, "group",
                                                      gene_sets = d$sets, ...))
  expect_error(run(sequencing_depth = unname(d$depth)), "named by unique sample IDs")
  expect_error(run(sequencing_depth = d$depth[-1]), "no value for 1 sample")
  bad <- d$depth; bad[3] <- 0
  expect_error(run(sequencing_depth = bad), "positive, finite")
  expect_error(run(sequencing_depth = setNames(rep(1e4, 16), colnames(d$x))),
               "constant across samples")
  expect_error(run(sequencing_depth = "missing_column"), "not found in metadata")
  expect_error(run(sequencing_depth = "reads", covariates = "reads"),
               "must not also be the group or a covariate")
  expect_error(run(sequencing_depth = TRUE), "metadata column name or a named numeric")
  clash <- d$metadata
  clash$log_sequencing_depth <- 1
  expect_error(suppressMessages(pathway_gsea(d$x, clash, "group", gene_sets = d$sets,
                                             sequencing_depth = d$depth)),
               "already contains a column named 'log_sequencing_depth'")
})

test_that("sequencing_depth is rejected for preranked methods", {
  d <- depth_fixture()
  expect_error(
    pathway_gsea(d$x, d$metadata, "group", method = "fgsea", gene_sets = d$sets,
                 sequencing_depth = d$depth),
    "only supported for 'camera' and 'fry'"
  )
})
