test_that("logCPM public analyses match direct limma and per-sample rescaling", {
  skip_if_not_installed("limma")
  set.seed(6172)
  x <- matrix(rexp(240 * 20), nrow = 240,
              dimnames = list(sprintf("K%05d", 1:240), paste0("s", 1:20)))
  metadata <- data.frame(sample_name = colnames(x), group = factor(rep(c("A", "B"), each = 10)),
                         z = rnorm(20), row.names = colnames(x))
  sets <- split(rownames(x), rep(paste0("set", 1:12), each = 20))
  local_mocked_bindings(prepare_gene_sets = function(...) stop("Unexpected bundled reference load"))
  y <- log2(1e6 * sweep(x, 2, colSums(x), "/") + 0.5)
  design <- build_design_matrix(metadata, "group", "z")
  contrast <- resolve_limma_contrast(design, metadata, "group", NULL)
  indices <- lapply(sets, function(ids) which(rownames(x) %in% ids))
  for (method in c("camera", "fry")) {
    actual <- suppressMessages(pathway_gsea(x, metadata, "group", method = method,
      covariates = "z", transformation = "logCPM", inter.gene.cor = NA_real_, gene_sets = sets))
    direct <- if (method == "camera") {
      limma::camera(y, indices, design, contrast = contrast, inter.gene.cor = NA_real_, trend.var = TRUE)
    } else {
      limma::fry(y, indices, design, contrast = contrast, trend = TRUE)
    }
    expected <- direct[actual$pathway_id, "PValue"]
    expect_equal(actual$pvalue, expected, tolerance = 1e-12)
    expect_equal(actual$p.adjust, p.adjust(expected, "BH"), tolerance = 1e-12)
    scaled <- suppressMessages(pathway_gsea(sweep(x, 2, 10^seq(-3, 3, length.out = 20), "*"),
      metadata, "group", method = method, covariates = "z", transformation = "logCPM",
      inter.gene.cor = NA_real_, gene_sets = sets))
    expect_equal(scaled$pvalue[match(actual$pathway_id, scaled$pathway_id)], actual$pvalue, tolerance = 1e-11)
    legacy <- suppressMessages(pathway_gsea(x, metadata, "group", method = method, covariates = "z", gene_sets = sets))
    voom_fit <- limma::voom(x, design, plot = FALSE)
    native_legacy <- if (method == "camera") {
      limma::camera(voom_fit, indices, design, contrast = contrast)
    } else {
      limma::fry(voom_fit, indices, design, contrast = contrast)
    }
    expect_equal(legacy$pvalue, native_legacy[legacy$pathway_id, "PValue"],
                 tolerance = 1e-12)
    explicit <- suppressMessages(pathway_gsea(x, metadata, "group", method = method,
      covariates = "z", transformation = "voom", gene_sets = sets))
    expect_equal(legacy, explicit)
  }
})

test_that("logCPM detects a specified signal and reverses direction with the contrast", {
  skip_if_not_installed("limma")
  set.seed(90210)
  x <- matrix(exp(rnorm(600 * 24, 3, 0.6)), 600,
              dimnames = list(sprintf("K%05d", 1:600), paste0("s", 1:24)))
  x[1:30, 13:24] <- x[1:30, 13:24] * 4
  metadata <- data.frame(sample = colnames(x), group = rep(c("A", "B"), each = 12))
  sets <- split(rownames(x), rep(paste0("set", 1:20), each = 30))
  for (method in c("camera", "fry")) {
    forward <- suppressMessages(pathway_gsea(x, metadata, "group", method = method,
      gene_sets = sets, transformation = "logCPM"))
    reverse <- suppressMessages(pathway_gsea(x, metadata, "group", method = method,
      gene_sets = sets, transformation = "logCPM", contrast = c(0, -1)))
    reverse <- reverse[match(forward$pathway_id, reverse$pathway_id), ]
    expect_identical(forward$direction[forward$pathway_id == "set1"], "Up")
    expect_lt(forward$p.adjust[forward$pathway_id == "set1"], 0.05)
    expect_equal(forward$pvalue, reverse$pvalue, tolerance = 1e-10)
    expect_true(all(forward$direction != reverse$direction))
  }
})

test_that("logCPM rejects zero-total samples and unsupported analysis modes", {
  skip_if_not_installed("limma")
  x <- matrix(1, 40, 6, dimnames = list(sprintf("K%05d", 1:40), paste0("s", 1:6)))
  x[, 1] <- 0
  metadata <- data.frame(sample_name = colnames(x), group = rep(c("A", "B"), each = 3))
  local_mocked_bindings(prepare_gene_sets = function(...) list(target = rownames(x)[1:20]))
  expect_error(suppressMessages(pathway_gsea(x, metadata, "group", transformation = "logCPM")),
               "total abundance of 0")
  expect_error(pathway_gsea(x + 1, metadata, "group", method = "fgsea", transformation = "logCPM"),
               "only to camera/fry")
  expect_error(pathway_gsea(x, metadata, "group", transformation = "invalid"), "arg")
  expect_error(pathway_gsea(x + 1, metadata, "group", gene_sets = list(target = rownames(x)),
                            go_category = "all"), "cannot be combined")
  expect_error(suppressMessages(pathway_gsea(x + 1, metadata, "group", gene_sets = list(rownames(x)))),
               "nam")
})
