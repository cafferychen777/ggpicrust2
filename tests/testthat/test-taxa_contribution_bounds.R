alloc <- function(sample, labels, values, function_id = "K00001") {
  data.frame(sample = sample, function_id = function_id, taxon_label = labels,
             contribution = values, stringsAsFactors = FALSE)
}

test_that("compare_taxa_contributions gives the documented bounds", {
  x <- alloc("S1", c("A", "B", "Unclassified"), c(5, 3, 2))
  y <- alloc("S1", c("A", "B"), c(5, 5))
  res <- compare_taxa_contributions(x, y)
  expect_equal(res$x_unassigned, 0.2)
  expect_equal(res$y_unassigned, 0)
  expect_equal(res$tv_lower, 0)
  expect_equal(res$tv_upper, 0.2)
  expect_equal(res$tv_unassigned_as_category, 0.2)

  same <- compare_taxa_contributions(y, y)
  expect_equal(unlist(same[, c("tv_lower", "tv_upper", "tv_unassigned_as_category")]),
               c(tv_lower = 0, tv_upper = 0, tv_unassigned_as_category = 0))

  disjoint <- compare_taxa_contributions(alloc("S1", "A", 1), alloc("S1", "B", 1))
  expect_equal(c(disjoint$tv_lower, disjoint$tv_upper), c(1, 1))
})

test_that("bounds bracket the collapsed distance on random allocations", {
  set.seed(31)
  for (i in 1:200) {
    labels <- c(sample(LETTERS[1:6], sample(1:5, 1)), "Unclassified")
    x <- alloc("S1", labels, rexp(length(labels)))
    labels_y <- c(sample(LETTERS[1:6], sample(1:5, 1)), "Unclassified")
    y <- alloc("S1", labels_y, rexp(length(labels_y)))
    res <- compare_taxa_contributions(x, y)
    expect_true(res$tv_lower <= res$tv_unassigned_as_category + 1e-12)
    expect_true(res$tv_unassigned_as_category <= res$tv_upper + 1e-12)
    expect_true(res$tv_lower >= 0 && res$tv_upper <= 1)
  }
})

test_that("compare_taxa_contributions handles zero totals, filters and lumped taxa", {
  x <- rbind(alloc("S1", c("A", "B"), c(1, 1)),
             alloc("S2", "A", 0),
             alloc("S1", "C", 10, function_id = "K00002"))
  y <- rbind(alloc("S1", c("A", "B"), c(1, 3)),
             alloc("S2", "A", 2))
  res <- compare_taxa_contributions(x, y, function_ids = "K00001")
  expect_equal(res$sample, c("S1", "S2"))
  expect_equal(res$tv_lower[1], 0.25)
  expect_true(is.na(res$tv_lower[2]))

  lumped <- alloc("S1", c("A", "Other"), c(1, 1))
  expect_error(compare_taxa_contributions(lumped, y), "Other")
  expect_error(compare_taxa_contributions(x, alloc("S9", "A", 1)), "share no samples")
})

test_that("aggregate_taxa_contributions keeps unassigned mass out of Other", {
  contrib <- data.frame(
    sample = "S1", function_id = "K00001",
    taxon = c("ASV1", "ASV2", "ASV3", "ASV4"),
    taxon_function_abun = c(10, 5, 3, 50)
  )
  taxonomy <- data.frame(
    Feature.ID = c("ASV1", "ASV2", "ASV3", "ASV4"),
    Taxon = c("k__Bacteria; g__Alpha", "k__Bacteria; g__Beta",
              "k__Bacteria; g__Gamma", "k__Bacteria")
  )
  agg <- aggregate_taxa_contributions(contrib, taxonomy = taxonomy,
                                      tax_level = "Genus", top_n = 1)
  totals <- stats::setNames(agg$contribution, agg$taxon_label)
  expect_equal(totals[["Unclassified"]], 50)
  expect_equal(totals[["Alpha"]], 10)
  expect_equal(totals[["Other"]], 8)
})
