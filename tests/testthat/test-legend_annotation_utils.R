test_that("p-value formatting helpers reject impossible probabilities and thresholds", {
  expect_error(
    format_pvalue_smart(c(0.01, 1.2)),
    "between 0 and 1"
  )
  expect_error(
    get_significance_stars(c(0.01, 0.02), thresholds = c(0.05, -0.01), symbols = c("*", "**")),
    "thresholds"
  )
  expect_error(
    get_significance_colors(c(0.01, 0.02), colors = c("red", "bad_color", "blue")),
    "valid R color"
  )
})

test_that("p-value formatting helpers handle NA p-values explicitly", {
  expect_equal(
    get_significance_stars(c(0.0005, 0.02, NA)),
    c("***", "*", "")
  )
  expect_equal(
    get_significance_colors(
      c(0.0005, 0.02, NA),
      colors = c("red", "orange", "yellow"),
      default_color = "black"
    ),
    c("red", "yellow", "black")
  )
})

test_that("significance mappings are independent of threshold order", {
  thresholds <- c(0.05, 0.001, 0.01)
  symbols <- c("*", "***", "**")
  colors <- c("gold", "red", "orange")
  p_values <- c(0.0005, 0.005, 0.02, 0.2)

  expect_equal(
    get_significance_stars(p_values, thresholds, symbols),
    c("***", "**", "*", "")
  )
  expect_equal(
    get_significance_colors(p_values, thresholds, colors, "grey"),
    c("red", "orange", "gold", "grey")
  )
})

test_that("format_pvalue_smart validates the stars flag", {
  expect_error(
    format_pvalue_smart(0.01, stars = NA),
    "stars.*TRUE or FALSE"
  )
})

test_that("smart text sizing validates scalar counts and coherent bounds", {
  expect_error(
    calculate_smart_text_size(NA_integer_),
    "n_items.*finite integer"
  )
  expect_error(
    calculate_smart_text_size(c(2, 3)),
    "n_items.*single finite integer"
  )
  expect_error(
    calculate_smart_text_size(2, min_size = 14, max_size = 8),
    "min_size.*less than or equal"
  )
  expect_equal(calculate_smart_text_size(30), 8)
})

test_that("annotation overlap resolution rejects undefined geometry", {
  expect_error(
    resolve_annotation_overlaps(c("a", "b"), c(1, NA_real_)),
    "positions.*finite numeric"
  )
  expect_error(
    resolve_annotation_overlaps(c("a", "b"), c(1, 2), min_distance = -1),
    "min_distance.*non-negative"
  )
  expect_error(
    resolve_annotation_overlaps(1:2, c(1, 2)),
    "labels.*character"
  )
  expect_equal(
    resolve_annotation_overlaps(c("a", "b"), c(1, 1.5), min_distance = 1),
    c(1, 2)
  )
})

test_that("legend and pathway annotation themes validate visual parameters", {
  expect_error(
    create_legend_theme(key_size = 0),
    "key_size.*positive"
  )
  expect_error(
    create_legend_theme(ncol = 1.5),
    "ncol.*finite integer"
  )
  expect_error(
    create_pathway_class_theme(text_color = "not-a-color"),
    "text_color.*invalid R color"
  )
  expect_error(
    create_pathway_class_theme(text_hjust = 1.5),
    "text_hjust.*range"
  )
  expect_error(
    create_pathway_class_theme(text_size = "large"),
    "text_size.*auto"
  )
})
