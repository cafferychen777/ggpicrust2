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
