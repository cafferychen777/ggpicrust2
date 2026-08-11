test_that("create_gradient_colors returns exactly the requested size", {
  for (n_colors in 1:12) {
    colors <- create_gradient_colors(n_colors = n_colors)
    expect_length(colors, n_colors)
    expect_false(anyNA(colors))
    expect_true(all(vapply(
      colors,
      function(color) !inherits(try(grDevices::col2rgb(color), silent = TRUE),
                                  "try-error"),
      logical(1)
    )))
  }
})

test_that("create_gradient_colors preserves diverging endpoints", {
  theme <- get_color_theme("default")

  for (n_colors in c(2, 4, 5, 6, 11)) {
    colors <- create_gradient_colors("default", n_colors = n_colors)
    expect_equal(tolower(colors[1]),
                 tolower(unname(theme$fold_change_colors["negative"])))
    expect_equal(tolower(colors[n_colors]),
                 tolower(unname(theme$fold_change_colors["positive"])))
  }

  odd_colors <- create_gradient_colors("default", n_colors = 5)
  expect_equal(
    tolower(odd_colors[3]),
    tolower(unname(theme$fold_change_colors["neutral"]))
  )
})

test_that("color theme counts and gradient flags are validated", {
  expect_error(get_color_theme(n_colors = 0), "n_colors.*positive")
  expect_error(create_gradient_colors(n_colors = 2.5), "n_colors.*integer")
  expect_error(create_gradient_colors(diverging = NA), "diverging.*TRUE or FALSE")
  expect_error(
    get_color_theme("defualt"),
    "Unknown color theme 'defualt'.*Available themes"
  )
})

test_that("color themes expand categorical palettes without recycling colors", {
  theme <- get_color_theme("default", n_colors = 12)

  expect_length(theme$group_colors, 12)
  expect_length(unique(tolower(theme$group_colors)), 12)
  expect_true(all(vapply(
    theme$group_colors,
    function(color) !inherits(try(grDevices::col2rgb(color), silent = TRUE),
                                "try-error"),
    logical(1)
  )))
})

test_that("smart color selection validates every decision input", {
  expect_error(smart_color_selection(0), "n_groups.*positive")
  expect_error(
    smart_color_selection(2, has_pathway_class = NA),
    "has_pathway_class.*TRUE or FALSE"
  )
  expect_error(
    smart_color_selection(2, accessibility_mode = NA),
    "accessibility_mode.*TRUE or FALSE"
  )
  expect_error(
    smart_color_selection(2, data_type = "unknown"),
    "data_type.*must be one of"
  )

  selection <- smart_color_selection(
    2,
    has_pathway_class = TRUE,
    accessibility_mode = TRUE
  )
  expect_equal(selection$theme_name, "colorblind_friendly")
  expect_match(selection$reason, "pathway-class annotations")

  expect_equal(
    smart_color_selection(
      2,
      data_type = "pvalue",
      accessibility_mode = TRUE
    )$theme_name,
    "colorblind_friendly"
  )
  expect_equal(
    smart_color_selection(
      2,
      data_type = "foldchange",
      accessibility_mode = TRUE
    )$theme_name,
    "colorblind_friendly"
  )
})

test_that("preview_color_theme validates save controls before side effects", {
  expect_error(
    preview_color_theme(save_plot = NA),
    "save_plot.*TRUE or FALSE"
  )
  expect_error(
    preview_color_theme(save_plot = TRUE, filename = ""),
    "filename.*non-empty"
  )
  expect_s3_class(preview_color_theme(save_plot = "FALSE"), "ggplot")
})
