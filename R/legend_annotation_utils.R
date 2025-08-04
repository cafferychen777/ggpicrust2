#' Legend and Annotation Utilities for ggpicrust2
#'
#' This module provides enhanced legend and annotation functionality for ggpicrust2
#' visualizations, including intelligent p-value formatting, significance marking,
#' and customizable legend styling.
#'
#' @name legend_annotation_utils
NULL

#' Smart P-value Formatting
#'
#' @param p_values Numeric vector of p-values
#' @param format Character string specifying format type
#' @param stars Logical, whether to include star symbols
#' @param thresholds Numeric vector of significance thresholds
#' @param star_symbols Character vector of star symbols
#' @return Character vector of formatted p-values
#' @export
format_pvalue_smart <- function(p_values, 
                               format = "smart", 
                               stars = TRUE,
                               thresholds = c(0.001, 0.01, 0.05),
                               star_symbols = c("***", "**", "*")) {
  
  # Input validation
  if (!is.numeric(p_values)) {
    stop("p_values must be numeric")
  }
  
  format <- match.arg(format, c("numeric", "scientific", "smart", "stars_only", "combined"))
  
  # Base formatting
  formatted <- switch(format,
    "numeric" = sprintf("%.3f", p_values),
    "scientific" = sprintf("%.2e", p_values),
    "smart" = ifelse(p_values < 0.001, 
                    sprintf("p < 0.001"), 
                    ifelse(p_values < 0.01,
                          sprintf("p = %.3f", p_values),
                          sprintf("p = %.2f", p_values))),
    "stars_only" = rep("", length(p_values)),
    "combined" = sprintf("%.3f", p_values)
  )
  
  # Add significance stars if requested
  if (stars && format != "numeric") {
    stars_vec <- get_significance_stars(p_values, thresholds, star_symbols)
    
    if (format == "stars_only") {
      formatted <- stars_vec
    } else if (format == "combined") {
      formatted <- paste(formatted, stars_vec)
    } else {
      # For smart format, add stars after the p-value
      formatted <- ifelse(stars_vec != "", 
                         paste(formatted, stars_vec), 
                         formatted)
    }
  }
  
  return(formatted)
}

#' Get Significance Stars
#'
#' @param p_values Numeric vector of p-values
#' @param thresholds Numeric vector of significance thresholds
#' @param symbols Character vector of star symbols
#' @return Character vector of star symbols
#' @export
get_significance_stars <- function(p_values, 
                                 thresholds = c(0.001, 0.01, 0.05),
                                 symbols = c("***", "**", "*")) {
  
  if (length(thresholds) != length(symbols)) {
    stop("thresholds and symbols must have the same length")
  }
  
  stars <- character(length(p_values))
  
  for (i in seq_along(thresholds)) {
    stars[p_values < thresholds[i]] <- symbols[i]
  }
  
  return(stars)
}

#' Get Significance Colors
#'
#' @param p_values Numeric vector of p-values
#' @param thresholds Numeric vector of significance thresholds  
#' @param colors Character vector of colors for each significance level
#' @param default_color Default color for non-significant values
#' @return Character vector of colors
#' @export
get_significance_colors <- function(p_values,
                                  thresholds = c(0.001, 0.01, 0.05),
                                  colors = c("#d73027", "#fc8d59", "#fee08b"),
                                  default_color = "#999999") {
  
  if (length(thresholds) != length(colors)) {
    stop("thresholds and colors must have the same length")
  }
  
  result_colors <- rep(default_color, length(p_values))
  
  for (i in seq_along(thresholds)) {
    result_colors[p_values < thresholds[i]] <- colors[i]
  }
  
  return(result_colors)
}

#' Create Enhanced Legend Theme
#'
#' @param position Legend position ("top", "bottom", "left", "right", "none")
#' @param direction Legend direction ("horizontal", "vertical")
#' @param title Legend title
#' @param title_size Title font size
#' @param text_size Text font size
#' @param key_size Key size in cm
#' @param key_width Key width
#' @param key_height Key height
#' @param ncol Number of columns
#' @param nrow Number of rows
#' @param box_just Legend box justification
#' @param margin Legend margin
#' @return ggplot2 theme elements
#' @export
create_legend_theme <- function(position = "top",
                               direction = "horizontal", 
                               title = NULL,
                               title_size = 12,
                               text_size = 10,
                               key_size = 0.8,
                               key_width = NULL,
                               key_height = NULL,
                               ncol = NULL,
                               nrow = NULL,
                               box_just = "center",
                               margin = ggplot2::margin(0, 0, 0, 0)) {
  
  # Validate inputs
  position <- match.arg(position, c("top", "bottom", "left", "right", "none"))
  direction <- match.arg(direction, c("horizontal", "vertical"))
  box_just <- match.arg(box_just, c("center", "top", "bottom", "left", "right"))
  
  # Auto-adjust direction based on position if not explicitly set
  if (missing(direction)) {
    direction <- ifelse(position %in% c("top", "bottom"), "horizontal", "vertical")
  }
  
  # Create legend theme
  legend_theme <- ggplot2::theme(
    legend.position = position,
    legend.direction = direction,
    legend.title = if (!is.null(title)) {
      ggplot2::element_text(size = title_size, face = "bold")
    } else {
      ggplot2::element_blank()
    },
    legend.text = ggplot2::element_text(size = text_size),
    legend.key.size = ggplot2::unit(key_size, "cm"),
    legend.box.just = box_just,
    legend.margin = margin
  )
  
  # Add key dimensions if specified
  if (!is.null(key_width)) {
    legend_theme$legend.key.width <- ggplot2::unit(key_width, "cm")
  }
  if (!is.null(key_height)) {
    legend_theme$legend.key.height <- ggplot2::unit(key_height, "cm")
  }
  
  # Add guide specifications for multi-column/row layouts
  guide_specs <- list()
  if (!is.null(ncol)) {
    guide_specs$ncol <- ncol
  }
  if (!is.null(nrow)) {
    guide_specs$nrow <- nrow
  }
  
  return(list(theme = legend_theme, guide_specs = guide_specs))
}

#' Smart Text Size Calculator
#'
#' @param n_items Number of items to display
#' @param base_size Base text size
#' @param min_size Minimum text size
#' @param max_size Maximum text size
#' @return Calculated text size
#' @export
calculate_smart_text_size <- function(n_items, base_size = 10, min_size = 8, max_size = 14) {
  
  if (n_items <= 5) {
    size <- base_size
  } else if (n_items <= 10) {
    size <- base_size * 0.9
  } else if (n_items <= 20) {
    size <- base_size * 0.8
  } else {
    size <- base_size * 0.7
  }
  
  # Ensure within bounds
  size <- max(min_size, min(max_size, size))
  
  return(size)
}

#' Create Pathway Class Annotation Theme
#'
#' @param text_size Text size
#' @param text_color Text color
#' @param text_face Text face ("plain", "bold", "italic")
#' @param text_family Text family
#' @param text_angle Text angle in degrees
#' @param text_hjust Horizontal justification (0-1)
#' @param text_vjust Vertical justification (0-1)
#' @param bg_color Background color
#' @param bg_alpha Background alpha
#' @param position Annotation position ("left", "right", "none")
#' @return List of annotation styling parameters
#' @export
create_pathway_class_theme <- function(text_size = "auto",
                                     text_color = "black",
                                     text_face = "bold",
                                     text_family = "sans",
                                     text_angle = 0,
                                     text_hjust = 0.5,
                                     text_vjust = 0.5,
                                     bg_color = NULL,
                                     bg_alpha = 0.2,
                                     position = "left") {
  
  # Validate inputs
  text_face <- match.arg(text_face, c("plain", "bold", "italic"))
  position <- match.arg(position, c("left", "right", "none"))
  
  return(list(
    text_size = text_size,
    text_color = text_color,
    text_face = text_face,
    text_family = text_family,
    text_angle = text_angle,
    text_hjust = text_hjust,
    text_vjust = text_vjust,
    bg_color = bg_color,
    bg_alpha = bg_alpha,
    position = position
  ))
}

#' Detect and Resolve Annotation Overlaps
#'
#' @param labels Character vector of labels
#' @param positions Numeric vector of positions
#' @param min_distance Minimum distance between labels
#' @return Adjusted positions
#' @export
resolve_annotation_overlaps <- function(labels, positions, min_distance = 1) {
  
  if (length(labels) != length(positions)) {
    stop("labels and positions must have the same length")
  }
  
  if (length(positions) <= 1) {
    return(positions)
  }
  
  # Sort by position
  order_idx <- order(positions)
  sorted_positions <- positions[order_idx]
  
  # Adjust overlapping positions
  adjusted_positions <- sorted_positions
  
  for (i in 2:length(adjusted_positions)) {
    if (adjusted_positions[i] - adjusted_positions[i-1] < min_distance) {
      adjusted_positions[i] <- adjusted_positions[i-1] + min_distance
    }
  }
  
  # Restore original order
  result_positions <- numeric(length(positions))
  result_positions[order_idx] <- adjusted_positions
  
  return(result_positions)
}

# Export utility functions
utils::globalVariables(c("p_values", "thresholds", "symbols", "colors"))