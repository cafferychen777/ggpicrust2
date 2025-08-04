#' Color Theme System for ggpicrust2
#'
#' This module provides a comprehensive color theme system for ggpicrust2 visualizations,
#' including journal-specific themes, colorblind-friendly palettes, and intelligent
#' color selection based on data characteristics.
#'
#' @name color_themes
NULL

#' Get Available Color Themes
#' 
#' @return A character vector of available theme names
#' @export
get_available_themes <- function() {
  return(c("default", "nature", "science", "cell", "nejm", "lancet", 
           "colorblind_friendly", "viridis", "plasma", "minimal", 
           "high_contrast", "pastel", "bold"))
}

#' Get Color Theme
#'
#' @param theme_name Character string specifying the theme name
#' @param n_colors Integer specifying the number of colors needed
#' @return A list containing theme colors and settings
#' @export
get_color_theme <- function(theme_name = "default", n_colors = 8) {
  
  # Define base themes
  themes <- list(
    
    # Default theme (current ggpicrust2 colors)
    default = list(
      group_colors = c("#d93c3e", "#3685bc", "#6faa3e", "#e8a825", "#c973e6", "#ee6b3d", "#2db0a7", "#f25292"),
      pathway_class_colors = c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B",
                              "#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4",
                              "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
                              "#6E4B9E","#0C727C", "#7E1416","#D8A767","#3D3D3D"),
      fold_change_colors = c(negative = "#2166ac", neutral = "#f7f7f7", positive = "#b2182b"),
      fold_change_single = "#87ceeb",
      description = "Default ggpicrust2 color scheme"
    ),
    
    # Nature journal style
    nature = list(
      group_colors = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"),
      pathway_class_colors = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
                              "#B09C85", "#F39B7F", "#FFDC91", "#8491B4", "#91D1C2", "#7E6148", "#B2DF8A", "#33A02C",
                              "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
      fold_change_colors = c(negative = "#3C5488", neutral = "#F7F7F7", positive = "#E64B35"),
      fold_change_single = "#4DBBD5",
      description = "Nature journal inspired colors - professional and clear"
    ),
    
    # Science journal style  
    science = list(
      group_colors = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"),
      pathway_class_colors = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
                              "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94",
                              "#f7b6d3", "#c7c7c7", "#dbdb8d", "#9edae5"),
      fold_change_colors = c(negative = "#1f77b4", neutral = "#f0f0f0", positive = "#d62728"),
      fold_change_single = "#2ca02c",
      description = "Science journal style - classic and authoritative"
    ),
    
    # Cell journal style
    cell = list(
      group_colors = c("#c51b7d", "#4d9221", "#2166ac", "#762a83", "#f1a340", "#998ec3", "#1b7837", "#e08214"),
      pathway_class_colors = c("#c51b7d", "#4d9221", "#2166ac", "#762a83", "#f1a340", "#998ec3", "#1b7837", "#e08214",
                              "#f7f7f7", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837", "#00441b", "#f4a582", "#d6604d",
                              "#b2182b", "#67001f", "#053061", "#2166ac"),
      fold_change_colors = c(negative = "#2166ac", neutral = "#f7f7f7", positive = "#c51b7d"),
      fold_change_single = "#4d9221",
      description = "Cell journal style - vibrant and distinctive"
    ),
    
    # NEJM style
    nejm = list(
      group_colors = c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91", "#EE4C97"),
      pathway_class_colors = c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91", "#EE4C97",
                              "#B09C85", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B2DF8A", "#33A02C",
                              "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
      fold_change_colors = c(negative = "#0072B5", neutral = "#F5F5F5", positive = "#BC3C29"),
      fold_change_single = "#20854E",
      description = "New England Journal of Medicine style - medical professional"
    ),
    
    # Lancet style
    lancet = list(
      group_colors = c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6"),
      pathway_class_colors = c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6",
                              "#1B1919", "#FF6A6A", "#90EE90", "#87CEEB", "#DDA0DD", "#F0E68C", "#FFB6C1", "#D3D3D3",
                              "#8FBC8F", "#20B2AA", "#9370DB", "#CD853F"),
      fold_change_colors = c(negative = "#00468B", neutral = "#F8F8FF", positive = "#ED0000"),
      fold_change_single = "#42B540",
      description = "The Lancet style - authoritative medical colors"
    ),
    
    # Colorblind friendly (using ColorBrewer)
    colorblind_friendly = list(
      group_colors = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"),
      pathway_class_colors = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666",
                              "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                              "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"),
      fold_change_colors = c(negative = "#7570b3", neutral = "#f7f7f7", positive = "#d95f02"),
      fold_change_single = "#1b9e77",
      description = "Colorblind friendly palette - accessible to all users"
    ),
    
    # Viridis-inspired
    viridis = list(
      group_colors = c("#440154", "#31688e", "#35b779", "#fde725", "#482777", "#287d8e", "#95d840", "#c2df23"),
      pathway_class_colors = c("#440154", "#482777", "#3f4a8a", "#31688e", "#26838f", "#1f9d8a", "#6cce5a", "#b6de2b",
                              "#fee825", "#35b779", "#1fa187", "#4fbc73", "#85d54a", "#c2df23", "#f0f921", "#fcffa4",
                              "#2d6a8c", "#3d4a8a", "#482777", "#440154"),
      fold_change_colors = c(negative = "#440154", neutral = "#35b779", positive = "#fde725"),
      fold_change_single = "#31688e",
      description = "Viridis-inspired - perceptually uniform colors"
    ),
    
    # Plasma-inspired
    plasma = list(
      group_colors = c("#0d0887", "#6a00a8", "#b12a90", "#e16462", "#fca636", "#f0f921", "#9c179e", "#5302a3"),
      pathway_class_colors = c("#0d0887", "#2e049c", "#4c02a1", "#6a00a8", "#8b0aa5", "#a82296", "#c43c75", "#dd513a",
                              "#f1605d", "#feaa3a", "#feca57", "#f0f921", "#b12a90", "#e16462", "#fca636", "#fcffa4",
                              "#9c179e", "#5302a3", "#320597", "#0d0887"),
      fold_change_colors = c(negative = "#0d0887", neutral = "#dd513a", positive = "#f0f921"),
      fold_change_single = "#b12a90",
      description = "Plasma-inspired - vibrant and energetic"
    ),
    
    # Minimal/clean style
    minimal = list(
      group_colors = c("#2c3e50", "#e74c3c", "#3498db", "#2ecc71", "#f39c12", "#9b59b6", "#1abc9c", "#34495e"),
      pathway_class_colors = c("#2c3e50", "#34495e", "#3498db", "#2980b9", "#2ecc71", "#27ae60", "#e74c3c", "#c0392b",
                              "#f39c12", "#e67e22", "#9b59b6", "#8e44ad", "#1abc9c", "#16a085", "#f1c40f", "#95a5a6",
                              "#ecf0f1", "#bdc3c7", "#7f8c8d", "#2c3e50"),
      fold_change_colors = c(negative = "#3498db", neutral = "#ecf0f1", positive = "#e74c3c"),
      fold_change_single = "#2ecc71",
      description = "Minimal clean style - modern and professional"
    ),
    
    # High contrast
    high_contrast = list(
      group_colors = c("#000000", "#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928", "#A6CEE3"),
      pathway_class_colors = c("#000000", "#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928", "#A6CEE3",
                              "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FFFFFF", "#999999",
                              "#FB9A99", "#E31A1C", "#A6CEE3", "#1F78B4"),
      fold_change_colors = c(negative = "#1F78B4", neutral = "#FFFFFF", positive = "#E31A1C"),
      fold_change_single = "#33A02C",
      description = "High contrast - maximum visibility and clarity"
    ),
    
    # Pastel theme
    pastel = list(
      group_colors = c("#FFB3BA", "#BAFFC9", "#BAE1FF", "#FFFFBA", "#FFD1FF", "#E1BAFF", "#BAFFE1", "#FFBABA"),
      pathway_class_colors = c("#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFBA", "#BAFFFF", "#BABDFF", "#DBBAFF", "#FFBAD1",
                              "#FFC9BA", "#D1FFBA", "#BAFFE1", "#BAF7FF", "#BAD1FF", "#E1BAFF", "#FFBADB", "#F7FFBA",
                              "#C9BAFF", "#FFBAC9", "#BADBFF", "#BAFFBA"),
      fold_change_colors = c(negative = "#BAE1FF", neutral = "#FFFFBA", positive = "#FFB3BA"),
      fold_change_single = "#BAFFC9",
      description = "Pastel theme - soft and gentle colors"
    ),
    
    # Bold theme
    bold = list(
      group_colors = c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#FF8000", "#8000FF"),
      pathway_class_colors = c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#FF8000", "#8000FF",
                              "#FF4000", "#80FF00", "#0080FF", "#FF0080", "#8000FF", "#FF8040", "#40FF80", "#4080FF",
                              "#FF4080", "#80FF40", "#8040FF", "#FF0040"),
      fold_change_colors = c(negative = "#0000FF", neutral = "#808080", positive = "#FF0000"),
      fold_change_single = "#00FF00",
      description = "Bold theme - maximum impact and visibility"
    )
  )
  
  # Get the requested theme
  if (!theme_name %in% names(themes)) {
    warning(paste("Theme", theme_name, "not found. Using default theme."))
    theme_name <- "default"
  }
  
  theme <- themes[[theme_name]]
  
  # Adjust colors to requested number
  if (n_colors > 0) {
    theme$group_colors <- rep(theme$group_colors, length.out = n_colors)
    theme$pathway_class_colors <- rep(theme$pathway_class_colors, length.out = max(n_colors, 20))
  }
  
  return(theme)
}

#' Smart Color Selection
#'
#' Intelligently selects colors based on data characteristics
#'
#' @param n_groups Number of groups in the data
#' @param has_pathway_class Whether pathway class information is available
#' @param data_type Type of data ("abundance", "pvalue", "foldchange")
#' @param accessibility_mode Whether to use accessibility-friendly colors
#' @return A list with suggested theme and colors
#' @export
smart_color_selection <- function(n_groups, has_pathway_class = FALSE, 
                                 data_type = "abundance", accessibility_mode = FALSE) {
  
  # Base recommendation logic
  if (accessibility_mode) {
    recommended_theme <- "colorblind_friendly"
  } else if (n_groups <= 3) {
    recommended_theme <- "nature"  # Clean and professional for few groups
  } else if (n_groups <= 6) {
    recommended_theme <- "science"  # Good balance for medium number of groups
  } else {
    recommended_theme <- "viridis"  # Better for many groups
  }
  
  # Adjust based on data type
  if (data_type == "pvalue") {
    recommended_theme <- "high_contrast"  # Better for significance visualization
  } else if (data_type == "foldchange") {
    recommended_theme <- "cell"  # Good diverging colors for fold changes
  }
  
  # Get the theme
  theme <- get_color_theme(recommended_theme, n_groups)
  
  # Add recommendation explanation
  theme$recommendation_reason <- paste(
    "Selected", recommended_theme, "theme for", n_groups, "groups,",
    "data type:", data_type,
    if (accessibility_mode) ", with accessibility considerations" else ""
  )
  
  return(list(
    theme_name = recommended_theme,
    theme = theme,
    reason = theme$recommendation_reason
  ))
}

#' Create Gradient Colors
#'
#' Creates gradient colors for fold change visualization
#'
#' @param theme_name Character string specifying the theme
#' @param n_colors Number of colors in the gradient
#' @param diverging Whether to create a diverging gradient (for fold changes)
#' @return A vector of colors
#' @export
create_gradient_colors <- function(theme_name = "default", n_colors = 11, diverging = TRUE) {
  
  theme <- get_color_theme(theme_name)
  
  if (n_colors == 1) {
    # Special case for single color
    colors <- theme$group_colors[1]
  } else if (diverging) {
    # Create diverging gradient for fold changes
    low_color <- theme$fold_change_colors["negative"]
    mid_color <- theme$fold_change_colors["neutral"] 
    high_color <- theme$fold_change_colors["positive"]
    
    if (n_colors == 2) {
      colors <- c(low_color, high_color)
    } else {
      # Create gradient
      n_half <- (n_colors - 1) / 2
      low_gradient <- colorRampPalette(c(low_color, mid_color))(n_half + 1)
      high_gradient <- colorRampPalette(c(mid_color, high_color))(n_half + 1)
      
      colors <- c(low_gradient[1:n_half], mid_color, high_gradient[2:(n_half + 1)])
    }
  } else {
    # Create single-direction gradient
    colors <- colorRampPalette(theme$group_colors[1:2])(n_colors)
  }
  
  return(colors)
}

#' Preview Color Theme
#'
#' Creates a visual preview of a color theme
#'
#' @param theme_name Character string specifying the theme name
#' @param save_plot Whether to save the preview plot
#' @param filename Filename for saved plot
#' @return A ggplot object showing the color preview
#' @export
preview_color_theme <- function(theme_name = "default", save_plot = FALSE, filename = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for theme preview")
  }
  
  theme <- get_color_theme(theme_name, 8)
  
  # Create preview data
  preview_data <- data.frame(
    type = rep(c("Group Colors", "Pathway Class Colors", "Fold Change Colors"), 
               c(8, 20, 3)),
    color = c(theme$group_colors[1:8], 
              theme$pathway_class_colors[1:20],
              unlist(theme$fold_change_colors)),
    index = c(1:8, 1:20, 1:3),
    label = c(paste("Group", 1:8),
              paste("Class", 1:20), 
              c("Negative", "Neutral", "Positive"))
  )
  
  # Create the plot
  p <- ggplot2::ggplot(preview_data, ggplot2::aes(x = index, y = type, fill = color)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3, angle = 45) +
    ggplot2::labs(
      title = paste("Color Theme Preview:", theme_name),
      subtitle = theme$description,
      x = "Color Index",
      y = "Color Category"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12)
    )
  
  if (save_plot) {
    if (is.null(filename)) {
      filename <- paste0("color_theme_preview_", theme_name, ".pdf")
    }
    ggplot2::ggsave(filename, p, width = 12, height = 6)
    message(paste("Theme preview saved as", filename))
  }
  
  return(p)
}

# Export all functions
utils::globalVariables(c("type", "color", "index", "label"))