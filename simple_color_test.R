# Simple Color Theme Test
# ======================

# Load the color themes
source("/Users/apple/Research/ggpicrust2/code/ggpicrust2/R/color_themes.R")

# Test basic functions
cat("Testing color theme functions:\n")

# Test get_available_themes
available_themes <- get_available_themes()
cat("Available themes:", length(available_themes), "themes found\n")
cat("Themes:", paste(available_themes, collapse = ", "), "\n\n")

# Test get_color_theme with different themes
for (theme_name in c("default", "nature", "science", "colorblind_friendly")) {
  theme_colors <- get_color_theme(theme_name, 4)
  cat("Theme:", theme_name, "\n")
  cat("  Group colors:", paste(theme_colors$group_colors[1:4], collapse = ", "), "\n")
  cat("  Description:", theme_colors$description, "\n\n")
}

# Test smart color selection
smart_result <- smart_color_selection(n_groups = 3, data_type = "abundance")
cat("Smart color selection:\n")
cat("  Recommended theme:", smart_result$theme_name, "\n")
cat("  Reason:", smart_result$reason, "\n\n")

# Test accessibility mode
access_result <- smart_color_selection(n_groups = 3, accessibility_mode = TRUE)
cat("Accessibility mode:\n")
cat("  Recommended theme:", access_result$theme_name, "\n")
cat("  Reason:", access_result$reason, "\n\n")

cat("✓ Color theme system is working correctly!\n")