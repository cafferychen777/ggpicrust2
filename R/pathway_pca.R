#' Perform Principal Component Analysis (PCA) on functional pathway abundance data
#'
#' This function performs PCA analysis on pathway abundance data and creates an informative visualization
#' that includes a scatter plot of the first two principal components (PC1 vs PC2) with density plots
#' for both PCs. The plot helps to visualize the clustering patterns and distribution of samples
#' across different groups.
#'
#' @param abundance A numeric matrix or data frame containing pathway abundance data.
#'        Rows represent pathways, columns represent samples.
#'        Data frames may also provide a leading non-numeric feature ID column
#'        (for example \code{#NAME}, \code{feature}, or \code{pathway}); it is
#'        converted to row names before sample alignment.
#'        Column names must match the sample names in metadata.
#'        Values must be numeric and cannot contain missing values (NA).
#'
#' @param metadata A data frame containing sample information.
#'        Must include a column for grouping samples (specified by the 'group' parameter).
#'        Sample identifiers are auto-detected from columns named sample_name, Sample_ID,
#'        SampleID, etc., or from rownames.
#'
#' @param group A character string specifying the column name in metadata that contains
#'        group information for samples (e.g., "treatment", "condition", "group").
#'        Every sample retained after abundance/metadata alignment must have a
#'        non-missing, non-empty group label, and at least two groups must
#'        remain.
#'
#' @param colors Optional. A character vector of colors for different groups.
#'        Length must match the number of unique groups.
#'        If NULL, default colors will be used.
#'
#' @param show_marginal Logical. Whether to show marginal density plots for PC1 and PC2.
#'        Default is TRUE. Set to FALSE to show only the PCA scatter plot.
#'
#' @return A ggplot object showing:
#'        \itemize{
#'          \item Center: PCA scatter plot with confidence ellipses (95%)
#'          \item Top: Density plot for PC1
#'          \item Right: Density plot for PC2
#'        }
#'
#' @details
#' The function automatically aligns samples between abundance data and metadata,
#' supporting various sample identifier formats. Pathways with zero variance
#' across samples are filtered before PCA because they are variables in
#' \code{prcomp(t(abundance))} and cannot be scaled. Sample profiles are kept
#' as observations, even if a sample has zero variance across pathways. PCA
#' confidence ellipses are drawn only for groups with at least four samples;
#' smaller groups remain in the scatter plot but are skipped for ellipse
#' estimation. Marginal densities require at least two samples per group;
#' singleton groups are omitted from those panels, and the scatter plot is
#' returned alone when no group has enough observations for density estimation.
#'
#' @examples
#' # Create example abundance data
#' abundance_data <- matrix(rnorm(30), nrow = 3, ncol = 10)
#' colnames(abundance_data) <- paste0("Sample", 1:10)
#' rownames(abundance_data) <- c("PathwayA", "PathwayB", "PathwayC")
#'
#' # Create example metadata
#' metadata <- data.frame(
#'   sample_name = paste0("Sample", 1:10),
#'   group = factor(rep(c("Control", "Treatment"), each = 5))
#' )
#'
#' # Basic PCA plot with default colors
#' pca_plot <- pathway_pca(abundance_data, metadata, "group")
#'
#' # PCA plot with custom colors
#' pca_plot <- pathway_pca(
#'   abundance_data,
#'   metadata,
#'   "group",
#'   colors = c("blue", "red")  # One color per group
#' )
#'
#' # PCA plot without marginal density plots
#' pca_plot <- pathway_pca(
#'   abundance_data,
#'   metadata,
#'   "group",
#'   show_marginal = FALSE
#' )
#'
#' \donttest{
#' # Example with real data
#' data("metacyc_abundance")  # Load example pathway abundance data
#' data("metadata")          # Load example metadata
#'
#' # Generate PCA plot
#' # Prepare abundance data
#' abundance_data <- as.data.frame(metacyc_abundance)
#' rownames(abundance_data) <- abundance_data$pathway
#' abundance_data <- abundance_data[, -which(names(abundance_data) == "pathway")]
#' 
#' # Create PCA plot
#' pathway_pca(
#'   abundance_data,
#'   metadata,
#'   "Environment",
#'   colors = c("green", "purple")
#' )
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual stat_ellipse
#'             labs theme_classic theme element_line element_text element_blank
#'             geom_vline geom_hline geom_density scale_fill_manual coord_flip
#' @importFrom aplot insert_top insert_right
#' @importFrom ggplotify as.ggplot
#'
#' @export
pathway_pca <- function(abundance,
                        metadata,
                        group,
                        colors = NULL,
                        show_marginal = TRUE) {
  # Input validation using unified functions
  abundance <- normalize_abundance_feature_ids(
    abundance,
    context = "pathway_pca() abundance"
  )
  validate_abundance(abundance, min_samples = 3, check_zero_columns = FALSE)
  validate_metadata(metadata)
  validate_group(metadata, group, min_groups = 2)
  show_marginal <- normalize_logical_flag(show_marginal, "show_marginal")

  # Align samples between abundance and metadata
  aligned <- align_samples(abundance, metadata)
  abundance <- as.matrix(aligned$abundance)
  metadata <- aligned$metadata

  if (aligned$n_samples < 3) {
    stop(sprintf("PCA requires at least 3 matching samples, found %d", aligned$n_samples))
  }

  # PCA-specific: require complete numeric data
  if (!is.numeric(abundance)) {
    stop("Abundance data must contain only numeric values")
  }
  if (any(is.na(abundance))) {
    stop("Abundance matrix contains missing values (NA). PCA requires complete data.")
  }
  if (any(!is.finite(abundance))) {
    stop("Abundance matrix contains non-finite values. PCA requires finite data.")
  }
  if (nrow(abundance) < 2) {
    stop("PCA requires at least 2 pathways (rows)")
  }

  # Filter out pathways with zero variance
  pathway_var <- apply(abundance, 1, var)
  zero_var_pathways <- pathway_var == 0
  if (any(zero_var_pathways)) {
    warning(paste("Removing", sum(zero_var_pathways), "pathway(s) with zero variance"))
    abundance <- abundance[!zero_var_pathways, , drop = FALSE]
  }

  # Post-filter dimension checks
  if (nrow(abundance) < 2) {
    stop("After filtering, less than 2 pathways remain. PCA requires at least 2 pathways.")
  }
  if (ncol(abundance) < 3) {
    stop("After filtering, less than 3 samples remain. PCA requires at least 3 samples.")
  }

  group_values <- validate_group_vector_for_summary(
    metadata[[group]],
    context = paste0("metadata column '", group, "'"),
    sample_ids = colnames(abundance),
    min_groups = 2
  )
  metadata[[group]] <- droplevels(factor(group_values))
  n_groups <- nlevels(metadata[[group]])

  # Validate colors if provided
  if (!is.null(colors)) {
    colors <- normalize_level_colors(
      colors,
      levels(metadata[[group]]),
      "colors"
    )
  }

  # Perform PCA on the abundance data
  pca_result <- tryCatch({
    stats::prcomp(t(abundance), center = TRUE, scale = TRUE)
  }, error = function(e) {
    if (grepl("cannot rescale a constant/zero column", e$message)) {
      stop("PCA failed: Some pathways have zero variance. ",
           "This has been checked earlier, but there might still be near-zero variance columns. ",
           "Consider using a more stringent filtering or transforming your data.")
    } else {
      stop("PCA failed: ", e$message)
    }
  })
  
  # Keep the first two principal components
  pca_axis <- pca_result$x[, 1:2, drop = FALSE]

  # Calculate the proportion of total variance explained by each PC
  # Note: variance = sdev^2, so we need to square the standard deviations
  pca_proportion <- (pca_result$sdev[1:2]^2) / sum(pca_result$sdev^2) * 100

  # Keep internal plotting columns independent of user metadata names. A group
  # column named "PC1", "PC2", or "Group" previously created duplicate names
  # and failed only when ggplot evaluated its aesthetics.
  pca_data <- data.frame(
    PC1 = pca_axis[, 1],
    PC2 = pca_axis[, 2],
    Group = metadata[[group]],
    row.names = rownames(pca_axis),
    check.names = FALSE
  )

  # Set default colors if colors are not provided
  if (is.null(colors)) {
    base_colors <- c(
      "#d93c3e", "#3685bc", "#208A42", "#89288F", "#F47D2B",
      "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", "#90D5E4",
      "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8",
      "#6E4B9E", "#0C727C", "#7E1416", "#D8A767", "#3D3D3D"
    )
    colors <- if (n_groups <= length(base_colors)) {
      base_colors[seq_len(n_groups)]
    } else {
      grDevices::colorRampPalette(base_colors)(n_groups)
    }
    colors <- normalize_level_colors(
      colors,
      levels(metadata[[group]]),
      "colors"
    )
  }

  ellipse_counts <- table(pca_data$Group)
  ellipse_groups <- names(ellipse_counts)[ellipse_counts >= 4]
  skipped_ellipse_groups <- names(ellipse_counts)[ellipse_counts < 4]
  if (length(skipped_ellipse_groups) > 0) {
    warning(
      "Skipping PCA confidence ellipse(s) for group(s) with fewer than 4 samples: ",
      paste(paste0(skipped_ellipse_groups, "=", ellipse_counts[skipped_ellipse_groups]),
            collapse = ", "),
      ". ggplot2::stat_ellipse() requires at least 4 points for a two-dimensional ellipse.",
      call. = FALSE
    )
  }

  # Create a ggplot object for the PCA scatter plot
  pca_plot <- ggplot2::ggplot(
    pca_data,
    ggplot2::aes(x = .data$PC1, y = .data$PC2)
  ) +
    ggplot2::geom_point(
      size = 4,
      ggplot2::aes(color = .data$Group),
      show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      x = paste0("PC1(", round(pca_proportion[1], 1), "%)"),
      y = paste0("PC2(", round(pca_proportion[2], 1), "%)"),
      color = group
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black"),
      axis.title = ggplot2::element_text(color = "black", face = "bold"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(
        color = "black",
        size = 10,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(size = 16, face = "bold"),
      legend.title = ggplot2::element_text(size = 16, face = "bold")
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black")

  if (length(ellipse_groups) > 0) {
    ellipse_data <- pca_data[pca_data$Group %in% ellipse_groups, , drop = FALSE]
    pca_plot <- pca_plot +
      ggplot2::stat_ellipse(
        data = ellipse_data,
        ggplot2::aes(color = .data$Group),
        fill = "white",
        geom = "polygon",
        level = 0.95,
        alpha = 0.01,
        show.legend = FALSE
      )
  }

  if (!show_marginal) {
    return(pca_plot)
  }

  density_groups <- names(ellipse_counts)[ellipse_counts >= 2]
  skipped_density_groups <- names(ellipse_counts)[ellipse_counts < 2]
  if (length(skipped_density_groups) > 0) {
    warning(
      "Skipping PCA marginal density for group(s) with fewer than 2 samples: ",
      paste(
        paste0(skipped_density_groups, "=", ellipse_counts[skipped_density_groups]),
        collapse = ", "
      ),
      ". Density estimation requires at least 2 observations.",
      call. = FALSE
    )
  }
  if (length(density_groups) == 0) {
    return(pca_plot)
  }
  density_data <- pca_data[pca_data$Group %in% density_groups, , drop = FALSE]

  # Marginal density for PC1. geom_density() produces a continuous y
  # aesthetic, so its matching position scale must also be continuous.
  pc1_density <-
    ggplot2::ggplot(density_data) +
    ggplot2::geom_density(
      ggplot2::aes(
        x = .data$PC1,
        group = .data$Group,
        fill = .data$Group
      ),
      color = "black",
      alpha = 1,
      position = "identity",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = c(0, 0.001)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  # Marginal density for PC2 (horizontal after coord_flip()). Same
  # continuous-scale reasoning applies.
  pc2_density <-
    ggplot2::ggplot(density_data) +
    ggplot2::geom_density(
      ggplot2::aes(
        x = .data$PC2,
        group = .data$Group,
        fill = .data$Group
      ),
      color = "black",
      alpha = 1,
      position = "identity",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = c(0, 0.001)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    ) +
    ggplot2::coord_flip()

  pca_plot %>%
    aplot::insert_top(pc1_density, height = 0.3) %>%
    aplot::insert_right(pc2_density, width = 0.3) %>%
    ggplotify::as.ggplot()
}
