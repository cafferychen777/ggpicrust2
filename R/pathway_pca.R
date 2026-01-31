#' Perform Principal Component Analysis (PCA) on functional pathway abundance data
#'
#' This function performs PCA analysis on pathway abundance data and creates an informative visualization
#' that includes a scatter plot of the first two principal components (PC1 vs PC2) with density plots
#' for both PCs. The plot helps to visualize the clustering patterns and distribution of samples
#' across different groups.
#'
#' @param abundance A numeric matrix or data frame containing pathway abundance data.
#'        Rows represent pathways, columns represent samples.
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
#' supporting various sample identifier formats. Samples and pathways with zero
#' variance are filtered before PCA.
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
  validate_abundance(abundance, min_samples = 3)
  validate_metadata(metadata)
  validate_group(metadata, group, min_groups = 1)

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
  if (nrow(abundance) < 2) {
    stop("PCA requires at least 2 pathways (rows)")
  }

  # Filter out samples with zero variance
  sample_var <- apply(abundance, 2, var)
  zero_var_samples <- sample_var == 0
  if (any(zero_var_samples)) {
    warning(paste("Removing", sum(zero_var_samples), "sample(s) with zero variance"))
    abundance <- abundance[, !zero_var_samples, drop = FALSE]
    metadata <- metadata[!zero_var_samples, , drop = FALSE]
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

  # Convert group to factor if needed
  if (!is.factor(metadata[[group]])) {
    metadata[[group]] <- factor(metadata[[group]])
  }

  # Validate colors if provided
  if (!is.null(colors)) {
    n_groups <- length(levels(metadata[[group]]))
    if (!is.character(colors)) {
      stop("Colors must be a character vector")
    }
    if (!all(sapply(colors, function(x) tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE)))) {
      stop("Invalid color names provided")
    }
    if (length(colors) != n_groups) {
      stop(sprintf("Number of colors (%d) does not match number of groups (%d)",
                   length(colors), n_groups))
    }
  }

  # due to NSE notes in R CMD check
  PC1 = PC2 = Group = NULL
  
  # Perform PCA on the abundance data
  pca_result <- tryCatch({
    stats::prcomp(t(abundance), center = TRUE, scale = TRUE)
  }, error = function(e) {
    if (grepl("cannot rescale a constant/zero column", e$message)) {
      stop("PCA failed: Some samples or pathways have zero variance. ",
           "This has been checked earlier, but there might still be near-zero variance columns. ",
           "Consider using a more stringent filtering or transforming your data.")
    } else {
      stop("PCA failed: ", e$message)
    }
  })
  
  # Keep the first two principal components
  pca_axis <- pca_result$x[,1:2]

  # Calculate the proportion of total variance explained by each PC
  # Note: variance = sdev^2, so we need to square the standard deviations
  pca_proportion <- (pca_result$sdev[1:2]^2) / sum(pca_result$sdev^2) * 100

  # Combine the PCA results with the metadata information
  pca <- cbind(pca_axis, metadata %>% select(all_of(c(group))))
  pca$Group <- pca[,group]

  levels <- length(levels(factor(pca$Group)))

  # Set default colors if colors are not provided
  if (is.null(colors)) {
    colors <- c("#d93c3e", "#3685bc","#208A42","#89288F","#F47D2B",
                "#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4",
                "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
                "#6E4B9E","#0C727C", "#7E1416","#D8A767","#3D3D3D")[1:levels]
  }

  # Ensure the number of colors matches the number of levels in Group
  if (length(colors) != levels) {
    stop("The length of colors vector must match the number of levels in Group")
  }


  # Create a ggplot object for the PCA scatter plot
  Fig1a.taxa.pca <- ggplot2::ggplot(pca,ggplot2::aes(PC1,PC2))+
    ggplot2::geom_point(size=4,ggplot2::aes(color=Group),show.legend = T)+
    ggplot2::scale_color_manual(values=colors)+
    ggplot2::stat_ellipse(ggplot2::aes(color = Group),fill="white",geom = "polygon",
                 level=0.95,alpha = 0.01,show.legend = F)+
    ggplot2::labs(x=paste0("PC1(",round(pca_proportion[1],1),"%)"),y=paste0("PC2(",round(pca_proportion[2],1),"%)"),color = group)+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.line=ggplot2::element_line(colour = "black"),
          axis.title=ggplot2::element_text(color="black",face = "bold"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(color="black",size=10,face = "bold"),
          legend.text = ggplot2::element_text(size = 16, face = "bold"),
          legend.title = ggplot2::element_text(size = 16, face = "bold"))+
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black")

  # Create a ggplot object for the density plot of PC1
  # Plot the density of PC1 colored by Group
  Fig1a.taxa.pc1.density <-
    ggplot2::ggplot(pca) +
    ggplot2::geom_density(ggplot2::aes(x=PC1, group=Group, fill=Group), # Plot the density of PC1
                 color="black", alpha=1, position = 'identity',show.legend = F) +
    ggplot2::scale_fill_manual(values=colors) + # Manually set the fill color for each group
    ggplot2::theme_classic() + # Set the classic theme
    ggplot2::scale_y_discrete(expand = c(0,0.001)) + # Scale the y-axis with a small expansion to improve appearance
    ggplot2::labs(x=NULL, y=NULL) + # Remove x and y labels
    ggplot2::theme(axis.text.x=ggplot2::element_blank(), # Remove x axis text
          axis.ticks.x = ggplot2::element_blank()) # Remove x axis ticks

  # Plot the density of PC2 colored by Group
  Fig1a.taxa.pc2.density <-
    ggplot2::ggplot(pca) +
    ggplot2::geom_density(ggplot2::aes(x=PC2, group=Group, fill=Group), # Plot the density of PC2
                 color="black", alpha=1, position = 'identity',show.legend = F) +
    ggplot2::scale_fill_manual(values=colors) + # Manually set the fill color for each group
    ggplot2::theme_classic() + # Set the classic theme
    ggplot2::scale_y_discrete(expand = c(0,0.001)) + # Scale the y-axis with a small expansion to improve appearance
    ggplot2::labs(x=NULL, y=NULL) + # Remove x and y labels
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), # Remove y axis text
          axis.ticks.y = ggplot2::element_blank()) + # Remove y axis ticks
    ggplot2::coord_flip() # Flip the plot to change it from a horizontal to a vertical orientation

  # Return plot with or without marginal density plots
  if (show_marginal) {
    Fig1a.taxa.pca %>%
      aplot::insert_top(Fig1a.taxa.pc1.density, height = 0.3) %>%
      aplot::insert_right(Fig1a.taxa.pc2.density, width=0.3) %>%
      ggplotify::as.ggplot()
  } else {
    Fig1a.taxa.pca
  }
}
