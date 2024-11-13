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
#'        Must include:
#'        \itemize{
#'          \item A column named "sample_name" matching the column names in abundance
#'          \item A column for grouping samples (specified by the 'group' parameter)
#'        }
#'
#' @param group A character string specifying the column name in metadata that contains
#'        group information for samples (e.g., "treatment", "condition", "group").
#'
#' @param colors Optional. A character vector of colors for different groups.
#'        Length must match the number of unique groups.
#'        If NULL, default colors will be used.
#'
#' @return A ggplot object showing:
#'        \itemize{
#'          \item Center: PCA scatter plot with confidence ellipses (95%)
#'          \item Top: Density plot for PC1
#'          \item Right: Density plot for PC2
#'        }
#'
#' @details
#' The function performs several validations on input data:
#' \itemize{
#'   \item Abundance matrix must have at least 2 pathways and 3 samples
#'   \item All values in abundance matrix must be numeric
#'   \item Sample names must match between abundance and metadata
#'   \item Group column must exist in metadata
#'   \item If custom colors are provided, they must be valid color names or codes
#' }
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
#' \donttest{
#' # Example with real data
#' data("metacyc_abundance")  # Load example pathway abundance data
#' data("metadata")          # Load example metadata
#'
#' # Generate PCA plot
#' pathway_pca(
#'   metacyc_abundance %>% column_to_rownames("pathway"),
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
                        colors = NULL) {
  # Input validation
  # Check if inputs are missing
  if (missing(abundance)) {
    stop("Abundance matrix is required")
  }
  if (missing(metadata)) {
    stop("Metadata is required")
  }
  if (missing(group)) {
    stop("Group variable name is required")
  }
  
  # Check abundance matrix
  if (!is.matrix(abundance) && !is.data.frame(abundance)) {
    stop("Abundance must be a matrix or data frame")
  }
  if (any(is.na(abundance))) {
    stop("Abundance matrix contains missing values (NA)")
  }
  if (!all(apply(abundance, 2, is.numeric))) {
    stop("Abundance matrix must contain only numeric values")
  }
  if (nrow(abundance) < 2) {
    stop("Abundance matrix must contain at least 2 pathways")
  }
  if (ncol(abundance) < 3) {
    stop("Abundance matrix must contain at least 3 samples")
  }
  
  # Check metadata
  if (!is.data.frame(metadata)) {
    stop("Metadata must be a data frame")
  }
  if (!"sample_name" %in% colnames(metadata)) {
    stop("Metadata must contain a 'sample_name' column")
  }
  if (!group %in% colnames(metadata)) {
    stop(sprintf("Group column '%s' not found in metadata", group))
  }
  
  # Check sample names match
  if (length(unique(metadata$sample_name)) != ncol(abundance)) {
    stop("Number of unique samples in metadata does not match abundance matrix")
  }
  if (!all(colnames(abundance) %in% metadata$sample_name)) {
    stop("Some sample names in abundance matrix are not found in metadata")
  }
  
  # Check group variable
  group_values <- metadata[[group]]
  if (length(unique(group_values)) < 1) {
    stop("Group variable must have at least one level")
  }
  if (!is.factor(group_values)) {
    warning("Converting group variable to factor")
    metadata[[group]] <- factor(group_values)
  }
  
  # Check colors if provided
  if (!is.null(colors)) {
    if (!is.vector(colors) || !is.character(colors)) {
      stop("Colors must be provided as a character vector")
    }
    if (!all(sapply(colors, function(x) tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE)))) {
      stop("Invalid color names provided")
    }
    levels <- length(levels(factor(metadata[[group]])))
    if (length(colors) != levels) {
      stop(sprintf("Number of colors (%d) does not match number of groups (%d)", 
                  length(colors), levels))
    }
  }

  # due to NSE notes in R CMD check
  PC1 = PC2 = Group = NULL
  
  # Perform PCA on the abundance data, keeping the first two principal components
  pca_axis <- stats::prcomp(t(abundance), center = TRUE, scale = TRUE)$x[,1:2]

  # Calculate the proportion of total variation explained by each PC
  pca_proportion <- stats::prcomp(t(abundance), center = TRUE, scale = TRUE)$sdev[1:2]/sum(stats::prcomp(t(abundance), center = TRUE, scale = TRUE)$sdev)*100

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

  # Combine the two plots into a single plot with the PC1 plot on top and the PC2 plot on the right
  Fig1a.taxa.pca %>%
    aplot::insert_top(Fig1a.taxa.pc1.density, height = 0.3) %>% # Insert the PC1 plot on top with a height of 0.3
    aplot::insert_right(Fig1a.taxa.pc2.density, width=0.3) %>% # Insert the PC2 plot on the right with a width of 0.3
    ggplotify::as.ggplot() # Convert the combined plot to a ggplot object
}
