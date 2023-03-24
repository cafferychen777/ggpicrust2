#' Perform Principal Component Analysis (PCA) on functional pathway abundance data and create visualizations of the PCA results.
#'
#' @param abundance A data frame, predicted functional pathway abundance.
#' @param metadata A tibble, consisting of sample information.
#' @param group A character, group name.
#' @return A ggplot object showing the PCA results.
#' @export
#'
#' @examples
#' # Create example functional pathway abundance data
#' abundance_example <- data.frame(A = rnorm(10), B = rnorm(10), C = rnorm(10))
#'
#' # Create example metadata
#' metadata_example <- tibble::tibble(sample_id = 1:10,
#'                                    group = factor(rep(c("Control", "Treatment"), each = 5)))
#'
#' # Perform PCA and create visualizations
#' pca_plot <- pathway_pca(t(abundance_example), metadata_example, "group")
#' print(pca_plot)
pathway_pca <- function(abundance, metadata, group){
  # due to NSE notes in R CMD check
  PC1 = PC2 = Group = NULL
  # Perform PCA on the abundance data, keeping the first two principal components
  pca_axis <- stats::prcomp(t(abundance), center = TRUE, scale = TRUE)$x[,1:2]

  # Calculate the proportion of total variation explained by each PC
  pca_proportion <- stats::prcomp(t(abundance), center = TRUE, scale = TRUE)$sdev[1:2]/sum(stats::prcomp(t(abundance), center = TRUE, scale = TRUE)$sdev)*100

  # Combine the PCA results with the metadata information
  pca <- cbind(pca_axis, metadata[,group])
  pca$Group <- pca[,group]

  levels <- length(levels(factor(pca$Group)))



  colors_choices <- c("#d93c3e", "#3685bc","#208A42","#89288F","#F47D2B",
             "#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4",
             "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
             "#6E4B9E","#0C727C", "#7E1416","#D8A767","#3D3D3D")[1:levels]


  # Create a ggplot object for the PCA scatter plot
  Fig1a.taxa.pca <- ggplot2::ggplot(pca,ggplot2::aes(PC1,PC2))+
    ggplot2::geom_point(size=2,ggplot2::aes(color=Group),show.legend = T)+
    ggplot2::scale_color_manual(values=colors_choices)+
    ggplot2::stat_ellipse(ggplot2::aes(color = Group),fill="white",geom = "polygon",
                 level=0.95,alpha = 0.01,show.legend = F)+
    ggplot2::labs(x=paste0("PC1(",round(pca_proportion[1],1),"%)"),y=paste0("PC1(",round(pca_proportion[2],1),"%)"),color = group)+
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
    ggplot2::scale_fill_manual(values=colors_choices) + # Manually set the fill color for each group
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
    ggplot2::scale_fill_manual(values=colors_choices) + # Manually set the fill color for each group
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
