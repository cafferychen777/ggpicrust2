#' Create pathway heatmap
#'
#' This function creates a heatmap of the predicted functional pathway abundance data. The function first makes the abundance data relative, then converts the abundance data to a long format and orders the samples based on the environment information. The heatmap is then created using the `ggplot2` library. The color palette, appearance and the color bar of the heatmap can be customized using the `scale_fill_gradientn`, `theme` and `guides` functions respectively.
#'
#' @param abundance A matrix or data frame of pathway abundance data, where rows correspond to samples and columns correspond to pathways.
#' @param metadata A data frame of metadata, where each row corresponds to a sample and each column corresponds to a metadata variable.
#' @param group A character string specifying the column name in the metadata data frame that contains the group variable.
#'
#' @return A ggplot heatmap object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @import tidyr

pathway_heatmap <- function(abundance, metadata, group) {
  # Make the abundance matrix relative
  rel_abundance <- make_relative(abundance)

  # Convert the abundance matrix to a data frame
  rel_df <- as.data.frame(rel_abundance)

  # Order the samples based on the environment information
  ordered_metadata <- metadata[order(metadata[, group]),]
  order <- ordered_metadata$sample_name

  # Convert the abundance data frame to a long format
  long_df <- rel_df %>%
    tibble::rownames_to_column() %>%
    tidyr::pivot_longer(cols = -rowname,
                        names_to = "Sample",
                        values_to = "Value")

  # Set the order of the samples in the heatmap
  long_df$Sample <- factor(long_df$Sample, levels = order)

  # Create the heatmap using ggplot
  p <-
    ggplot2::ggplot(data = long_df,
                    mapping = ggplot2::aes(x = Sample, y = rowname, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = c("#4575b4",
                                                       "#f7f7f7",
                                                       "#d73027")) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::scale_y_discrete(expand = c(0, 0), position = "left") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    # Customize the appearance of the heatmap
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(
        color = "black",
        size = 10,
        face = "bold"
      ),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank()
    ) +
    # Add a color bar to the heatmap
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(
        direction = "vertical",
        reverse = F,
        barwidth = unit(0.6, "cm"),
        barheight = unit(18.5, "cm")
      )
    )

  return(p)
}
