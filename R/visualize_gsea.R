#' Visualize GSEA results
#'
#' This function creates various visualizations for Gene Set Enrichment Analysis (GSEA) results.
#'
#' @param gsea_results A data frame containing GSEA results from the pathway_gsea function
#' @param plot_type A character string specifying the visualization type: "enrichment_plot", "dotplot", "barplot", "network", or "heatmap"
#' @param n_pathways An integer specifying the number of pathways to display
#' @param sort_by A character string specifying the sorting criterion: "NES", "pvalue", or "p.adjust"
#' @param colors A vector of colors for the visualization
#' @param abundance A data frame containing the original abundance data (required for heatmap visualization)
#' @param metadata A data frame containing sample metadata (required for heatmap visualization)
#' @param group A character string specifying the column name in metadata that contains the grouping variable (required for heatmap visualization)
#' @param network_params A list of parameters for network visualization
#' @param heatmap_params A list of parameters for heatmap visualization
#'
#' @return A ggplot2 object or ComplexHeatmap object
#' @export
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(ko_abundance)
#' data(metadata)
#'
#' # Prepare abundance data
#' abundance_data <- as.data.frame(ko_abundance)
#' rownames(abundance_data) <- abundance_data[, "#NAME"]
#' abundance_data <- abundance_data[, -1]
#'
#' # Run GSEA analysis
#' gsea_results <- pathway_gsea(
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment",
#'   pathway_type = "KEGG",
#'   method = "fgsea"
#' )
#'
#' # Create enrichment plot
#' visualize_gsea(gsea_results, plot_type = "enrichment_plot", n_pathways = 10)
#'
#' # Create dotplot
#' visualize_gsea(gsea_results, plot_type = "dotplot", n_pathways = 20)
#'
#' # Create barplot
#' visualize_gsea(gsea_results, plot_type = "barplot", n_pathways = 15)
#'
#' # Create network plot
#' visualize_gsea(gsea_results, plot_type = "network", n_pathways = 15)
#'
#' # Create heatmap
#' visualize_gsea(
#'   gsea_results,
#'   plot_type = "heatmap",
#'   n_pathways = 15,
#'   abundance = abundance_data,
#'   metadata = metadata,
#'   group = "Environment"
#' )
#' }
visualize_gsea <- function(gsea_results,
                          plot_type = "enrichment_plot",
                          n_pathways = 20,
                          sort_by = "p.adjust",
                          colors = NULL,
                          abundance = NULL,
                          metadata = NULL,
                          group = NULL,
                          network_params = list(),
                          heatmap_params = list()) {

  # Input validation
  if (!is.data.frame(gsea_results)) {
    stop("'gsea_results' must be a data frame")
  }

  if (!plot_type %in% c("enrichment_plot", "dotplot", "barplot", "network", "heatmap")) {
    stop("plot_type must be one of 'enrichment_plot', 'dotplot', 'barplot', 'network', or 'heatmap'")
  }

  if (!sort_by %in% c("NES", "pvalue", "p.adjust")) {
    stop("sort_by must be one of 'NES', 'pvalue', or 'p.adjust'")
  }

  if (!is.null(colors) && !is.character(colors)) {
    stop("colors must be NULL or a character vector")
  }

  # Check if required packages are installed
  if (plot_type == "enrichment_plot" || plot_type == "dotplot" || plot_type == "barplot") {
    if (!requireNamespace("enrichplot", quietly = TRUE)) {
      stop("Package 'enrichplot' is required. Please install it using BiocManager::install('enrichplot').")
    }
  }

  if (plot_type == "network") {
    if (!requireNamespace("igraph", quietly = TRUE) || !requireNamespace("ggraph", quietly = TRUE)) {
      stop("Packages 'igraph' and 'ggraph' are required for network plots. Please install them.")
    }
  }

  if (plot_type == "heatmap") {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || !requireNamespace("circlize", quietly = TRUE)) {
      stop("Packages 'ComplexHeatmap' and 'circlize' are required for heatmap plots. Please install them.")
    }

    # Check if required parameters are provided
    if (is.null(abundance) || is.null(metadata) || is.null(group)) {
      stop("For heatmap visualization, 'abundance', 'metadata', and 'group' parameters are required")
    }
  }

  # Set default colors if not provided
  if (is.null(colors)) {
    colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
  }

  # Sort results based on the specified criterion
  if (sort_by == "NES") {
    gsea_results <- gsea_results[order(abs(gsea_results$NES), decreasing = TRUE), ]
  } else if (sort_by == "pvalue") {
    gsea_results <- gsea_results[order(gsea_results$pvalue), ]
  } else if (sort_by == "p.adjust") {
    gsea_results <- gsea_results[order(gsea_results$p.adjust), ]
  }

  # Limit to top n_pathways
  if (nrow(gsea_results) > n_pathways) {
    gsea_results <- gsea_results[1:n_pathways, ]
  }

  # Create visualization based on plot_type
  if (plot_type == "enrichment_plot") {
    # Create enrichment plot
    # For this, we need to convert our results to a format compatible with enrichplot

    # Check if we have the necessary data
    if (!all(c("pathway_id", "NES", "pvalue", "p.adjust") %in% colnames(gsea_results))) {
      stop("GSEA results missing required columns for enrichment plot")
    }

    # Create a simple enrichment plot using ggplot2
    # In a real implementation, we would use enrichplot::gseaplot2
    # But for simplicity, we'll create a basic version

    # Sort by NES
    gsea_results <- gsea_results[order(gsea_results$NES), ]

    # Create a basic barplot of NES values
    p <- ggplot2::ggplot(gsea_results, ggplot2::aes(x = reorder(.data$pathway_name, .data$NES), y = .data$NES, fill = .data$p.adjust)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient(low = "red", high = "blue") +
      ggplot2::labs(
        title = "GSEA Enrichment Results",
        x = "Pathway",
        y = "Normalized Enrichment Score (NES)",
        fill = "Adjusted p-value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )

  } else if (plot_type == "dotplot") {
    # Create dotplot
    # Sort by NES
    gsea_results$pathway_name <- factor(gsea_results$pathway_name,
                                      levels = gsea_results$pathway_name[order(gsea_results$NES)])

    p <- ggplot2::ggplot(gsea_results,
                       ggplot2::aes(x = .data$NES, y = .data$pathway_name, color = .data$p.adjust, size = .data$size)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_gradient(low = "red", high = "blue") +
      ggplot2::labs(
        title = "GSEA Results",
        x = "Normalized Enrichment Score (NES)",
        y = "Pathway",
        color = "Adjusted p-value",
        size = "Gene Set Size"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )

  } else if (plot_type == "barplot") {
    # Create barplot
    # Sort by NES
    gsea_results$pathway_name <- factor(gsea_results$pathway_name,
                                      levels = gsea_results$pathway_name[order(gsea_results$NES)])

    # Add color based on NES direction
    gsea_results$direction <- ifelse(gsea_results$NES > 0, "Positive", "Negative")

    p <- ggplot2::ggplot(gsea_results,
                       ggplot2::aes(x = .data$pathway_name, y = .data$NES, fill = .data$direction)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8")) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "GSEA Results",
        x = "Pathway",
        y = "Normalized Enrichment Score (NES)",
        fill = "Direction"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )

  } else if (plot_type == "network") {
    # Set default network parameters
    default_params <- list(
      similarity_measure = "jaccard",
      similarity_cutoff = 0.3,
      layout = "fruchterman",
      node_color_by = "NES",
      edge_width_by = "similarity"
    )

    # Merge with user-provided parameters
    network_params <- utils::modifyList(default_params, network_params)

    # Create network plot
    p <- create_network_plot(
      gsea_results = gsea_results,
      n_pathways = n_pathways,
      similarity_measure = network_params$similarity_measure,
      similarity_cutoff = network_params$similarity_cutoff,
      layout = network_params$layout,
      node_color_by = network_params$node_color_by,
      edge_width_by = network_params$edge_width_by
    )

  } else if (plot_type == "heatmap") {
    # Set default heatmap parameters
    default_params <- list(
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_rownames = TRUE,
      annotation_colors = list(Group = stats::setNames(colors[seq_along(unique(metadata[[group]]))], unique(metadata[[group]])))
    )

    # Merge with user-provided parameters
    heatmap_params <- utils::modifyList(default_params, heatmap_params)

    # Create heatmap
    p <- create_heatmap_plot(
      gsea_results = gsea_results,
      abundance = abundance,
      metadata = metadata,
      group = group,
      n_pathways = n_pathways,
      cluster_rows = heatmap_params$cluster_rows,
      cluster_columns = heatmap_params$cluster_columns,
      show_rownames = heatmap_params$show_rownames,
      annotation_colors = heatmap_params$annotation_colors
    )
  }

  return(p)
}

#' Create network visualization of GSEA results
#'
#' @param gsea_results A data frame containing GSEA results from the pathway_gsea function
#' @param similarity_measure A character string specifying the similarity measure: "jaccard", "overlap", or "correlation"
#' @param similarity_cutoff A numeric value specifying the similarity threshold for filtering connections
#' @param n_pathways An integer specifying the number of pathways to display
#' @param layout A character string specifying the network layout algorithm: "fruchterman", "kamada", or "circle"
#' @param node_color_by A character string specifying the node color mapping: "NES", "pvalue", or "p.adjust"
#' @param edge_width_by A character string specifying the edge width mapping: "similarity" or "constant"
#'
#' @return A ggplot2 object
#' @keywords internal
create_network_plot <- function(gsea_results,
                               similarity_measure = "jaccard",
                               similarity_cutoff = 0.3,
                               n_pathways = 20,
                               layout = "fruchterman",
                               node_color_by = "NES",
                               edge_width_by = "similarity") {

  # Check required packages
  if (!requireNamespace("igraph", quietly = TRUE) ||
      !requireNamespace("ggraph", quietly = TRUE)) {
    stop("Packages 'igraph' and 'ggraph' are required for network plots. Please install them.")
  }

  # Limit number of pathways
  if (nrow(gsea_results) > n_pathways) {
    gsea_results <- gsea_results[order(gsea_results$p.adjust), ][1:n_pathways, ]
  }

  # Extract leading edge genes
  leading_edges <- strsplit(gsea_results$leading_edge, ";")
  names(leading_edges) <- gsea_results$pathway_id

  # Calculate pathway similarity
  n <- length(leading_edges)
  pathway_ids <- names(leading_edges)
  similarity_matrix <- matrix(0, nrow = n, ncol = n)
  rownames(similarity_matrix) <- pathway_ids
  colnames(similarity_matrix) <- pathway_ids

  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        set1 <- leading_edges[[i]]
        set2 <- leading_edges[[j]]

        # Handle empty sets
        if (length(set1) == 0 || length(set2) == 0) {
          similarity_matrix[i, j] <- 0
          next
        }

        if (similarity_measure == "jaccard") {
          # Jaccard similarity: |A∩B|/|A∪B|
          similarity_matrix[i, j] <- length(intersect(set1, set2)) / length(union(set1, set2))
        } else if (similarity_measure == "overlap") {
          # Overlap coefficient: |A∩B|/min(|A|,|B|)
          similarity_matrix[i, j] <- length(intersect(set1, set2)) / min(length(set1), length(set2))
        } else if (similarity_measure == "correlation") {
          # Simplified correlation measure
          similarity_matrix[i, j] <- length(intersect(set1, set2)) / sqrt(length(set1) * length(set2))
        }
      }
    }
  }

  # Apply similarity cutoff
  similarity_matrix[similarity_matrix < similarity_cutoff] <- 0

  # Check if there are any connections after applying cutoff
  if (sum(similarity_matrix) == 0) {
    # Create a simple plot with a message
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0,
                      label = "No significant connections found with current similarity cutoff") +
      ggplot2::theme_void()
    return(p)
  }

  # Create graph object
  graph <- igraph::graph_from_adjacency_matrix(
    similarity_matrix,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  # Check if graph is empty
  if (igraph::vcount(graph) == 0) {
    # Create a simple plot with a message
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0,
                      label = "No nodes found in the network") +
      ggplot2::theme_void()
    return(p)
  }

  # Add node attributes
  vertex_attr <- data.frame(
    name = pathway_ids,
    NES = gsea_results$NES[match(pathway_ids, gsea_results$pathway_id)],
    pvalue = gsea_results$pvalue[match(pathway_ids, gsea_results$pathway_id)],
    p.adjust = gsea_results$p.adjust[match(pathway_ids, gsea_results$pathway_id)],
    size = gsea_results$size[match(pathway_ids, gsea_results$pathway_id)],
    pathway_name = gsea_results$pathway_name[match(pathway_ids, gsea_results$pathway_id)],
    stringsAsFactors = FALSE
  )

  # Create a tidygraph object
  tbl_graph <- tidygraph::as_tbl_graph(graph) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(
      name = vertex_attr$name,
      NES = vertex_attr$NES,
      pvalue = vertex_attr$pvalue,
      p.adjust = vertex_attr$p.adjust,
      size = vertex_attr$size,
      pathway_name = vertex_attr$pathway_name
    )

  # Select layout algorithm
  if (layout == "fruchterman") {
    layout_name <- "fr"
  } else if (layout == "kamada") {
    layout_name <- "kk"
  } else if (layout == "circle") {
    layout_name <- "circle"
  } else {
    layout_name <- "fr"
  }

  # Create ggraph visualization
  p <- ggraph::ggraph(tbl_graph, layout = layout_name) +
    ggraph::geom_edge_link(ggplot2::aes(width = weight, alpha = weight)) +
    ggraph::geom_node_point(ggplot2::aes(color = .data[[node_color_by]], size = size)) +
    ggraph::geom_node_text(ggplot2::aes(label = pathway_name), repel = TRUE, size = 3) +
    ggplot2::scale_edge_width(range = c(0.1, 2)) +
    ggplot2::scale_edge_alpha(range = c(0.1, 0.8)) +
    ggplot2::scale_color_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      name = node_color_by
    ) +
    ggplot2::scale_size(range = c(2, 8), name = "Gene Set Size") +
    ggraph::theme_graph() +
    ggplot2::labs(
      title = "GSEA Pathway Network",
      subtitle = paste("Similarity measure:", similarity_measure, "| Cutoff:", similarity_cutoff)
    )

  return(p)
}

#' Create heatmap visualization of GSEA results
#'
#' @param gsea_results A data frame containing GSEA results from the pathway_gsea function
#' @param abundance A data frame containing the original abundance data
#' @param metadata A data frame containing sample metadata
#' @param group A character string specifying the column name in metadata that contains the grouping variable
#' @param n_pathways An integer specifying the number of pathways to display
#' @param cluster_rows A logical value indicating whether to cluster rows
#' @param cluster_columns A logical value indicating whether to cluster columns
#' @param show_rownames A logical value indicating whether to show row names
#' @param annotation_colors A list of colors for annotations
#'
#' @return A ComplexHeatmap object
#' @keywords internal
create_heatmap_plot <- function(gsea_results,
                               abundance,
                               metadata,
                               group,
                               n_pathways = 20,
                               cluster_rows = TRUE,
                               cluster_columns = TRUE,
                               show_rownames = TRUE,
                               annotation_colors = NULL) {

  # Check required packages
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    stop("Packages 'ComplexHeatmap' and 'circlize' are required for heatmap plots. Please install them.")
  }

  # Limit number of pathways
  if (nrow(gsea_results) > n_pathways) {
    gsea_results <- gsea_results[order(gsea_results$p.adjust), ][1:n_pathways, ]
  }

  # Extract leading edge genes
  leading_edges <- lapply(strsplit(gsea_results$leading_edge, ";"), function(x) x[x != ""])
  names(leading_edges) <- gsea_results$pathway_id

  # Create heatmap data matrix
  # For each pathway, calculate the average expression of leading edge genes
  heatmap_data <- matrix(0, nrow = length(leading_edges), ncol = ncol(abundance))
  rownames(heatmap_data) <- names(leading_edges)
  colnames(heatmap_data) <- colnames(abundance)

  for (i in seq_along(leading_edges)) {
    # Get pathway ID and genes
    current_pathway_id <- names(leading_edges)[i]
    genes <- leading_edges[[i]]

    # Ensure all genes are in abundance data
    genes <- genes[genes %in% rownames(abundance)]

    if (length(genes) > 0) {
      # Calculate average abundance
      heatmap_data[i, ] <- colMeans(abundance[genes, , drop = FALSE])
    }
  }

  # Scale data
  heatmap_data_scaled <- t(scale(t(heatmap_data)))

  # Prepare column annotation
  column_annotation <- metadata[[group]]
  names(column_annotation) <- rownames(metadata)

  # Ensure column annotation matches heatmap columns
  column_annotation <- column_annotation[colnames(heatmap_data)]

  # Create column annotation object
  ha <- ComplexHeatmap::HeatmapAnnotation(
    Group = column_annotation,
    col = list(Group = annotation_colors$Group),
    show_legend = TRUE
  )

  # Create row annotation (pathway enrichment scores)
  row_annotation <- data.frame(
    NES = gsea_results$NES,
    row.names = gsea_results$pathway_id
  )

  # Ensure row annotation matches heatmap rows
  row_annotation <- row_annotation[rownames(heatmap_data), , drop = FALSE]

  # Create row annotation object
  ra <- ComplexHeatmap::rowAnnotation(
    NES = row_annotation$NES,
    col = list(NES = circlize::colorRamp2(
      c(min(row_annotation$NES), 0, max(row_annotation$NES)),
      c("blue", "white", "red")
    )),
    show_legend = TRUE
  )

  # Create heatmap
  heatmap <- ComplexHeatmap::Heatmap(
    heatmap_data_scaled,
    name = "Z-score",
    col = circlize::colorRamp2(
      c(-2, 0, 2),
      c("blue", "white", "red")
    ),
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_names = show_rownames,
    row_names_gp = grid::gpar(fontsize = 8),
    top_annotation = ha,
    right_annotation = ra,
    row_title = "Pathways",
    column_title = "Samples",
    row_names_max_width = grid::unit(15, "cm")
  )

  return(heatmap)
}
