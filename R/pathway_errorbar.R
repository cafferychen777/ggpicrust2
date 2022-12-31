#' visualization of predictional functional pathway DA results
#'
#' @param abundance A data frame, pathway is its row name, sample is its row name
#' @param daa_results_df A data frame, output of pathway_daa
#' @param Group A dataframe or a vector, such as metadata$group
#' @param ko_to_kegg A charachter to control if converting ko abundance to kegg abundance
#' @param p_values_threshold A numerir
#' @param order A character to control the main plot rows' order
#' @param select A vector consisting of pathway names
#' @param p_value_bar A character to control if the main plot has the p_values_bar
#' @param colors A vector consisting of colors number
#' @param x_lab A character to control x_lab name
#'
#' @return A figure
#' @export A return
#'
#' @examples
pathway_errorbar <-
  function(abundance,
           daa_results_df,
           Group,
           ko_to_kegg = FALSE,
           p_values_threshold = 0.05,
           order = "group",
           select = NULL,
           p_value_bar = TRUE,
           colors = NULL,
           x_lab = NULL) {
    if (is.null(x_lab)){
      if (ko_to_kegg == TRUE){
        x_lab <- "pathway_name"
      }else{
        x_lab <- "description"
      }
    }
    if (is.null(colors)) {
      colors <- c("#d93c3e", "#3685bc", "#87ceeb")
    }
    errorbar_abundance_mat <- as.matrix(abundance)
    daa_results_filtered_df <-
      daa_results_df[daa_results_df$p_adjust < p_values_threshold,]
    if (!is.null(select)) {
      daa_results_filtered_sub_df <-
        daa_results_filtered_df[daa_results_filtered_df$feature %in% select, ]
    } else {
      daa_results_filtered_sub_df <- daa_results_filtered_df
    }
    if (nrow(daa_results_filtered_sub_df) > 30) {
      stop(
        paste0(
          "The feature with statistically significane are more than 30, the visualization will be terrible.\n Please use select to reduce the number.\n Now you have ",
          paste(daa_results_filtered_sub_df$feature, collapse = ",")
        )
      )
    }
    errorbar_sub_abundance_mat <-
      errorbar_abundance_mat[rownames(errorbar_abundance_mat) %in% daa_results_filtered_sub_df$feature,]
    errorbar_sub_relative_abundance_mat <-
      make_relative(errorbar_sub_abundance_mat)
    error_bar_matrix <-
      cbind(
        sample = colnames(errorbar_sub_relative_abundance_mat),
        group = Group,
        t(errorbar_sub_relative_abundance_mat)
      )
    error_bar_df <- as.data.frame(error_bar_matrix)
    error_bar_pivot_longer_df <-
      pivot_longer(error_bar_df,-c(sample, group))
    error_bar_pivot_longer_tibble <-
      mutate(error_bar_pivot_longer_df, group = as.factor(group))
    error_bar_pivot_longer_tibble$sample <-
      factor(error_bar_pivot_longer_tibble$sample)
    error_bar_pivot_longer_tibble$name <-
      factor(error_bar_pivot_longer_tibble$name)
    error_bar_pivot_longer_tibble$value <-
      as.numeric(error_bar_pivot_longer_tibble$value)
    error_bar_pivot_longer_tibble_summarised <-
      error_bar_pivot_longer_tibble %>% group_by(name, group) %>%
      summarise(mean = mean(value), sd = sd(value))
    error_bar_pivot_longer_tibble_summarised <-
      error_bar_pivot_longer_tibble_summarised %>% mutate(group2 = "nonsense")
    switch(
      order,
      "p_values" = {
        order <- order(daa_results_filtered_sub_df$p_adjust)
      },
      "name" = {
        order <- order(daa_results_filtered_sub_df$feature)
      },
      "group" = {
        daa_results_filtered_sub_df$pro <- 1
        for (i in levels(error_bar_pivot_longer_tibble_summarised$name)) {
          error_bar_pivot_longer_tibble_summarised_sub <-
            error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name ==
                                                       i,]
          pro_group <-
            error_bar_pivot_longer_tibble_summarised_sub[error_bar_pivot_longer_tibble_summarised_sub$mean ==
                                                           max(error_bar_pivot_longer_tibble_summarised_sub$mean),]$group
          pro_group <- as.vector(pro_group)
          daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature ==
                                        i,]$pro <- pro_group
        }
        order <-
          order(daa_results_filtered_sub_df$pro,
                daa_results_filtered_sub_df$p_adjust)
      },
      "pathway_class" = {
        if (!"pathway_class" %in% colnames(daa_results_filtered_sub_df)) {
          stop("Please use pathway_annotation function to annotation the pathway_daa results")
        }
        order <- order(
          daa_results_filtered_sub_df$pathway_class,
          daa_results_filtered_sub_df$p_adjust
        )
      },
      {
        order <- order
      }
    )
    daa_results_filtered_sub_df <-
      daa_results_filtered_sub_df[order,]
    error_bar_pivot_longer_tibble_summarised_ordered <-
      data.frame(
        name = NULL,
        group = NULL,
        mean = NULL,
        sd = NULL
      )
    for (i in daa_results_filtered_sub_df$feature) {
      error_bar_pivot_longer_tibble_summarised_ordered <-
        rbind(
          error_bar_pivot_longer_tibble_summarised_ordered,
          error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name ==
                                                     i,]
        )
    }
    if (ko_to_kegg == FALSE){
      error_bar_pivot_longer_tibble_summarised_ordered[, x_lab] <-
        rep(daa_results_filtered_sub_df[, x_lab], each = length(levels(
          factor(error_bar_pivot_longer_tibble_summarised_ordered$group)
        )))
    }


    # levels(error_bar_pivot_longer_tibble_summarised_ordered$name) <-
    #   rev(daa_results_filtered_sub_df$feature)
    #
    if (ko_to_kegg == TRUE) {
      error_bar_pivot_longer_tibble_summarised_ordered$pathway_class <-
        rep(daa_results_filtered_sub_df$pathway_class,
            each = length(levels(
              factor(error_bar_pivot_longer_tibble_summarised_ordered$group)
            )))
    }
    error_bar_pivot_longer_tibble_summarised_ordered$name <- factor(error_bar_pivot_longer_tibble_summarised_ordered$name, levels = rev(daa_results_filtered_sub_df$feature))

    #error_bar_pivot_longer_tibble_summarised_ordered$order <- rep(0:(nrow(daa_results_filtered_sub_df)-1),each=2)
    bar_errorbar <-
      ggplot(error_bar_pivot_longer_tibble_summarised_ordered,
             aes(mean, name, fill = group)) +
      geom_errorbar(
        aes(xmax = mean + sd, xmin = 0),
        position = position_dodge(width = 0.8),
        width = 0.5,
        size = 0.5,
        color = "black"
      ) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 0.8),
               width = 0.8) +
      geom_stripped_cols() +
      scale_fill_manual(values = c(colors[1], colors[2])) +
      scale_color_manual(values = c(colors[1], colors[2])) +
      theme_prism() +
      scale_x_continuous(expand = c(0, 0),
                         guide = "prism_offset_minor",) +
      scale_y_discrete(labels = rev(daa_results_filtered_sub_df[, x_lab])) +
      labs(x = "Relative Abundance(%)", y = NULL) +
      theme(
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.5),
        axis.ticks.x = element_line(size = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(margin = margin(r = 0)),
        axis.text.y = element_text(
          size = 10,
          color = "black",
          margin = margin(b = 6)
        ),
        axis.title.x =  element_text(
          size = 10,
          color = "black",
          hjust = 0.5
        ),
        legend.position = "top",
        legend.key.size = unit(0.1, "cm"),
        legend.direction = "vertical",
        legend.justification = "left",
        legend.text = element_text(size = 8, face = "bold"),
        legend.box.just = "right",
        plot.margin = margin(0, 0.5, 0.5, 0, unit = "cm")
      ) + coord_cartesian(clip = "off")








    if (ko_to_kegg == TRUE) {
      bar_errorbar_data <- bar_errorbar$data
      bar_errorbar_aes_x <-
        ggiraphExtra::getMapping(object$mapping, "x")
      bar_errorbar_aes_y <-
        ggiraphExtra::getMapping(object$mapping, "y")
      bar_errorbar_data_x <- data[, c(bar_errorbar_aes_x)]
      bar_errorbar_data_y <- data[, c(bar_errorbar_aes_y)]
      pathway_class_group <-
        daa_results_filtered_sub_df$pathway_class %>%
        table() %>%
        data.frame()
      start <-
        c(1, rev(pathway_class_group$Freq)[1:(length(pathway_class_group$Freq) - 1)]) %>%
        cumsum()
      end <- cumsum(rev(pathway_class_group$Freq))
      ymin <- start - 1 / 2
      ymax <- end + 1 / 2
      nPoints <- length(start)
      pCol <- c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B",
                 "#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4",
                 "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
                 "#6E4B9E","#0C727C", "#7E1416","#D8A767","#3D3D3D")[1:nPoints]
      pFill <- pCol
      for (i in 1:nPoints)  {
        bar_errorbar <- bar_errorbar +
          ggplot2::annotation_custom(
            grob = grid::rectGrob(
              gp = grid::gpar(
                col = pCol[i],
                fill = pFill[i],
                lty = NULL,
                lwd = NULL,
                alpha = 0.2
              )
            ),
            xmin = ggplot2::unit(-2, 'native'),
            xmax = ggplot2::unit(0, 'native'),
            ymin = ggplot2::unit(ymin[i], 'native'),
            ymax = ggplot2::unit(ymax[i], 'native')
          )
      }
    }





    daa_results_filtered_sub_df <-
      cbind(
        daa_results_filtered_sub_df,
        negative_log10_p = -log10(daa_results_filtered_sub_df$p_adjust),
        group_nonsense = "nonsense"
      )


    p_values_bar <- daa_results_filtered_sub_df %>%
      ggplot(aes(feature, negative_log10_p, fill = group_nonsense)) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 0.8),
               width = 0.8) +
      labs(y = "P Value (-log 10)", x = NULL) +
      geom_stripped_cols() +
      scale_fill_manual(values = colors[3]) +
      scale_color_manual(values = colors[3]) +
      geom_hline(aes(yintercept = 1.30103),
                 linetype = 'dashed',
                 color = 'black') +
      theme_prism() +
      scale_y_continuous(expand = c(0, 0),
                         guide = "prism_offset_minor") +
      labs(y = "P Values (-log10)", x = NULL) +
      theme(
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.5),
        axis.ticks.x = element_line(size = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.text.y = element_blank(),
        axis.text.x = element_text(
          size = 10,
          color = "black",
          margin = margin(b = 6)
        ),
        axis.title.x =  element_text(
          size = 11,
          color = "black",
          hjust = 0.5
        ),
        legend.position = "non"
      ) +
      coord_flip()




    if (ko_to_kegg == TRUE) {
      pathway_class_y <- (ymax + ymin) / 2 - 0.5
      pathway_class_plot_df <-
        as.data.frame(
          cbind(
            nonsense = "nonsense",
            pathway_class_y = pathway_class_y,
            pathway_class = rev(unique(
              daa_results_filtered_sub_df$pathway_class
            ))
          )
        )
      pathway_class_plot_df$pathway_class_y <-
        as.numeric(pathway_class_plot_df$pathway_class_y)
      pathway_class_annotation <-
        pathway_class_plot_df %>% ggplot(aes(nonsense, pathway_class_y)) + geom_text(
          aes(nonsense, pathway_class_y, label = pathway_class),
          size = 3.5,
          color = "black",
          fontface = "bold",
          family = "sans"
        ) +
        scale_y_discrete(position = "right") +
        theme_prism() +
        theme(
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          plot.margin = unit(c(0, 0.2, 0, 0), "cm"),
          axis.title.y =  element_blank(),
          axis.title.x = element_blank(),
          legend.position = "non"
        )
    }








    daa_results_filtered_sub_df$p_adjust <-
      as.character(daa_results_filtered_sub_df$p_adjust)
    daa_results_filtered_sub_df$unique <-
      nrow(daa_results_filtered_sub_df) - seq_len(nrow(daa_results_filtered_sub_df)) + 1
    daa_results_filtered_sub_df$p_adjust <-
      substr(daa_results_filtered_sub_df$p_adjust, 1, 5)
    p_annotation <- daa_results_filtered_sub_df %>%
      ggplot(aes(group_nonsense, p_adjust)) +
      geom_text(
        aes(group_nonsense, unique, label = p_adjust),
        size = 3.5,
        color = "black",
        fontface = "bold",
        family = "sans"
      ) +
      labs(y = "p-value (corrected)") +
      scale_y_discrete(position = "right") +
      theme_prism() +
      theme(
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        plot.margin = unit(c(0, 0.2, 0, 0), "cm"),
        axis.title.y =  element_text(
          size = 11,
          color = "black",
          hjust = 0.5
        ),
        axis.title.x = element_blank(),
        legend.position = "non"
      )
    if (p_value_bar == TRUE) {
      if (ko_to_kegg == TRUE) {
        combination_bar_plot <-
          pathway_class_annotation + bar_errorbar + p_values_bar + p_annotation + plot_layout(ncol = 4, width =
                                                                                                c(1, 1.5, 0.5, 0.1))
      }
      else{
        combination_bar_plot <-
          bar_errorbar + p_values_bar + p_annotation + plot_layout(ncol = 3, width = c(2.3, 0.7, 0.3))
      }
    }else{
      combination_bar_plot <-
        bar_errorbar + p_annotation + plot_layout(ncol = 2, width = c(2.5,  0.2))
    }
    return(combination_bar_plot)
  }
