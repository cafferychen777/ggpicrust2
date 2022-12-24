pathway_errorbar <-
  function(abundance,
           daa_results_df,
           Group,
           p_values_threshold = 0.05,
           order = ,
           select = NULL,
           boxplot = TRUE) {
    errorbar_abundance_mat <- as.matrix(abundance)
    daa_resuls_filtered_df <-
      daa_results_df[daa_results_df$p_adjust < p_values_threshold, ]
    errorbar_relative_abundance_mat <-
      make_relative(errorbar_abundance_mat)
    errorbar_relative_sub_abundance_mat <-
      errorbar_relative_abundance_mat[rownames(errorbar_relative_abundance_mat) %in% daa_resuls_filtered_df$feature, ]
    #这里可能要加一个排序
    error_bar_matrix <-
      cbind(
        sample = colnames(errorbar_relative_sub_abundance_mat),
        group = Group,
        t(errorbar_relative_sub_abundance_mat)
      )
    error_bar_df <- as.data.frame(error_bar_matrix)
    error_bar_pivot_longer_df <-
      pivot_longer(error_bar_df, -c(sample, group))
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

    bar_errorbar <-
      ggplot(error_bar_pivot_longer_tibble_summarised,
             aes(name, mean, fill = group)) +
      geom_errorbar(
        aes(ymax = mean + sd, ymin = 0),
        position = position_dodge(width = 0.8),
        width = 0.5,
        size = 0.5,
        color = "black"
      ) +
      geom_bar(stat = "identity",
               position = position_dodge(width = 0.8),
               width = 0.8) +
      geom_stripped_cols() +
      scale_fill_manual(values = c("#fdb462", "#7fb1d3")) +
      scale_color_manual(values = c("#fdb462", "#7fb1d3")) +
      theme_prism() +
      scale_y_continuous(expand = c(0, 0),
                         #limits = c(0, 0.1),
                         guide = "prism_offset_minor") +
      labs(y = "Relative Abundance(%)", x = NULL) +
      theme(
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.5),
        axis.ticks.x = element_line(size = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = 10, color = "black"),
        axis.text.y = element_text(margin = margin(r = 0)),
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
    ggsave("bar_errorbar.pdf",
           bar_errorbar,
           width = 8,
           height = 6)



    daa_resuls_filtered_df <-
      cbind(
        daa_resuls_filtered_df,
        negative_log10_p = -log10(daa_resuls_filtered_df$p_adjust),
        group_nonsense = "nonsense"
      )


    p_values_bar <- daa_resuls_filtered_df %>%
      ggplot(aes(feature,negative_log10_p,fill=group_nonsense))+
      geom_bar(stat="identity",position=position_dodge(width=0.8),width = 0.8)+
      labs(y="P Value (-log 10)",x=NULL)+
      geom_stripped_cols() +
      scale_fill_manual(values = c("#7fb1d3")) +
      scale_color_manual(values = c("#7fb1d3")) +
      geom_hline(aes(yintercept = 1.30103),linetype='dashed',color = 'black') +
      theme_prism() +
      scale_y_continuous(expand = c(0, 0),
                         #limits = c(0, 0.1),
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
    ggsave("p_values_bar.pdf", p_values_bar, width = 8, height = 6)


    daa_resuls_filtered_df$p_adjust <- as.character(daa_resuls_filtered_df$p_adjust)
    daa_resuls_filtered_df$unique <- nrow(daa_resuls_filtered_df)-order(daa_resuls_filtered_df$feature)+1
    daa_resuls_filtered_df$p_adjust <- substr(daa_resuls_filtered_df$p_adjust,1,5)
    p_annotation <- daa_resuls_filtered_df %>%
      ggplot(aes(group_nonsense,p_adjust))+
      geom_text(aes(group2,unique,label=p_adjust),size=3.5,color="black") +
      labs(y="P_value")+
      scale_y_discrete(position = "right")+
      #scale_y_continuous(breaks = daa_resuls_filtered_df$p_adjust)+
      theme_prism()+
      theme(axis.ticks=element_blank(),
            axis.line = element_blank(),
            panel.grid.major.y =element_blank(),
            panel.grid.major.x = element_blank(),
            panel.background = element_blank(),
            axis.text=element_blank(),
            plot.margin = unit(c(0,0.2,0,0),"cm"),
            axis.title.y =  element_text(size=11,color = "black",hjust = 0.5),
            axis.title.x=element_blank(),
            legend.position = "non")




    bar_errorbar + p_values_bar + p_annotation




    lefserPlot(Lefser_results,
               colors = c("#7fb1d3", "#fdb462"),
               trim.names = FALSE)
  }
