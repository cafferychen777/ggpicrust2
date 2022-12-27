ggpicrust2 <-
  function(file,
           metadata,
           group,
           pathway,
           daa_method = "ALDEx2",
           ko_to_kegg = FALSE,
           p.adjust = "BH",
           order = "group",
           p_values_bar = TRUE,
           x_lab = "pathway_name",
           select = NULL,
           reference = NULL,
           p_values_bar = TRUE,
           colors = NULL)
  {
    switch(ko_to_kegg,
           "TRUE" = {
             if (ko_to_kegg = TRUE) {
               abundance <- ko2kegg_abundance(file)
             }
             daa_results_df <-
               pathway_daa(
                 abundance = abundance,
                 metadata = metadata,
                 group = group,
                 daa_method = daa_method,
                 select = select,
                 p.adjust = p.adjust,
                 reference = reference,
               )
    if (x_lab == "pathway_name") {
      daa_results_df  <-
        pathway_annotation(daa_results_df = daa_results_df)
    }
    j <- 1
    for (i in unique(daa_results_df$method)) {
      daa_sub_method_results_df <- daa_results_df[daa_results_df[,"method"] == i,]
      combination_bar_plot <-
        pathway_errorbar(
          abundance,
          daa_sub_method_results_df,
          metadata[, group],
          ko_to_kegg = ko_to_kegg,
          order = "pathway_class",
          colors = colors,
          x_lab = x_lab
        )
      print(combination_bar_plot)
      message(paste0("No.", j, " plot is method ", i))
    }
           },
    "FALSE" = {
      abundance <-
        read_delim(
          file,
          delim = "\t",
          escape_double = FALSE,
          trim_ws = TRUE
        )
      abundance <- column_to_rownames(abundance, var = "function")
      daa_results_df <-
        pathway_daa(abundance = abundance,
                    metadata = metadata,
                    group = group,
                    daa_method = daa_method,
                    select = select,
                    p.adjust = p.adjust,
                    reference = reference)
      daa_results_df <- pathway_annotation(pathway = pathway, ko_to_kegg = FALSE, daa_results_df = daa_results_df)
      for (i in unique(daa_results_df$method)) {
        daa_sub_method_results_df <- daa_results_df[daa_results_df[,"method"] == i,]
        combination_bar_plot <-
          pathway_errorbar(
            abundance,
            daa_sub_method_results_df,
            metadata[, group],
            ko_to_kegg = ko_to_kegg,
            order = order,
            colors = colors,
            x_lab = x_lab
          )
        print(combination_bar_plot)
        message(paste0("No.", j, " plot is method ", i))
      }
      pathway_errorbar(abundance = abundance, metadata = metadata, Group = metadata[,group], order = order, )
    })
  }
