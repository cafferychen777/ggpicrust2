#' Title
#'
#' @param abundance A data frame, predicted functional pathway abundance
#' @param metadata A tibble, consisting of samples information
#' @param group A character, group name
#'
#' @return
#' @export
#'
#' @examples
pathway_heatmap <- function(abundance, metadata, group){
  abundance <- make_relative(as.matrix(abundance))
  abundance <- as.data.frame(abundance)
  order <- metadata[order(metadata$Enviroment),]$sample_name
  new_abundance <- abundance %>% rownames_to_column()%>% pivot_longer(-rowname)
  new_abundance$name <- factor(new_abundance$name, levels = order)
  new_abundance %>% ggplot(aes(name, rowname,fill=value))+
    geom_tile()+
    labs(x=NULL,y=NULL)+
    scale_y_discrete(expand=c(0,0),position="left")+
    scale_x_discrete(expand=c(0,0))+
    scale_fill_gradientn(colours = c("#273b68", "#2190bc","#32b7d2","#62c3c3","#95d1b6","#cbdea7","#fbefa6"))+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(color="black",size=10,face = "bold"),
          legend.title = element_blank(),
          legend.text = element_blank())+
    guides(fill=guide_colorbar(direction = "vertical",reverse = F,barwidth = unit(.6, "cm"),
                               barheight = unit(18.5,"cm")))
}
