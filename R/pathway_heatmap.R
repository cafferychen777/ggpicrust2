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
  order <- metadata[order(metadata$Enviroment),]$sample_name
  new_abundance <- abundance %>% rownames_to_column()%>% pivot_longer(-rowname)
  new_abundance$name <- factor(new_abundance$name, levels = order)
  new_abundance %>% ggplot(aes(name, rowname,fill=value))+
    geom_tile()+
    labs(x=NULL,y=NULL)+
    scale_y_discrete(expand=c(0,0),position="left")+
    scale_x_discrete(expand=c(0,0))+
    scale_fill_gradientn(colours = c("#273b68", "#2190bc","#32b7d2","#62c3c3","#95d1b6","#cbdea7","#fbefa6"))+
    theme(axis.text.x=element_text(color="black",angle=90,vjust=0.5,size=8,hjust=1),
          axis.text.y=element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(color="black",size=10,face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size=8,color="black"))+
    guides(fill=guide_colorbar(direction = "vertical",reverse = F,barwidth = unit(.6, "cm"),
                               barheight = unit(18.5,"cm")))
}
