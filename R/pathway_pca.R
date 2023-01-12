#' Title
#'
#' @param abundance
#' @param metadata
#' @param group
#'
#' @return
#' @export
#'
#' @examples
pathway_pca <- function(abundance, metadata, group){
  pca_axis <- prcomp(t(abundance), center = TRUE, scale = TRUE)$x[,1:2]
  pca_proportion <- prcomp(t(abundance), center = TRUE, scale = TRUE)$sdev[1:2]/sum(prcomp(t(abundance), center = TRUE, scale = TRUE)$sdev)*100
  pca <- cbind(pca_axis, metadata[,group])
  pca$Group <- pca[,group]
  Fig1a.taxa.pca <- ggplot(pca,aes(PC1,PC2))+
    geom_point(size=2,aes(color=Group),show.legend = T)+
    scale_color_manual(values=c("#d93c3e", "#3685bc"))+
    stat_ellipse(aes(color = Group),fill="white",geom = "polygon",
                 level=0.95,alpha = 0.01,show.legend = F)+
    labs(x=paste0("PC1(",round(pca_proportion[1],1),"%)"),y=paste0("PC1(",round(pca_proportion[2],1),"%)"),color = group)+
    theme_classic()+
    theme(axis.line=element_line(colour = "black"),
          axis.title=element_text(color="black",face = "bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text = element_text(color="black",size=10,face = "bold"),
          legend.text = element_text(size = 8, face = "bold"),
          legend.title = element_text(size = 10, face = "bold"))

  Fig1a.taxa.pc1.density <-
    ggplot(pca) +
    geom_density(aes(x=PC1, group=Group, fill=Group),
                 color="black", alpha=1,position = 'identity',
                 show.legend = F) +
    scale_fill_manual(values=c("#d93c3e", "#3685bc")) +
    scale_y_discrete(expand = c(0,0.001))+
    labs(x=NULL,y=NULL)+
    theme_classic() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank())

  Fig1a.taxa.pc2.density <-
    ggplot(pca) +
    geom_density(aes(x=PC2, group=Group, fill=Group),
                 color="black", alpha=1, position = 'identity',show.legend = F) +
    scale_fill_manual(values=c("#d93c3e", "#3685bc")) +
    theme_classic()+
    scale_y_discrete(expand = c(0,0.001))+
    labs(x=NULL,y=NULL)+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    coord_flip()

  Fig1a.taxa.pca %>%
    insert_top(Fig1a.taxa.pc1.density,height = 0.3) %>%
    insert_right(Fig1a.taxa.pc2.density,width=0.3) %>%
    as.ggplot()
}
