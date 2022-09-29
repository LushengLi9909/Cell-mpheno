library(ggplot2)
signif_results_h <- read.delim("data/signif_results_h.txt",header = TRUE,sep=",")

signif_results <- read.delim("data/signif_results_Mondo.txt",header = TRUE,sep=",")
signif_results_m = signif_results[(signif_results$MONDO_ID %in% signif_results_h$MONDO),]

mapping = data.frame()
for( d in unique(signif_results_m$MONDO_ID)){
  cell_m = unique(signif_results_m[signif_results_m$MONDO_ID ==d,"CellType"])
  cell_h = unique(signif_results_h[signif_results_h$MONDO ==d,"CellType"])
  same_cell = which(cell_m %in% cell_h)
  mapping = rbind(mapping, data.frame("MONDO_ID"=d,"s_cell"=length(same_cell),"n_cell_m"=length(cell_m),
                                      "n_cell_h"=length(cell_h),"p_cell_m" = sprintf("%.4f", (length(same_cell)/length(cell_m))),
                                      "p_cell_h" = sprintf("%.4f", (length(same_cell)/length(cell_h)))))
}
mapping$p_cell_m <- as.numeric(mapping$p_cell_m)
mapping$p_cell_h <- as.numeric(mapping$p_cell_h)
mapping_plot = data.frame()
mapping_plot = rbind(mapping_plot, data.frame("PO"=c(rep("MPO",length(mapping$p_cell_m)),rep("HPO",length(mapping$p_cell_h))),
                                    "percent" = c(mapping$p_cell_m,mapping$p_cell_h)))
mapping_plot$PO = factor(mapping_plot$PO,levels = c("MPO","HPO"))

mapping_plt = ggplot(mapping_plot, aes(x = PO, y = percent,fill= PO))+
  geom_boxplot(width=0.6,outlier.size=0.2)+
  stat_boxplot(geom = "errorbar",
               lwd=0.5,
               width=0.2)+
  scale_y_continuous(labels = scales::percent)+
  geom_jitter(color="black", size=0.05, alpha=0.9,width = 0.3)+
  cowplot::theme_cowplot()+
  theme(axis.title.x = element_blank(),legend.position ="none")+
  ylab("Percent of the same significant cell type")

