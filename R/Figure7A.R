library(ggplot2)
library(cowplot)
all_results_merged_mpo <- readRDS("data/mpo_results_merged.rds")
all_results_merged_hpo <- readRDS("data/results_2022-07-17 23-06-35.rds")
all_results_merged_hpo$Phenotype=gsub("-SLASH-","/",all_results_merged_hpo$Phenotype)
all_results_merged_hpo$CellType = gsub("_"," ",all_results_merged_hpo$CellType)

mapping_all <- read.delim("data/ph_mp_mapping.txt",header = TRUE,sep=",")
gene_data_hpo = read.delim("data/phenotype_to_genes.txt", header = FALSE, skip =1)
colnames(gene_data_hpo) = c("HPID", "Phenotype", "EntrezID", "Gene", "Additional", "Source", "LinkID")
gene_data_hpo = gene_data_hpo[,c("Phenotype","Gene")]
gene_data_hpo = gene_data_hpo[!duplicated(gene_data_hpo),]
gene_data_mpo <- read.delim("data/Phenotype_MPO.txt",header = TRUE,sep=",")
mapping_all$overlapping = NA
mapping_all$gene_mpo =NA
mapping_all$gene_hpo =NA
mapping_all$MPO =NA
mapping_all$HPO =NA
for(i in 1:length(mapping_all$HPterm)){
  Gene_hpo = gene_data_hpo[gene_data_hpo$Phenotype == mapping_all$HPterm[i] , "Gene"]
  Gene_mpo = gene_data_mpo[gene_data_mpo$Phenotype == mapping_all$MPterm[i] , "Gene"]
  mapping_all$overlapping[i] = sum((Gene_hpo %in% Gene_mpo) == TRUE)
  mapping_all$gene_mpo[i] = length(Gene_mpo)
  mapping_all$gene_hpo[i] = length(Gene_hpo)
  mapping_all$MPO[i] = sprintf("%.4f", (mapping_all$overlapping[i]/mapping_all$gene_mpo[i]))
  mapping_all$HPO[i] = sprintf("%.4f", (mapping_all$overlapping[i]/mapping_all$gene_hpo[i]))
}
mapping_all= mapping_all[which((mapping_all$HPterm %in% all_results_merged_hpo$Phenotype) == TRUE),]
mapping_all= mapping_all[which((mapping_all$MPterm %in% all_results_merged_mpo$Phenotype) == TRUE),]

mapping_all$MPO <- as.numeric(mapping_all$MPO)
mapping_all$HPO <- as.numeric(mapping_all$HPO)


mapping = data.frame()
mapping = rbind(mapping, data.frame("PO"=c(rep("MPO",length(mapping_all$MPO)),rep("HPO",length(mapping_all$MPO))),
                                           "percent" = c(mapping_all$MPO,mapping_all$HPO)))
mapping$PO = factor(mapping$PO,levels = c("MPO","HPO"))
mapping_gene_plot = ggplot(mapping, aes(x = PO, y = percent,fill= PO))+
  geom_boxplot(width=0.6,outlier.size=0.2)+
  stat_boxplot(geom = "errorbar",
               lwd=0.5,
               width=0.2)+
  scale_y_continuous(labels = scales::percent)+
  geom_jitter(color="black", size=0.05, alpha=0.9,width = 0.3)+
  cowplot::theme_cowplot()+
  theme(axis.title.x = element_blank(),legend.position ="none")+
  ylab("Percent of overlapping gene")

