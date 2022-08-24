library(ontologyIndex)
library(ggplot2)
library(cowplot)
all_results_merged <- readRDS("data/mpo_results_merged.rds")
mpo <- get_ontology("data/MPheno_OBO.ontology")
ctd <- readRDS("data/TabulaMuris_n.rds")

branch = "abnormal consumption behavior"
branch_id = paste(mpo$id[match(branch,mpo$name)])
branch_descendants = paste(get_descendants(mpo,branch_id))
branch_descendants_names = paste(mpo$name[branch_descendants])

all_results_merged$cur_branch = paste("Other")
all_results_merged$cur_branch[all_results_merged$Phenotype %in% branch_descendants_names] = branch
branch_signif_counts = data.frame()
for (c in unique(all_results_merged$CellType)) {
  n_signif = length(all_results_merged[all_results_merged$CellType == c & all_results_merged$cur_branch == branch & all_results_merged$q < 0.05, ]$q)
  branch_signif_counts = rbind(branch_signif_counts,
                               data.frame("branch"=branch,
                                          "CellType"=c,
                                          "n_signif"=n_signif))
}

cell_order = factor(ctd[[1]]$plotting$cell_ordering,levels = ctd[[1]]$plotting$cell_ordering)
branch_signif_counts$CellType = factor(branch_signif_counts$CellType,levels = cell_order)

branch_signif_counts$labels = branch_signif_counts$n_signif
branch_signif_counts$labels[branch_signif_counts$labels == 0] = ""
branch_plt <- ggplot(branch_signif_counts[branch_signif_counts$branch==branch,], aes(x=CellType,y=n_signif)) +
  geom_col( fill = "#69b3a2", color = "#e9ecef") +
  geom_text(mapping= aes(label = labels, y = n_signif + 2))+
  theme_cowplot()+
  ylab("N phenotypes") +
  scale_y_continuous(expand=c(0,0), limits= c(0,max(branch_signif_counts$n_signif)+3))+
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),legend.position="none") +
  #coord_flip() +
  ggtitle("Significant enrichments per cell: abnormal consumption behavior")

branch_descendants = all_results_merged[(all_results_merged$CellType == "Satellite cells" | all_results_merged$CellType == "Pericytes" | all_results_merged$CellType == "Epithelial cells") & all_results_merged$cur_branch == branch & all_results_merged$q < 0.05, "Phenotype"]
branch_descendants_names = unique(branch_descendants)
branch_descendants_names = branch_descendants_names[-which(branch_descendants_names ==branch)]
branch_descendants = all_results_merged[(all_results_merged$CellType == "Satellite cells") & all_results_merged$cur_branch == branch & all_results_merged$q < 0.05, "Phenotype"]

pheno_df = all_results_merged[all_results_merged$Phenotype %in% branch_descendants_names, ]
pheno_df$signif_asterics = ""
pheno_df$signif_asterics[pheno_df$q<0.05] = "*"
pheno_df$signif_asterics[pheno_df$q<0.005] = "**"
pheno_df$signif_asterics[pheno_df$q<0.0005] = "***"
pheno_df$signif_asterics[pheno_df$q<0.00005] = "****"
pheno_df$CellType = factor(pheno_df$CellType,levels = cell_order)
consumpt_fold_plt <- ggplot(pheno_df, aes(x = CellType, y= fold_change, fill = Phenotype)) +
  ggtitle(paste0('Descendants of MPO branch "',branch,'"')) +
  geom_col() +
  geom_text(label = pheno_df$signif_asterics, mapping = aes(y = fold_change + 0.1))+
  theme_cowplot() +
#  scale_y_continuous(expand=c(0,0), limits= c(0,max(branch_signif_counts$n_signif)+2))+
  ylab("Fold change") +
#  xlab("Cell type") +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust = 0.3),axis.title.x = element_blank(),legend.position = "none") +
  facet_wrap(~Phenotype, ncol=1)

all_branch_plt = cowplot::plot_grid(branch_plt,consumpt_fold_plt,axis = "lr",
                                align = "v", ncol = 1,
                                rel_heights = c(0.12,1),labels=c("A","B"))
