library(ontologyIndex)
library(ggplot2)
###MPO
blank_x <-   theme(axis.title.x = element_blank())
gene_data <- read.delim("data/Phenotype_MPO.txt",header = TRUE,sep=",")
ontology <- get_ontology("data/MPheno_OBO.ontology")
all_results_merged <- readRDS("data/mpo_results_merged.rds")

ont_levels <- data.frame("phenotype"=unique(gene_data$Phenotype),
                         "ontology_id"=ontology$id[match(unique(gene_data$Phenotype),ontology$name)])

get_ont_level = function(ontology,term_ids) {
  children = unique(setdiff(unlist(ontology$children[term_ids]), term_ids))
  if (length(children) == 0) {
    return(0)
  } else {
    return(1 + get_ont_level(ontology,children)) #<- recursion..
  }
}
ont_levels <- ont_levels[complete.cases(ont_levels),]
lvls <- c()
for (id in ont_levels$ontology_id) {
  lvls = append(lvls, get_ont_level(ontology,id))
}
ont_levels$ont_lev <- lvls
rm(lvls)
n_associated_cells <- c()
for (p in ont_levels$phenotype) {
  n_associated_cells <- append(n_associated_cells, length(all_results_merged$CellType[all_results_merged$Phenotype == p & all_results_merged$q < 0.05 & all_results_merged$fold_change > 1]))
}
ont_levels$n_associated_cells <- n_associated_cells
rm(n_associated_cells)
pal <- wesanderson::wes_palette("Darjeeling2", n=2)
ontlvl_ncells_plt <- ggplot(ont_levels, mapping=aes(x= factor(ont_lev), y=n_associated_cells)) +
  geom_jitter(color = pal[1]) +
  geom_violin(fill = NA) +
  geom_smooth(color = pal[2], method="loess",mapping = aes(x=ont_lev)) +
  cowplot::theme_cowplot() +
  labs(x="Ontology level", y="Associated cells/phenotype (n)",title = "n associated cells/phenotype")+
  blank_x
#ontlvl_ncells_plt
# ont level facet
signif_res <- all_results_merged[all_results_merged$q < 0.05,]
MammalianPhenotype <- read.delim("data/MammalianPhenotype.txt",header = FALSE)
colnames(MammalianPhenotype) = c("MPID","Phenotype")
signif_res <- merge(signif_res, MammalianPhenotype, by = "Phenotype", all.x = T)
signif_res[8880:8881,"MPID"] = "MP:0010574"
signif_res[10889:10890,"MPID"] = "MP:0012676"
signif_res[10898,"MPID"] = "MP:0013566"
signif_res[10924,"MPID"] = "MP:0000602"
signif_res[10928:10931,"MPID"] = "MP:0002705"
signif_res[13867,"MPID"] = "MP:0001183"

ontlevz <- c()
for (p in unique(signif_res$MPID)) {
  ontlevz[p] <- get_ont_level(ontology,p)
}
signif_res$ontlvl <- ontlevz[signif_res$MPID]
signif_res <- signif_res[complete.cases(signif_res),]
pal <- wesanderson::wes_palette("Darjeeling2", n=2)
ontlvl_fold_plt <- ggplot(signif_res, aes(x = factor(ontlvl), y = fold_change)) +
  geom_jitter(color = pal[1]) +
  geom_violin(fill = NA) +
  geom_smooth(color = pal[2], method="loess",mapping = aes(x=ontlvl)) +
  cowplot::theme_cowplot() +
  labs(x="Ontology level", y="Fold change",title = "Fold change in specific expression") +
  blank_x

# ngenes plt (stat smooth takes too long/doesnt work with all results so just use signif ?)
ngenes<-c()
for(p in unique(signif_res$Phenotype)) {
  ngenes[p] <- length(unique(gene_data$Gene[gene_data$Phenotype==p]))
}
signif_res$ngenes <- ngenes[signif_res$Phenotype]
signif_res <- signif_res[complete.cases(signif_res),]

pal <- wesanderson::wes_palette("Darjeeling2", n=2)
ontlvl_ngenes_plt <- ggplot(signif_res, aes(x = factor(ontlvl), y = ngenes)) +
  geom_jitter(color = pal[1]) +
  geom_violin(fill = NA) +
  geom_smooth(color = pal[2], method="loess",mapping = aes(x=ontlvl)) +
  cowplot::theme_cowplot() +
  labs(x="MP Ontology level", y="Genes (n)",title = "Number of genes")

###HPO 
gene_data = read.delim("data/phenotype_to_genes.txt", header = FALSE, skip =1)
colnames(gene_data) = c("HPID", "Phenotype", "EntrezID", "Gene", "Additional", "Source", "LinkID")
hpo <- get_ontology("data/hp.obo")
all_results_merged <- readRDS("data/results_2022-07-17 23-06-35.rds")
all_results_merged$Phenotype=gsub("-SLASH-","/",all_results_merged$Phenotype)
all_results_merged$CellType = gsub("_"," ",all_results_merged$CellType)
ont_levels_hpo <- data.frame("phenotype"=unique(gene_data$Phenotype),
                             "hpo_id"=hpo$id[match(unique(gene_data$Phenotype),hpo$name)])

get_ont_level = function(hpo,term_ids) {
  children = unique(setdiff(unlist(hpo$children[term_ids]), term_ids))
  if (length(children) == 0) {
    return(0)
  } else {
    return(1 + get_ont_level(hpo,children)) #<- recursion..
  }
}
ont_levels_hpo <- ont_levels_hpo[complete.cases(ont_levels_hpo),]
#ontlvl_ncells_plt
lvls <- c()
for (id in ont_levels_hpo$hpo_id) {
  lvls = append(lvls, get_ont_level(hpo,id))
}
ont_levels_hpo$ont_lev <- lvls
rm(lvls)
n_associated_cells <- c()
for (p in ont_levels_hpo$phenotype) {
  n_associated_cells <- append(n_associated_cells, length(all_results_merged$CellType[all_results_merged$Phenotype == p & all_results_merged$q < 0.05 & all_results_merged$fold_change > 1]))
}
ont_levels_hpo$n_associated_cells <- n_associated_cells
rm(n_associated_cells)
pal <- wesanderson::wes_palette("Darjeeling2", n=5)
ontlvl_ncells_plt1 <- ggplot(ont_levels_hpo, mapping=aes(x= factor(ont_lev), y=n_associated_cells)) +
  geom_jitter(color = pal[1]) +
  geom_violin(fill = NA) +
  geom_smooth(color = pal[2], method="loess",mapping = aes(x=ont_lev)) +
  cowplot::theme_cowplot() +
  labs(y="Associated cells/phenotype (n)",title = "n associated cells/phenotype")+
  blank_x
# ont level facet
signif_res_hpo <- all_results_merged[all_results_merged$q < 0.05,]
signif_res_hpo$HPID = NA
for( p in unique(signif_res_hpo$Phenotype)){
  signif_res_hpo[signif_res_hpo$Phenotype == p,"HPID"] = ont_levels_hpo[ont_levels_hpo$phenotype == p,"hpo_id"]
}

ontlevz <- c()
for (p in unique(signif_res_hpo$HPID)) {
  ontlevz[p] <- get_ont_level(hpo,p)
}
signif_res_hpo$ontlvl <- ontlevz[signif_res_hpo$HPID]
signif_res_hpo <- signif_res_hpo[complete.cases(signif_res_hpo),]
pal <- wesanderson::wes_palette("Darjeeling2", n=5)
ontlvl_fold_plt1 <- ggplot(signif_res_hpo, aes(x = factor(ontlvl), y = fold_change)) +
  geom_jitter(color = pal[1]) +
  geom_violin(fill = NA) +
  geom_smooth(color = pal[2], method="loess",mapping = aes(x=ontlvl)) +
  cowplot::theme_cowplot() +
  labs( y="Fold change",title = "Fold change in specific expression") +
  blank_x

# ngenes plt
ngenes<-c()
for(p in unique(signif_res_hpo$Phenotype)) {
  ngenes[p] <- length(unique(gene_data$Gene[gene_data$Phenotype==p]))
}
signif_res_hpo$ngenes <- ngenes[signif_res_hpo$Phenotype]
signif_res_hpo <- signif_res_hpo[complete.cases(signif_res_hpo),]

pal <- wesanderson::wes_palette("Darjeeling2", n=5)
ontlvl_ngenes_plt1 <- ggplot(signif_res_hpo, aes(x = factor(ontlvl), y = ngenes)) +
  geom_jitter(color = pal[1]) +
  geom_violin(fill = NA) +
  geom_smooth(color = pal[2], method="loess",mapping = aes(x=ontlvl)) +
  cowplot::theme_cowplot() +
  labs(x="HP ontology level", y="Genes (n)",title = "Number of genes")

ontlvlplots <- cowplot::plot_grid(plotlist=list(ontlvl_ncells_plt,ontlvl_ncells_plt1,ontlvl_fold_plt,ontlvl_fold_plt1,
                                                ontlvl_ngenes_plt,ontlvl_ngenes_plt1),align = "hv",nrow = 3,ncol = 2,labels=c("A","B","C","D","E","F"))



