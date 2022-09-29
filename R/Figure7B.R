library(ggplot2)
all_results_merged_mpo <- readRDS("data/mpo_results_merged.rds")
all_results_merged_hpo <- readRDS("data/results_2022-07-17 23-06-35.rds")
all_results_merged_hpo$Phenotype=gsub("-SLASH-","/",all_results_merged_hpo$Phenotype)
all_results_merged_hpo$CellType = gsub("_"," ",all_results_merged_hpo$CellType)

mapping_all <- read.delim("data/ph_mp_mapping.txt",header = TRUE,sep=",")

gene_data_hpo = read.delim("data/phenotype_to_genes.txt", header = FALSE, skip =1)
colnames(gene_data_hpo) = c("HPID", "Phenotype", "EntrezID", "Gene", "Additional", "Source", "LinkID")
gene_data_mpo <- read.delim("data/Phenotype_MPO.txt",header = TRUE,sep=",")

mapping_all$label = NA
for(i in 1:length(mapping_all$HPterm)){
pheno_m = all_results_merged_mpo[all_results_merged_mpo$Phenotype==mapping_all$MPterm[i] & all_results_merged_mpo$q <= 0.05,]
min_m = pheno_m[which.min(pheno_m$q),"q"]
pheno_min_m = pheno_m[pheno_m$q==min_m,"CellType"]

pheno_h = all_results_merged_hpo[all_results_merged_hpo$Phenotype==mapping_all$HPterm[i] & all_results_merged_hpo$q <= 0.05,]
min_h = pheno_h[which.min(pheno_h$q),"q"]
pheno_min_h = pheno_h[pheno_h$q==min_h,"CellType"]
if(sum((pheno_min_m %in% pheno_min_h) == TRUE)>0){
  mapping_all$label[i] = 1
} else{
  mapping_all$label[i] = 0
}
}
#HPOsterms have a cell type association, as the equivalent MPO term's most significant cell type
mapping_all$label_m = NA
for(i in 1:length(mapping_all$HPterm)){
  pheno_m = all_results_merged_mpo[all_results_merged_mpo$Phenotype==mapping_all$MPterm[i] & all_results_merged_mpo$q <= 0.05,]
  
  pheno_min_m = pheno_m[pheno_m$q<= 0.05,"CellType"]
  
  pheno_h = all_results_merged_hpo[all_results_merged_hpo$Phenotype==mapping_all$HPterm[i] & all_results_merged_hpo$q <= 0.05,]
  min_h = pheno_h[which.min(pheno_h$q),"q"]
  pheno_min_h = pheno_h[pheno_h$q==min_h,"CellType"]
  if(sum((pheno_min_m %in% pheno_min_h) == TRUE)>0){
    mapping_all$label_m[i] = 1
  } else{
    mapping_all$label_m[i] = 0
  }
}
mapping_all$label_h = NA
for(i in 1:length(mapping_all$HPterm)){
  pheno_m = all_results_merged_mpo[all_results_merged_mpo$Phenotype==mapping_all$MPterm[i] & all_results_merged_mpo$q <= 0.05,]
  min_m = pheno_m[which.min(pheno_m$q),"q"]
  pheno_min_m = pheno_m[pheno_m$q==min_m,"CellType"]
  
  pheno_h = all_results_merged_hpo[all_results_merged_hpo$Phenotype==mapping_all$HPterm[i] & all_results_merged_hpo$q <= 0.05,]
  
  pheno_min_h = pheno_h[pheno_h$q<=0.05,"CellType"]
  if(sum((pheno_min_m %in% pheno_min_h) == TRUE)>0){
    mapping_all$label_h[i] = 1
  } else{
    mapping_all$label_h[i] = 0
  }
}

#significant cell type association
mapping_all$label_a = NA
for(i in 1:length(mapping_all$HPterm)){
  pheno_m = all_results_merged_mpo[all_results_merged_mpo$Phenotype==mapping_all$MPterm[i] & all_results_merged_mpo$q <= 0.05,]
  
  pheno_min_m = pheno_m[pheno_m$q<= 0.05,"CellType"]
  
  pheno_h = all_results_merged_hpo[all_results_merged_hpo$Phenotype==mapping_all$HPterm[i] & all_results_merged_hpo$q <= 0.05,]
  
  pheno_min_h = pheno_h[pheno_h$q<=0.05,"CellType"]
  if(sum((pheno_min_m %in% pheno_min_h) == TRUE)>0){
    mapping_all$label_a[i] = 1
  } else{
    mapping_all$label_a[i] = 0
  }
}


same_most_sign = length(which(mapping_all$label ==1))
most_sign_m = length(which(mapping_all$label_m ==1))
most_sign_h = length(which(mapping_all$label_h ==1))
same_sign = length(which(mapping_all$label_a ==1))
Analysed_pheno_m = length(unique(all_results_merged_mpo$Phenotype))
Analysed_pheno_mapping_m = length(unique(mapping_all$MPterm)) - length(which((mapping_all$MPterm %in% all_results_merged_mpo$Phenotype) == FALSE))
Analysed_pheno_h = length(unique(all_results_merged_hpo$Phenotype))
Analysed_pheno_mapping_h = length(unique(mapping_all$HPterm)) - length(which((mapping_all$HPterm %in% all_results_merged_hpo$Phenotype) == FALSE))
per_mapping = data.frame()
per_mapping = rbind(per_mapping, data.frame("PO"=c("MPO","HPO","MPO","HPO","MPO","HPO"),
                                              "percent" = c(sprintf("%.4f", (same_most_sign/Analysed_pheno_mapping_m)),
                                                            sprintf("%.4f", (same_most_sign/Analysed_pheno_mapping_h)),
                                                            sprintf("%.4f", (most_sign_m/Analysed_pheno_mapping_m)),
                                                            sprintf("%.4f", (most_sign_h/Analysed_pheno_mapping_h)),
                                                            sprintf("%.4f", (same_sign/Analysed_pheno_mapping_m)),
                                                            sprintf("%.4f", (same_sign/Analysed_pheno_mapping_h))),
                                            "Type"= c("most significant cell type","most significant cell type",
                                                      "most/ significant cell type","most/ significant cell type",
                                                      "significant cell type","significant cell type")))

per_mapping$percent <- as.numeric(per_mapping$percent)

per_mapping$labels = paste0(per_mapping$percent*100, "%")
per_mapping$PO = factor(per_mapping$PO,levels = c("MPO","HPO"))
values = c("#DB423E","#008ECA")
mapping_plot = ggplot(data=per_mapping,aes(PO,percent,fill=Type))+  
  #geom_bar(stat="identity",color="black", width=0.6,size=0.25,fill=c("#F8766D","#00BFC4"))+
    geom_bar(stat="identity", position="dodge",color="black",width=0.8,size=0.25)+  
#scale_fill_manual(values = c("#DB423E","#008ECA"))+  
  scale_y_continuous(labels = scales::percent, limits=c(0,1))+
 ylab("Percent of the same significant cell type association")+
  geom_text(aes(label=labels,y = percent + 0.01),position = position_dodge(0.9),vjust = 0,size=3)+
  theme_classic()+
  theme(    axis.title.x = element_blank(), axis.text = element_text(size=12,face="plain",color="black"),
            legend.position = c(0.6, 0.9),legend.text = element_text( size=12 ))#+
#guides(fill=guide_legend(title=NULL)) ,    legend.title=element_text(size=12,face="plain",color="black") , legend.title= element_blank(),  legend.position = "top"  
plot_grid(mapping_gene_plot,mapping_plot,labels = c("A","B")) 
