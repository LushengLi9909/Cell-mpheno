library(ggplot2)
library(cowplot)
library(ontologyIndex)

ctd <- readRDS("data/TabulaMuris_n.rds")
cell_mappings = read.csv("data/mammalian_celltype_mapping.csv")
colnames(cell_mappings) = c("tissue", "level1", "level2") 
cell_mappings$level1 = gsub("\\."," ",cell_mappings$level1)

all_results_merged_mpo <- readRDS("data/mpo_results_merged.rds")
all_results_merged_hpo <- readRDS("data/results_2022-07-17 23-06-35.rds")
all_results_merged_hpo$Phenotype=gsub("-SLASH-","/",all_results_merged_hpo$Phenotype)
all_results_merged_hpo$CellType = gsub("_"," ",all_results_merged_hpo$CellType)

ontology_hpo <- get_ontology("data/hp.obo")
ontology <- get_ontology("data/MPheno_OBO.ontology")
plot_branches = c("nervous system phenotype","Abnormality of the nervous system",
                  "cardiovascular system phenotype","Abnormality of the cardiovascular system",
                  "immune system phenotype","Abnormality of the immune system")

plot_phenos_per_cell <- function (all_results_merged, cell_mappings, fold = 1, q_val = 0.05) 
{
  phenos_per_cell = data.frame()
  for (c in unique(all_results_merged$CellType)) {
    cur_cell = c
    n_phenos = length(all_results_merged[all_results_merged$CellType == 
                                           c & all_results_merged$q < q_val & all_results_merged$fold_change > 
                                           fold, "Phenotype"])
    phenos_per_cell = rbind(phenos_per_cell, data.frame(Cell = cur_cell, 
                                                        n_phenos = n_phenos))
  }

  n_pheno_text = phenos_per_cell
  phenos_per_cell$Tissue = rep(NA, length(phenos_per_cell$Cell))
  
    phenos_per_cell_tissue = data.frame()
    for (i in seq(1:length(phenos_per_cell$Cell))) {
      cur = cell_mappings[cell_mappings$level1 == paste(phenos_per_cell$Cell[i]), 
      ]$tissue
      for (t in unique(cur)) {
        if (!is.na(t)) {
          phenos_per_cell$Tissue[i] = t
          phenos_per_cell_tissue = rbind(phenos_per_cell_tissue, 
                                         phenos_per_cell[i, ])
        }
      }
    }
    
    phenos_per_cell = phenos_per_cell_tissue
  
  
    for (i in seq(1, length(phenos_per_cell$Cell))) {
      cur = phenos_per_cell[phenos_per_cell$Cell == phenos_per_cell$Cell[i], ]
      phenos_per_cell$n_phenos[i] = phenos_per_cell$n_phenos[i]/length(unique(cur$Tissue))
    }

    
    ctd[[1]]$plotting$cell_ordering = factor(ctd[[1]]$plotting$cell_ordering,levels = ctd[[1]]$plotting$cell_ordering)
    phenos_per_cell$Cell = factor(phenos_per_cell$Cell,levels = ctd[[1]]$plotting$cell_ordering)
    
  explan_plt <- ggplot(phenos_per_cell, aes(x=Cell, y=n_phenos)) + 
    scale_y_continuous("", expand = c(0, 0), limits = c(0, max(n_pheno_text$n_phenos) + 55)) + theme_classic() + 
    scale_x_discrete(breaks = NULL) +
    labs(title = paste0("Significantly enriched phenotypes per cell (q < ", q_val, ", fold change > ", fold,")" ))
    
    explan_plt = explan_plt + geom_col(aes(fill = Tissue), size = 2, position = "stack")

  explan_plt = explan_plt + geom_text(data = n_pheno_text, 
                                      aes(label = n_phenos, x = Cell, y = n_phenos + 26), size = 4, 
                                      color = "black")  + theme(axis.title.x = element_blank(), axis.text.y = element_text(colour = "black"),legend.key.size = unit(16,"pt"))+
    labs(x = NULL) 
  return(explan_plt)
}

plot_n_signif_phenos_per_cell_by_branch <- function(all_results_merged_mpo,all_results_merged_hpo,
                                                    plot_branches = c("nervous system phenotype",
                                                                      "cardiovascular system phenotype",
                                                                      "immune system phenotype"),
                                                    q_threshold = 0.05) {
  
  # Get main mpo branches (child nodes of PHenotypic abnormality)
  
  mpo_branches = ontology$children["MP:0000001"]
  mpo_branches_names = c()
  for (b in mpo_branches){
    mpo_branches_names = c(mpo_branches_names, ontology$name[b])
  }
  main_branch_df = data.frame("id"=mpo_branches,"phenotype"=mpo_branches_names)
  colnames(main_branch_df) = c("id","phenotype")
  hpo_branches = ontology_hpo$children["HP:0000118"]
  hpo_branches_names = c()
  for (b in hpo_branches){
    hpo_branches_names = c(hpo_branches_names, ontology_hpo$name[b])
  }
  main_branch_df_hpo = data.frame("id"=hpo_branches,"phenotype"=hpo_branches_names)
  colnames(main_branch_df_hpo) = c("id","phenotype")
  main_branch_df = rbind(main_branch_df,main_branch_df_hpo)
  descendants = list ()
  for (b in mpo_branches[[1]] ) {
    print(b)
    print(ontology$name[b])
    descendants[[ontology$name[b]]] = get_descendants(ontology,b)
  }
  descendants_hpo = list ()
  for (b in hpo_branches[[1]] ) {
    print(b)
    print(ontology_hpo$name[b])
    descendants_hpo[[ontology_hpo$name[b]]] = get_descendants(ontology_hpo,b)
  }
  # Label all results based on branch
  if (! "branch" %in% colnames(all_results_merged_mpo)) {
    phenotypes = unique(all_results_merged_mpo$Phenotype)
    all_results_merged_mpo$branch = NA
    pheno_branches = c()
    for (p in phenotypes){
      print(p)
      index = match(p,ontology$name)
      if (!is.na(index)){
        id = ontology$id[[index]]
        for (i in seq(1,length(descendants))) {
          if (id %in% descendants[[i]]) {
            pheno_branches[p] = names(descendants[i])
          }
        }
      }
    }
    
    all_results_merged_mpo$branch = NA
    #remaining = length(unique(all_results_merged$Phenotype))
    #cat("Assigning branch to phenotypes\nN remaining:\n")
    
    #cat("\r",remaining,"      ")
    #remaining = remaining - 1
    
    for (p in unique(all_results_merged_mpo$Phenotype)) {
      all_results_merged_mpo$branch[all_results_merged_mpo$Phenotype == p] = pheno_branches[p]
    }
  }
  
  if (! "branch" %in% colnames(all_results_merged_hpo)) {
    phenotypes = unique(all_results_merged_hpo$Phenotype)
    all_results_merged_hpo$branch = NA
    pheno_branches_hpo = c()
    for (p in phenotypes){
      print(p)
      index = match(p,ontology_hpo$name)
      if (!is.na(index)){
        id = ontology_hpo$id[[index]]
        for (i in seq(1,length(descendants_hpo))) {
          if (id %in% descendants_hpo[[i]]) {
            pheno_branches_hpo[p] = names(descendants_hpo[i])
          }
        }
      }
    }
    
    all_results_merged_hpo$branch = NA
    #remaining = length(unique(all_results_merged$Phenotype))
    #cat("Assigning branch to phenotypes\nN remaining:\n")
    
    #cat("\r",remaining,"      ")
    #remaining = remaining - 1
    
    for (p in unique(all_results_merged_hpo$Phenotype)) {
      all_results_merged_hpo$branch[all_results_merged_hpo$Phenotype == p] = pheno_branches_hpo[p]
    }
  }
  #for(i in 1:length(all_results_merged$Phenotype)) {
  #  if(all_results_merged[i,8]=="mammalian phenotype"){
  #    all_results_merged[i,9] <- " "
  # }}
  # plot it
  signif_results = all_results_merged_mpo[all_results_merged_mpo$q <= q_threshold & !is.na(all_results_merged_mpo$branch), ]
  n_signif_per_cell_by_branch = data.frame()
  
  plot_branches_mpo = plot_branches[plot_branches %in% mpo_branches_names]
  for (b in unique(plot_branches_mpo)) {
    for (c in unique(all_results_merged_mpo$CellType)){
      n = length(signif_results[signif_results$CellType == c & signif_results$branch == b,]$CellType)
      n_signif_per_cell_by_branch = rbind(n_signif_per_cell_by_branch,
                                          data.frame("branch" = b,"CellType"=c,"n_signif" = n))
    }
  }
  
  
  # Run hypergeometric tests
  # To show which cell types within a branch have more enriched phenotypes than expected
  n_signif_per_cell_by_branch$hypergeo_p = NA
  total_signif = length(signif_results$q)
  
  for (i in seq(1, length(n_signif_per_cell_by_branch$CellType))) {
    cell = n_signif_per_cell_by_branch$CellType[i]
    branch = n_signif_per_cell_by_branch$branch[i]
    group1 = sum(n_signif_per_cell_by_branch$n_signif[n_signif_per_cell_by_branch$CellType == cell])
    group2 = sum(n_signif_per_cell_by_branch$n_signif[n_signif_per_cell_by_branch$branch == branch])
    overlap = n_signif_per_cell_by_branch$n_signif[i]
    n_signif_per_cell_by_branch$hypergeo_p[i] = phyper(overlap - 1, group2, total_signif - group2, group1, lower.tail=FALSE)
  }
  n_signif_per_cell_by_branch$hypergeo_q = p.adjust(n_signif_per_cell_by_branch$hypergeo_p,method = "BH")
  
  n_signif_per_cell_by_branch$signif_asterics = ""
  n_signif_per_cell_by_branch$signif_asterics[n_signif_per_cell_by_branch$hypergeo_q<0.05] = "*"
  n_signif_per_cell_by_branch$signif_asterics[n_signif_per_cell_by_branch$hypergeo_q<0.005] = "**"
  n_signif_per_cell_by_branch$signif_asterics[n_signif_per_cell_by_branch$hypergeo_q<0.0005] = "***"
  n_signif_per_cell_by_branch$signif_asterics[n_signif_per_cell_by_branch$hypergeo_q<0.00005] = "****"
  
  plot_branches_hpo = plot_branches[plot_branches %in% hpo_branches_names]
  # plot it
  signif_results_hpo = all_results_merged_hpo[all_results_merged_hpo$q <= q_threshold & !is.na(all_results_merged_hpo$branch), ]
  n_signif_per_cell_by_branch_hpo = data.frame()
  for (b in unique(plot_branches_hpo)) {
    for (c in unique(all_results_merged_hpo$CellType)){
      n = length(signif_results_hpo[signif_results_hpo$CellType == c & signif_results_hpo$branch == b,]$CellType)
      n_signif_per_cell_by_branch_hpo = rbind(n_signif_per_cell_by_branch_hpo,
                                              data.frame("branch" = b,"CellType"=c,"n_signif" = n))
    }
  }
  
  
  # Run hypergeometric tests
  # To show which cell types within a branch have more enriched phenotypes than expected
  n_signif_per_cell_by_branch_hpo$hypergeo_p = NA
  total_signif = length(signif_results_hpo$q)
  
  for (i in seq(1, length(n_signif_per_cell_by_branch_hpo$CellType))) {
    cell = n_signif_per_cell_by_branch_hpo$CellType[i]
    branch = n_signif_per_cell_by_branch_hpo$branch[i]
    group1 = sum(n_signif_per_cell_by_branch_hpo$n_signif[n_signif_per_cell_by_branch_hpo$CellType == cell])
    group2 = sum(n_signif_per_cell_by_branch_hpo$n_signif[n_signif_per_cell_by_branch_hpo$branch == branch])
    overlap = n_signif_per_cell_by_branch_hpo$n_signif[i]
    n_signif_per_cell_by_branch_hpo$hypergeo_p[i] = phyper(overlap - 1, group2, total_signif - group2, group1, lower.tail=FALSE)
  }
  n_signif_per_cell_by_branch_hpo$hypergeo_q = p.adjust(n_signif_per_cell_by_branch_hpo$hypergeo_p,method = "BH")
  
  n_signif_per_cell_by_branch_hpo$signif_asterics = ""
  n_signif_per_cell_by_branch_hpo$signif_asterics[n_signif_per_cell_by_branch_hpo$hypergeo_q<0.05] = "*"
  n_signif_per_cell_by_branch_hpo$signif_asterics[n_signif_per_cell_by_branch_hpo$hypergeo_q<0.005] = "**"
  n_signif_per_cell_by_branch_hpo$signif_asterics[n_signif_per_cell_by_branch_hpo$hypergeo_q<0.0005] = "***"
  n_signif_per_cell_by_branch_hpo$signif_asterics[n_signif_per_cell_by_branch_hpo$hypergeo_q<0.00005] = "****"
  # Ordering cell types by dendrogram grouping
  
  cell_order = factor(ctd[[1]]$plotting$cell_ordering,levels = ctd[[1]]$plotting$cell_ordering)
  n_signif_per_cell_by_branch_new = rbind(n_signif_per_cell_by_branch, n_signif_per_cell_by_branch_hpo)
  n_signif_per_cell_by_branch_new$CellType = factor(n_signif_per_cell_by_branch_new$CellType,levels = cell_order)
  n_signif_per_cell_by_branch_new$branch = factor(n_signif_per_cell_by_branch_new$branch,levels =plot_branches)
  
  facet_branch_plt <- ggplot(n_signif_per_cell_by_branch_new, aes(x=CellType,y=n_signif,fill=branch)) +
    geom_col() +
    geom_text(mapping= aes(label = signif_asterics, y = n_signif + 5))+
    theme_cowplot()+
    ylab("Number of significantly enriched phenotypes") +labs(x = NULL)+
    theme(axis.text.x = element_text(angle = 90, hjust=1,vjust = 0.3),legend.position="none")+ 
    #scale_x_discrete(breaks =cell_order)+
    facet_wrap(~branch,ncol=1)  
  return(facet_branch_plt)
}


the_dendrogram = ctd[[1]]$plotting$ggdendro_horizontal +
  theme(plot.margin = unit(c(0, 0, 0, 0), units = "cm")) +
  scale_x_discrete(breaks = ctd[[1]]$plotting$cell_ordering)
# CHANGE ^: added scale_x_discrete to set the mapping of the dendro to the x axis scale

explan_plt1 = plot_phenos_per_cell(all_results_merged_mpo, cell_mappings, fold = 1, q_val = 0.05)
explan_plt2 = plot_phenos_per_cell(all_results_merged_hpo, cell_mappings, fold = 1, q_val = 0.05)
explan_plt2 = explan_plt2 + guides(fill=F)+ labs(title = NULL)

facet_branch_plt = plot_n_signif_phenos_per_cell_by_branch(all_results_merged_mpo,all_results_merged_hpo,
                                                           plot_branches = c("nervous system phenotype","Abnormality of the nervous system",
                                                                             "cardiovascular system phenotype","Abnormality of the cardiovascular system",
                                                                             "immune system phenotype","Abnormality of the immune system"),
                                                           q_threshold = 0.05)

figure1 = cowplot::plot_grid(the_dendrogram,explan_plt1,explan_plt2,facet_branch_plt,axis = "lr",
                                 align = "v", ncol = 1,
                                 rel_heights = c(0.13,1.2,1,2),labels=c("A","B","C","D"))
