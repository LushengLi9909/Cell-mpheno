library(ontologyIndex)
library(ggplot2)
#color_pal = wes_palette("Darjeeling1",5)
color_pal = RColorBrewer::brewer.pal(9,"Set1")
color = c()
color[1]= color_pal[2]
color[2]= color_pal[1]
color[3]= color_pal[3]

all_results_merged_mpo <- readRDS("data/mpo_results_merged.rds")
all_results_merged_mpo$minus_log_p = round(-log10(all_results_merged_mpo$p+0.000001),digits=0)

all_results_merged_hpo <- readRDS("data/results_2022-07-17 23-06-35.rds")
all_results_merged_hpo$Phenotype=gsub("-SLASH-","/",all_results_merged_hpo$Phenotype)
all_results_merged_hpo$CellType = gsub("_"," ",all_results_merged_hpo$CellType)
all_results_merged_hpo$minus_log_p = round(-log10(all_results_merged_hpo$p+0.000001),digits=0)

ontology_mpo <- get_ontology("data/MPheno_OBO.ontology")
ontology_hpo = get_ontology("data/hp.obo")
blank_x_axis = c(TRUE, TRUE, FALSE)
proportion_plots = list()

proportion_enrich = function(all_results_merged,ontology,plot_branches,expected_cells){
for (i in seq(1,length(plot_branches))) {
  target_cell = expected_cells[i]
  cur_branch = plot_branches[i]
  cur_id = ontology$id[match(cur_branch, ontology$name)]
  cur_descendants = get_descendants(ontology, cur_id)
  cur_descendants_names = c()
  for(d in cur_descendants) {
    cur_descendants_names = append(cur_descendants_names, paste(ontology$name[d]))
  }
  all_results_merged$cur_branch = "Other"
  all_results_merged[all_results_merged$Phenotype %in% cur_descendants_names,]$cur_branch = cur_branch
  all_results_merged$cur_branch = factor(all_results_merged$cur_branch)
  #all_results_merged$cur_branch = factor(all_results_merged$cur_branch, levels=c(cur_branch,"Other"))
  proportional_results = data.frame()
  for (logp in unique(all_results_merged$minus_log_p)) {
    cur = all_results_merged[all_results_merged$minus_log_p == logp, ]
    cur$CellType = gsub("mouse_TabulaMuris_","",cur$CellType)
    cur$CellType = gsub("_"," ",cur$CellType)
    prop_branch = 100 * (length(cur$cur_branch[cur$cur_branch==cur_branch & cur$CellType == target_cell])/length(cur$cur_branch[cur$CellType==target_cell]))
    prop_other = 100 * (length(cur$cur_branch[cur$cur_branch == "Other" & cur$CellType == target_cell])/length(cur$cur_branch[cur$CellType==target_cell]))
    proportional_results = rbind(proportional_results,
                                 data.frame("cur_branch" = cur_branch,
                                            "prop" = prop_branch,
                                            "minus_log_p" = logp))
    proportional_results = rbind(proportional_results,
                                 data.frame("cur_branch" = "Other",
                                            "prop" = prop_other,
                                            "minus_log_p" = logp))
  }
  prop_plot <- ggplot(proportional_results, aes(x = minus_log_p, 
                                                y = prop, color = cur_branch)) + geom_point(size = 3) + 
    geom_line(mapping = aes(linetype = cur_branch), size = 1) + 
    labs(color = "", linetype = "") + ylab(paste0("% Enrichments in \"", 
                                                  target_cell, "\"")) +
    scale_y_continuous(limits = c(0, 100)) + cowplot::theme_cowplot() + 
    scale_color_manual(values = c(color[i], 
                                  color_pal[9])) + theme(legend.position = c(0.5,0.95), 
                                                         legend.text = element_text(size = 10), title = element_text(size = 11), axis.title.y = element_text(size=12))
  if (! blank_x_axis[i]){
    prop_plot = prop_plot + xlab("-log10(p)")
  } else {
    prop_plot = prop_plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  } 
  proportion_plots[[i]] = prop_plot
  
}
  return(proportion_plots)}
prop = proportion_enrich(all_results_merged = all_results_merged_mpo,ontology = ontology_mpo,
                             plot_branches = c("nervous system phenotype", "cardiovascular system phenotype", 
                                               "immune system phenotype"),
                             expected_cells = c("Neurons", "Cardiac muscle cells","Natural killer cells"))
  proportion = list()
  proportion[[1]] = prop[[1]]
  proportion[[3]] = prop[[2]]
  proportion[[5]] = prop[[3]]
  prop = proportion_enrich(all_results_merged = all_results_merged_hpo,ontology = ontology_hpo,
                           plot_branches = c("Abnormality of the nervous system", "Abnormality of the cardiovascular system",
                                             "Abnormality of the immune system"),
                           expected_cells = c("Neurons", "Cardiac muscle cells","Natural killer cells"))  
  proportion[[2]] = prop[[1]]
  proportion[[4]] = prop[[2]]
  proportion[[6]] = prop[[3]]

print(cowplot::plot_grid(plotlist=proportion, align = "h",nrow = 3,ncol = 2))
#,labels=c("A","B","C")

