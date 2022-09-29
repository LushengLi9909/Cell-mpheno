library(wesanderson)
library(ontologyIndex)
color_pal = RColorBrewer::brewer.pal(9,"Set1")
hpo <- get_ontology("data/hp.obo")
plot_branches = c("Abnormality of the nervous system", "Abnormality of the cardiovascular system",
                  "Abnormality of the immune system")
plot_branches_ids = hpo$id[match(plot_branches, hpo$name)]
hpo_branches = hpo$children["HP:0000118"]
phenos_per_branch = data.frame()
for (b in hpo_branches[[1]]) {
  n = length(get_descendants(hpo,b))
  if (b %in% plot_branches_ids) {
    target_branch = b
  } else {
    target_branch = "Other"
  }
  phenos_per_branch = rbind(phenos_per_branch, data.frame("branch"=hpo$name[b],"n_phenos"=n, "target"=target_branch))
}
phenos_per_branch$branch = reorder(phenos_per_branch$branch, phenos_per_branch$n_phenos)
phenos_per_branch_plt_hpo <- ggplot(phenos_per_branch, aes(x=n_phenos, y=branch, color= target, fill=target))+
  geom_segment(mapping= aes(xend = 0, yend = branch),  size = 3)+
  xlab(element_blank()) +
  ylab(element_blank())+
  scale_color_manual(values=c(color_pal[1],color_pal[3],
                              color_pal[2],color_pal[9])) +
  theme(axis.line.x=element_blank(),panel.background=element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),axis.line.y=element_line(color = "black"),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),legend.position = "none")
phenos_per_branch_plt_hpo = phenos_per_branch_plt_hpo + geom_text(data=phenos_per_branch,aes(label = n_phenos, x = n_phenos+ 250, y = branch), size = 3, color = "black")
print(phenos_per_branch_plt_hpo)
