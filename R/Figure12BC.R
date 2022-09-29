library(ggplot2)
library(cowplot)
signif_results_h <- read.delim("data/signif_results_h.txt",header = TRUE,sep=",")
signif_results_m <- read.delim("data/signif_results_Mondo.txt",header = TRUE,sep=",")
#signif_results_m = signif_results_m[(signif_results_m$DO.Disease %in% signif_results_h$Disease),]
signif_results_m = signif_results_m[(signif_results_m$MONDO_ID %in% signif_results_h$MONDO),]
mapping_all <- read.delim("data/ph_mp_mapping.txt",header = TRUE,sep=",")

mapping = data.frame()
for(d in unique(signif_results_m$MONDO_ID)){
  id_m = unique(signif_results_m[signif_results_m$MONDO_ID == d, "MPID"])
  id_h = unique(signif_results_h[signif_results_h$MONDO == d, "HPID"])
  m =0
  for(i in id_m){
    map_h= mapping_all[mapping_all$MPID == i,"HPID"]
    if(length(map_h != 0)){
    if(sum((map_h %in% id_h) == TRUE)>0){
      m =m+1
    }}}
  h = 0
  for(i in id_h){
    map_m= mapping_all[mapping_all$HPID == i,"MPID"]
    if(length(map_m != 0)){
    if(sum((map_m %in% id_m) == TRUE)>0){
      h =h+1
    }}}
  mapping = rbind(mapping, data.frame("MONDO_ID"=d,"MP"=m,"HP" = h,"n_MP"=length(id_m),
                                      "n_HP"=length(id_h),"p_m" = sprintf("%.4f", (m/length(id_m))),
                                      "p_h" = sprintf("%.4f", (h/length(id_h)))))
}
mapping$p_m <- as.numeric(mapping$p_m)
mapping$p_h <- as.numeric(mapping$p_h)

mapping_plot = data.frame()
mapping_plot = rbind(mapping_plot, data.frame("PO"=c(rep("MPO",length(mapping$p_m)),rep("HPO",length(mapping$p_h))),
                                            "Number" = c(mapping$MP,mapping$HP), 
                                            "percent" = c(mapping$p_m,mapping$p_h)))

mapping_plot = mapping_plot[-which(mapping_plot$Number == 0),]
mapping_plot$PO = factor(mapping_plot$PO,levels = c("MPO","HPO"))

Per_cros = ggplot(mapping_plot, aes(x = Number, y = percent, fill = PO)) +
  geom_point(pch = 21, position = position_jitterdodge(0.6),size = 1.3,alpha=0.9,colour = "white",stroke = 0.001)+
  ylab("Percent of cross-species phenotypes") +
  xlab("Number of cross-species phenotypes")+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(breaks=mapping_plot$Number, labels = mapping_plot$Number)+
  cowplot::theme_cowplot()+
  theme(legend.title = element_blank())


hist =data.frame()
for( n in unique(mapping$MP)){
  hist = rbind(hist, data.frame("PO"= "MPO", "Number" = n, "percent" = sprintf("%.4f", (length(mapping[mapping$MP == n,"MP"])/length(mapping$MP)))))
}
for( n in unique(mapping$HP)){
  hist = rbind(hist, data.frame("PO"= "HPO", "Number" = n, "percent" = sprintf("%.4f", (length(mapping[mapping$HP == n,"HP"])/length(mapping$HP)))))
}
hist$percent <- as.numeric(hist$percent)
hist$Number <- as.character(hist$Number)
hist$PO = factor(hist$PO,levels = c("MPO","HPO"))
hist$Number = factor(hist$Number,levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14"))
per_hist = ggplot(data=hist,aes(PO,percent,fill=Number))+  
  scale_y_continuous(labels = scales::percent)+
  geom_bar(stat="identity", position="stack",width=0.6,size=0.25) + 
  cowplot::theme_cowplot()+
  theme(axis.title.x = element_blank())+
  ylab("Percent of cross-species phenotype number")

mapping_per_h = plot_grid(mapping_plt,per_hist,align = "h",nrow = 1,labels = c("A","B")) 

plot_grid(mapping_per_h,Per_cros,nrow = 2,rel_heights = c(1,1),labels = c("","C")) 
