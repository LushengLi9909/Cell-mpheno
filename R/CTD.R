library(stringr)
library(readxl)

FACS_URL = "https://ndownloader.figshare.com/files/10700143" 
annot_URL = "https://ndownloader.figshare.com/files/13088129"


if (!file.exists(file.path("data/FACS"))){ 
  download.file(FACS_URL, "data/FACS.zip", mode = "wb")
  facs = unzip("data/FACS.zip",  exdir=paste0(getwd(), "/data")); unlink("__MACOSX", recursive = TRUE)  #to delete an extra file in the tm link
}
if (!file.exists("data/annotations_facs.csv")){
  download.file(annot_URL, "data/annotations_facs.csv")}

f_names <- dir(path = "FACS/FACS/") #Requires filepath "FACS" containing all tm FACS data
annot <- read.csv("annotations_facs.csv") #Found on the original upload for the Tabula Muris data


levels(annot$cell_ontology_class) <- c(levels(annot$cell_ontology_class), "NA")
annot$cell_ontology_class[which(annot$cell_ontology_class == "")] <- "NA" #Unlabeled cells -> NA


for(i in 1:length(f_names)){
  cat("\r",paste("Files remaining:", length(f_names)-i, "    ")) 
  tissue <- str_sub(f_names[i], end = -12)
  annot1 <- annot[which(annot$tissue == tissue), ]
  data1 <- read.csv(paste("FACS/FACS/", f_names[i], sep = ""))
  rownames(data1) <- data1$X
  data1 <- data1[, -1]
  data2 <- data1[, as.character(annot1$cell)]
 
  
  uniquel2 <- unique(annot1$cell_ontology_class)
  x <- data.frame(data2[, 1:length(uniquel2)])
  names(x) <- paste(uniquel2, tissue, sep = "_")
  
  
  for(k in 1:length(uniquel2)){
    cid <- which(annot1$cell_ontology_class == uniquel2[k])
    x[, k] <- rowMeans(data2[, cid])
  }
  if(i == 1){
    y = x
  }
  else{
    y = cbind(y, x)
  }
}
L1_classifications_xlsx_path = "data/tm_level1classifications_full.xlsx"

Level1Data <- read_excel(L1_classifications_xlsx_path)

#Generating annotation labels by making a new annot column representing each individual cell's level 1 class.
annot$cell_ontology_class_with_tissue <- annot$cell_ontology_class
for(i in 1:length(annot$cell_ontology_class_with_tissue)){
  annot$cell_ontology_class_with_tissue[i] <- paste(annot$cell_ontology_class_with_tissue[i], annot$tissue[i], sep = "_")
}

annot$cell_ontology_class_l1 <- annot$cell_ontology_class_with_tissue

for(i in 1:length(Level1Data$Level1Classification)){
  currentCells <- Level1Data$L1Combinations[i]
  currentCells <- unlist(strsplit(currentCells, ", "))
  annot$cell_ontology_class_l1[which(annot$cell_ontology_class_l1 %in% currentCells)] <- Level1Data$Level1Classification[i]
}
annot_levels = annot[,c("cell_ontology_class_with_tissue","cell_ontology_class_l1")]
annot_levels = annotlevels[!duplicated(annotlevels),]
level2 = colnames(y)
level1 = c()
for( i in 1:length(level2)){
  if(level2[i] %in% annot_levels$cell_ontology_class_with_tissue ){
    level1 = append(level1,annot_levels$cell_ontology_class_l1[which(annot_levels$cell_ontology_class_with_tissue %in% level2[i])])
  }
}
annotLevels <- list(l1 = level1, l2 = level2)

ctd_result1 <- EWCE::generate_celltype_data(exp = y,
                                           groupName = "TabulaMuris",
                                           annotLevels = annotLevels,
                                           input_species = "mouse",
                                           output_species = "human",
                                           method = "homologene",
                                           as_sparse = TRUE,
                                           convert_orths = TRUE,
                                           dendrograms = TRUE,
                                           return_ctd = TRUE)
saveRDS(ctd_result1$ctd,file ="TabulaMuris_n.rds")

