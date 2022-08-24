library(EWCE)
set.seed(1234)
get_gene_list <- function(list_name,
                          gene_data,
                          list_name_column = "Phenotype",
                          gene_column = "Gene"){
  if (!list_name %in% unique(gene_data[,list_name_column])) {
    warning(paste("gene list", list_name, "is not present in" ,list_name_column, "column"))
  }
  return(paste(gene_data[,gene_column][gene_data[,list_name_column] == list_name]))
}
get_unfinished_list_names <- function (list_names, results_dir) {
  list_names_2 = c()
  for (l in list_names) {
    if (is_not_analysed(l,results_dir)) {
      list_names_2 <- append(list_names_2,l)
    }
  }
  return(list_names_2)
}
is_not_analysed <- function(list_name,results_dir) {
  list_name = gsub("/", "-SLASH-", list_name, fixed = TRUE)
  file_path <- paste0(results_dir,"/",list_name,".rds")
  if (file.exists(file_path)){
    return (FALSE)
  } else {
    return (TRUE)
  }
}
get_valid_gene_lists <- function(ctd,
                                 list_names,
                                 gene_data,
                                 list_name_column = "Phenotype",
                                 gene_column = "Gene"){
  ctd_genes = rownames(ctd[[1]]$specificity_quantiles)
  validLists = c()
  for (p in list_names) {
    if (sum(unique(get_gene_list(p,gene_data,list_name_column, gene_column)) %in% ctd_genes) >= 4) {
      validLists = append(validLists, p)
    }
  }
  return(validLists)
}

merge_results <- function(results_dir = "results", list_name_column = "phenotype") {
  results_merged <- data.frame()
  for (f in list.files(results_dir)) {
    cur = readRDS(paste0(results_dir,"/",f))$results
    descr <-strsplit(f,".rds")
    cur[,list_name_column] <- descr
    results_merged <- rbind(results_merged,cur)
  }
  return(results_merged)
}
gen_results <- function(ctd,
                        gene_data,
                        list_names,
                        background_genes,
                        list_name_column = "Phenotype",
                        gene_column = "Gene",
                        results_dir = "results",
                        overwrite_past_analysis = FALSE,
                        reps,
                        annotLevel,
                        genelistSpecies = "human",
                        sctSpecies = "human",
                        no_cores,
                        MergeResults = FALSE) {
  
  # remove gene lists that do not have enough valid genes (>= 4)
  validgenelists <- get_valid_gene_lists(ctd,
                                         list_names = list_names,
                                         gene_data = gene_data,
                                         list_name_column = list_name_column,
                                         gene_column = gene_column)
  
  # Create results directory and remove finished gene lists
  if (!file.exists(results_dir)) {
    dir.create(results_dir)
  }
  if (!overwrite_past_analysis) {
    validgenelists <- get_unfinished_list_names(validgenelists,results_dir)
  }
  
  # Run analysis
    
    results_list <- lapply(validgenelists,FUN=function(p,
                                                       gene_associations= gene_data,
                                                       lst_nm_col = list_name_column,
                                                       gn_col = gene_column,
                                                       results_directory = results_dir,
                                                       ctd_file = ctd,
                                                       background = background_genes,
                                                       bootstrap_reps=reps,
                                                       annotation_Level = annotLevel,
                                                       genes_Species = genelistSpecies,
                                                       ctd_Species = sctSpecies,
                                                       number_cores  = no_cores
                                                       ){
      print(p)
      genes = get_gene_list(p,gene_associations,lst_nm_col, gn_col)  
      p = gsub("/", "-SLASH-", p, fixed = TRUE)
      results = EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                  hits = genes,
                                                  bg = background,
                                                  reps = bootstrap_reps,
                                                  annotLevel = annotation_Level,
                                                  genelistSpecies=genes_Species,
                                                  sctSpecies=ctd_Species,
                                                  no_cores = number_cores)
      saveRDS(results, paste0(results_dir,"/",p, ".rds"))
      
    })
  
  # Combine results into a single dataframe
  if (MergeResults) {
    results_final <- merge_results(results_dir = results_dir,
                                   list_name_column = list_name_column)
    saveRDS(results_final,paste0("results_",stringr::str_replace_all(Sys.time(),":","-"),".rds"))
    return(results_final)
  } else {
    results_final <- merge_results(results_dir = results_dir,
                                   list_name_column = list_name_column)
    return(results_final)
  }
}
ctd <- readRDS("data/TabulaMuris_h.rds")

gene_data <- read.delim("data/phenotype_to_genes.txt", header = FALSE, skip =1)
colnames(gene_data) = c("HPID", "Phenotype", "EntrezID", "Gene", "Additional", "Source", "LinkID")

all_results <- gen_results(ctd,
                           gene_data,
                           list_names = unique(gene_data$Phenotype),
                           background_genes = unique(gene_data$Gene),
                           list_name_column = "Phenotype",
                           gene_column = "Gene",
                           results_dir = "results_test_level1_110722",
                           overwrite_past_analysis = FALSE,
                           reps = 100000,
                           annotLevel = 1,
                           genelistSpecies = "human",
                           sctSpecies = "human",
                           no_cores = 8,
                           MergeResults = TRUE)
