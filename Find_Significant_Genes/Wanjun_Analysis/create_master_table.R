rm(list = ls())

library(dplyr)
library(knitr)
if(!require("biomaRt")){
  if(!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("biomaRt")
}
library(biomaRt)

get_directory = function(){
  args <- commandArgs(trailingOnly = FALSE)
  file <- "--file="
  rstudio <- "RStudio"
  match <- grep(rstudio, args)
  if(length(match) > 0){
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }else{
    match <- grep(file, args)
    if (length(match) > 0) {
      return(dirname(normalizePath(sub(file, "", args[match]))))
    }else{
      return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
    }
  }
}

wd = get_directory()
dwd = paste0(wd, "/raw_data")
simwd = paste0(wd, "/simulation_result")

setwd(simwd)
critical_c1 = read.csv("critical_c1.csv")
critical_c5 = read.csv("critical_c5.csv")
critical_c7 = read.csv("critical_c7.csv")
critical_c9 = read.csv("critical_c9.csv")
critical_c10 = read.csv("critical_c10.csv")
critical_c12 = read.csv("critical_c12.csv")
critical_c15 = read.csv("critical_c15.csv")
critical_c18 = read.csv("critical_c18.csv")
critical_c19 = read.csv("critical_c19.csv")

setwd(dwd)
fl = list.files(pattern = "csv", recursive = TRUE)
fl = data.frame(comp = sapply(fl, function(x){unlist(strsplit(x, "_"))[1]}), loc = fl)
fl$comp = as.character(fl$comp)
fl$loc = as.character(fl$loc)

for(i in 1:dim(fl)[1]){
  f = read.csv(fl[i,2])
  f = mutate(f, pi_value = abs(logFC) * (-log10(P.Value)))
  if(exists(paste0("critical_c",i))){
    f = full_join(f,get(paste0("critical_c",i)), by = "ENSEMBL")
  }
  names(f)[names(f)=="logFC"] = paste0("logFC_C",i)
  names(f)[names(f)=="AveExpr"] = paste0("AveExpr_C",i)
  names(f)[names(f)=="P.Value"] = paste0("P_C",i)
  names(f)[names(f)=="pi_value"] = paste0("pi_C",i)
  names(f)[names(f)=="critical_Val"] = paste0("crit_C",i)
  f = dplyr::select(f, ENSEMBL:gene_type,
             starts_with("logFC"), starts_with("AveExpr"),
             starts_with("P"), starts_with("pi"),
             starts_with("crit"))
  assign(paste0("comparison_C",i), f)
}

comparison = comparison_C1
vl = ls(pattern = "comparison_")[-1]
for(v in vl){
  comparison = full_join(comparison, get(v), by = "ENSEMBL")
}
comparison = dplyr::select(comparison, ENSEMBL,ENTREZID,SYMBOL,gene_type,
                    starts_with("logFC"), starts_with("AveExpr"),
                    starts_with("P"), starts_with("pi"),
                    starts_with("crit"))

comparison = comparison %>%
  mutate(sig_pi_C1 = ifelse(pi_C1>=crit_C1,pi_C1,0)) %>%
  mutate(sig_pi_C5 = ifelse(pi_C5>=crit_C5,pi_C5,0)) %>%
  mutate(sig_pi_C7 = ifelse(pi_C7>=crit_C7,pi_C7,0)) %>%
  mutate(sig_pi_C9 = ifelse(pi_C9>=crit_C9,pi_C9,0)) %>%
  mutate(sig_pi_C10 = ifelse(pi_C10>=crit_C10,pi_C10,0)) %>%
  mutate(sig_pi_C12 = ifelse(pi_C12>=crit_C12,pi_C12,0)) %>%
  mutate(sig_pi_C15 = ifelse(pi_C15>=crit_C15,pi_C15,0)) %>%
  mutate(sig_pi_C18 = ifelse(pi_C18>=crit_C18,pi_C18,0)) %>%
  mutate(sig_pi_C19 = ifelse(pi_C19>=crit_C19,pi_C19,0)) %>%
  mutate(pi_C1_C5_C7 = (ifelse(sig_pi_C1>=1,sig_pi_C1,1) *
                            ifelse(sig_pi_C5>=1,sig_pi_C5,1) *
                            ifelse(sig_pi_C7>=1,sig_pi_C7,1))) %>%
  mutate(pi_C10_C18_C19 = (ifelse(sig_pi_C10>=1,sig_pi_C10,1) *
                               ifelse(sig_pi_C18>=1,sig_pi_C18,1) *
                               ifelse(sig_pi_C19>=1,sig_pi_C19,1))) %>%
  mutate(pi_C9_C12_C15 = (ifelse(sig_pi_C9>=1,sig_pi_C9,1) *
                              ifelse(sig_pi_C12>=1,sig_pi_C12,1) *
                              ifelse(sig_pi_C15>=1,sig_pi_C15,1)))

# listEnsemblArchives()
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl",
                  host = "grch37.ensembl.org",
                  path="/biomart/martservice", verbose = TRUE)
searchFilters(mart = ensembl, pattern = "ensembl")
get_chr_loc = function(){
  tryCatch({
    chromosome_info <<- getBM(attributes=c("ensembl_gene_id",
                                           "chromosome_name",
                                           "start_position",
                                           "end_position"), filter="ensembl_gene_id",
                              values = comparison$ENSEMBL, mart = ensembl)
  }, error = function(e){
    cat("\n Server failed to connect. Re-connecting... \n")
    get_chr_loc()
  })
}
get_chr_loc()
colnames(chromosome_info) = c("ENSEMBL","chr","start_pos","end_pos")
comparison = full_join(comparison, chromosome_info, by = "ENSEMBL")

setwd(wd)
names(comparison)
write.csv(comparison, file = "comparisons_master_table.csv", quote = FALSE, row.names = FALSE)
