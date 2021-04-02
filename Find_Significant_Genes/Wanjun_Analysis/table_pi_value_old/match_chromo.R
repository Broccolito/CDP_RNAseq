rm(list = ls())

library(dplyr)
library(knitr)

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
setwd(wd)

g = read.csv(file = "gene_list_comparison.csv")
chr = read.csv(file = "chromosome_info.csv")

names(chr)[names(chr) == "ensembl_gene_id"] = "ENSEMBL"

g_chr = left_join(g, chr, by = "ENSEMBL")
g_chr = select(g_chr,ENSEMBL:gene_type, chromosome_name:length, 
               pi_value_C1:pi_value_C9_C12_C15)
write.csv(g_chr, file = "gene_position_list_comparison.csv", row.names = FALSE, quote = FALSE)