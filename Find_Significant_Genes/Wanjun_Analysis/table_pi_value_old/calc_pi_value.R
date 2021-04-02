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

comp_df = data.frame(comparison = unlist(lapply(strsplit(
  as.vector(na.omit(unlist(lapply(strsplit(
    list.dirs(), "/"), function(x){x[2]})))), "_"), 
  function(x){x[1]})), dir = list.dirs()[-1])

rank_gene = function(comp_list = c("C1", "C5", "C7")){
  read_comparison = function(x){
    setwd(as.character(subset(comp_df, comparison == x)$dir))
    csv_file = list.files(pattern = "csv")
    data = read.csv(csv_file)
    data = mutate(data, pi_value = abs(logFC) * (-log10(P.Value)))
    data = mutate(data, comparison = x)
    setwd(wd)
    return(data)
  }
  comparison_mat = vector()
  for(i in comp_list){
    comparison_mat = rbind.data.frame(comparison_mat, read_comparison(i))
  }
  comparison_mat = arrange(comparison_mat,
                           desc(pi_value), desc(logFC), desc(P.Value))
  comparison_mat = select(comparison_mat, comparison, ENSEMBL:pi_value)
  return(comparison_mat)
}

rank_gene_interactive = function(rank_df, 
                                 threshold = 200, 
                                 lower_limit = 0.5){
  rank_df = rank_df[1:threshold,]
  aggregate_pi = function(ens, dff){
    pi_value_list = subset(dff, ENSEMBL == ens)$pi_value
    pi_value_list[pi_value_list < lower_limit] = NA
    return(prod(pi_value_list, na.rm = TRUE))
  }
  get_id = function(ens, dff){
    return(unique(subset(dff, ENSEMBL == ens)$ENTREZID))
  }
  get_symbol = function(ens, dff){
    return(unique(subset(dff, ENSEMBL == ens)$SYMBOL))
  }
  interactive_rank_df = with(rank_df, {
    data.frame(ENSEMBL = unique(as.character(ENSEMBL)), 
               ENTREZID = sapply(unique(ENSEMBL), get_id, rank_df),
               SYMBOL = sapply(unique(ENSEMBL), get_symbol, rank_df),
               overall_pi = sapply(unique(ENSEMBL), aggregate_pi, rank_df))
  })
  interactive_rank_df = arrange(interactive_rank_df, desc(overall_pi))
  return(interactive_rank_df)
}

# Assiociated first order comparison 
c157 = rank_gene(c("C1", "C5", "C7"))

# Secondary Comparison
c101819 = rank_gene(c("C10", "C18", "C19"))

# Fasting with CMS pre- post- hemodilution
c91215 = rank_gene(c("C9", "C12", "C15"))

c157_interactive_partial = rank_gene_interactive(c157)
c101819_interactive_partial = rank_gene_interactive(c101819)
c91215_interactive_partial = rank_gene_interactive(c91215)

kable(head(c157,20))
kable(head(c101819,20))
kable(head(c91215,20))

kable(head(c157_interactive_partial,20))
kable(head(c101819_interactive_partial,20))
kable(head(c91215_interactive_partial,20))

# This command takes about 20 minutes to run. 
# Please initiate multi-threading using 'jobRunScript()'
#
# c157_interactive = rank_gene_interactive(c157, threshold = dim(c157)[1])

