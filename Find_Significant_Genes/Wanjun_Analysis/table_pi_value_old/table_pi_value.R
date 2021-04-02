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

rank_df = rank_gene("C1")
for(comp in paste0("C", seq(2,19))){
  rank_df = full_join(rank_df, select(rank_gene(comp),ENSEMBL,pi_value),
                      by = "ENSEMBL", suffix = c("", paste0("_", comp)))
  rank_df = full_join(rank_df, select(rank_gene(comp),ENSEMBL,logFC),
                      by = "ENSEMBL", suffix = c("", paste0("_", comp)))
  rank_df = full_join(rank_df, select(rank_gene(comp),ENSEMBL,P.Value),
                      by = "ENSEMBL", suffix = c("", paste0("_", comp)))
}
names(rank_df)[names(rank_df) == "pi_value"] = "pi_value_C1"
names(rank_df)[names(rank_df) == "logFC"] = "logFC_C1"
names(rank_df)[names(rank_df) == "P.Value"] = "P.Value_C1"

rank_df = select(rank_df, ENSEMBL, ENTREZID,
                 SYMBOL, gene_type, starts_with("pi_value"),
                 starts_with("logFC"),
                 starts_with("P.Value"))

get_overall_pi = function(group = c("C1", "C5", "C7"),
                          threshold = 1){
  overall_pi = apply(select(select(rank_df, starts_with("pi_value")), ends_with(group)),
                     MARGIN = 1,
                     FUN = function(x){
                       x[x<= threshold] = NA
                       return(prod(x, na.rm = TRUE))
                     })
  return(overall_pi)
}

rank_df = mutate(rank_df, 
                 pi_value_C1_C5_C7 = get_overall_pi(group = c("C1", "C5", "C7")),
                 pi_value_C10_C18_C19 = get_overall_pi(group = c("C10", "C18", "C19")),
                 pi_value_C9_C12_C15 = get_overall_pi(group = c("C9", "C12", "C15")))
write.csv(rank_df, file = "gene_list_comparison.csv", quote = FALSE, row.names = FALSE)

pilot_sample_size = 500
# Exercise paired
C1_C5_C7 = head(arrange(select(rank_df, ENSEMBL:gene_type,
                               logFC_C1, logFC_C5, logFC_C7,
                               P.Value_C1,P.Value_C5,P.Value_C7,
                               pi_value_C1, pi_value_C5, pi_value_C7,
                               pi_value_C1_C5_C7), desc(pi_value_C1_C5_C7)),
                pilot_sample_size)
write.csv(C1_C5_C7, file = "exercise_paired.csv", quote = FALSE, row.names = FALSE)

# Exercise second order
C10_C18_C19 = head(arrange(select(rank_df, ENSEMBL:gene_type,
                                  logFC_C10, logFC_C18, logFC_C19,
                                  P.Value_C10,P.Value_C18,P.Value_C19,
                                  pi_value_C10, pi_value_C18, pi_value_C19,
                                  pi_value_C10_C18_C19), desc(pi_value_C10_C18_C19)),
                   pilot_sample_size)
write.csv(C10_C18_C19, file = "exercise_second_order.csv", quote = FALSE, row.names = FALSE)

# Fasting Rest
C9_C12_C15 = head(arrange(select(rank_df, ENSEMBL:gene_type,
                                 logFC_C9, logFC_C12, logFC_C15,
                                 P.Value_C9,P.Value_C12,P.Value_C15,
                                 pi_value_C9, pi_value_C12, pi_value_C15,
                                 pi_value_C9_C12_C15), desc(pi_value_C9_C12_C15)),
                  pilot_sample_size)
write.csv(C9_C12_C15, file = "fasting_rest.csv", quote = FALSE, row.names = FALSE)