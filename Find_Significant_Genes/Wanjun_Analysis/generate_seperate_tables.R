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

setwd(get_directory())

d = read.csv(file = "comparisons_master_table.csv")
d = dplyr::select(d, ENSEMBL:gene_type,pi_C1_C5_C7,pi_C10_C18_C19,pi_C9_C12_C15,
                  chr:end_pos,logFC_C1:sig_pi_C19)


fasting_rest = arrange(d, desc(pi_C9_C12_C15))[1:300,] %>%
  dplyr::select(ENSEMBL:SYMBOL,pi_C9_C12_C15,sig_pi_C9,sig_pi_C12,sig_pi_C15,chr,start_pos,end_pos)
exercise_second_order = arrange(d, desc(pi_C10_C18_C19))[1:300,] %>%
  dplyr::select(ENSEMBL:SYMBOL,pi_C10_C18_C19,sig_pi_C10,sig_pi_C18,sig_pi_C19,chr,start_pos,end_pos)
exercise_paired = arrange(d, desc(pi_C1_C5_C7))[1:300,] %>%
  dplyr::select(ENSEMBL:SYMBOL,pi_C1_C5_C7,sig_pi_C1,sig_pi_C5,sig_pi_C7,chr,start_pos,end_pos)

write.csv(fasting_rest, file = "fasting_rest.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(exercise_second_order, file = "exercise_second_order.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(exercise_paired, file = "exercise_paired.csv", 
          row.names = FALSE, quote = FALSE)

