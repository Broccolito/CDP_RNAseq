rm(list = ls())

library(dplyr)

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

load("fasting_rest_match.RData")
load("exercise_paired_match.RData")
load("exercise_second_order_match.RData")

cms = readLines("CDP_CompMultSign_list.txt")
cms = sapply(cms, function(x){
  s = gsub(" ","", unlist(strsplit(x, "  ")))
  s = s[nchar(s)>0]
  s = unlist(s)
  names(s) = NULL
  return(s)
})
names(cms) = NULL

for(i in 1:length(cms)){
  if(length(cms[[i]]) == 4){
    cms[[i]] = c(cms[[i]], "NONAME")
  }
}

cms = as.data.frame(matrix(unlist(cms), ncol = 5, byrow = TRUE))
colnames(cms) = c("chr", "pos", "score", "cms_score", "genes")
cms = select(cms, chr:pos, cms_score:genes) %>%
  mutate(pos = as.numeric(as.character(pos))) %>%
  mutate(cms_score = as.numeric(as.character(cms_score))) %>%
  mutate(genes = as.character(genes)) %>%
  mutate(genes = sapply(genes,function(x){gsub(","," & ",x)}))
cms = cms[cms$cms_score >= 6,] # Cutoff


d = read.csv("comparisons_master_table.csv")
d = mutate(d, chr = paste0("chr",chr))
fasting_rest = arrange(d, desc(pi_C1_C5_C7))[1:300,]
exercise_second_order = arrange(d, desc(pi_C10_C18_C19))[1:300,]
exercise_paired = arrange(d, desc(pi_C9_C12_C15))[1:300,]

exercise_paired_match = arrange(exercise_paired_match, seq_pos)
exercise_second_order_match = arrange(exercise_second_order_match, seq_pos)
fasting_rest_match = arrange(fasting_rest_match, seq_pos)

generate_match_mat = function(ddf = fasting_rest, 
                              ddf_match = fasting_rest_match){
  seq_df = select(ddf[ddf_match$seq_pos,],ENSEMBL,SYMBOL,
                  pi_C1_C5_C7,pi_C10_C18_C19,pi_C9_C12_C15,
                  chr:end_pos)
  cms_df = cms[ddf_match$cms_pos,] %>%
    mutate(variant_loc = paste(chr,pos, sep = ":")) %>%
    mutate(related_genes = genes) %>%
    select(variant_loc, cms_score, related_genes)
  res = cbind.data.frame(seq_df,cms_df)
  row.names(res) = NULL
  return(res)
}

fasting_rest_cms = generate_match_mat(fasting_rest, 
                                      fasting_rest_match)
exercise_paired_cms = generate_match_mat(exercise_paired, 
                                         exercise_paired_match)
exercise_second_order_cms = generate_match_mat(exercise_second_order, 
                                               exercise_second_order_match)

write.csv(fasting_rest_cms, file = "fasting_rest_cms.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(exercise_paired_cms, file = "exercise_paired_cms.csv", 
          row.names = FALSE, quote = FALSE)
write.csv(exercise_second_order_cms, file = "exercise_second_order_cms.csv", 
          row.names = FALSE, quote = FALSE)

