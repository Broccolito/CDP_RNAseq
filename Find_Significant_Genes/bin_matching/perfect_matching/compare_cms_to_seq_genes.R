rm(list = ls())

library(dplyr)
library(readxl)
library(rstudioapi)

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
  mutate(genes = as.character(genes))
cms = cms[cms$cms_score >= 6,] # Cutoff

d = read.csv("comparisons_master_table.csv")
d = mutate(d, chr = paste0("chr",chr))
fasting_rest = arrange(d, desc(pi_C1_C5_C7))[1:300,]
exercise_second_order = arrange(d, desc(pi_C10_C18_C19))[1:300,]
exercise_paired = arrange(d, desc(pi_C9_C12_C15))[1:300,]

find_match = function(seq_db = fasting_rest){
  matching_list = vector()
  for(i in 1:dim(cms)[1]){
    variant_from_cms = cms[i,]
    for(j in 1:dim(seq_db)[1]){
      variant_from_seq = seq_db[j,]
      variant_from_cms$chr = as.character(variant_from_cms$chr)
      variant_from_seq$chr = as.character(variant_from_seq$chr)
      if(variant_from_cms$chr == variant_from_seq$chr){
        if(variant_from_cms$pos >= variant_from_seq$start_pos &
           variant_from_cms$pos <= variant_from_seq$end_pos){
          matching_list = rbind(matching_list, c(j,i))
          cat("\nFind a match!")
        }
      }
    }
    if(i%%100==0){cat(paste0("\n", i, " out of ", dim(cms)[1], " done..."))}
  }
  return(data.frame(seq_pos = matching_list[,1], cms_pos = matching_list[,2]))
}

jobRunScript(path = "run1.R", importEnv = TRUE, workingDir = wd, name = "Fasting Rest")
jobRunScript(path = "run2.R", importEnv = TRUE, workingDir = wd, name = "Exercise Second Order")
jobRunScript(path = "run3.R", importEnv = TRUE, workingDir = wd, name = "Exercise Paired")