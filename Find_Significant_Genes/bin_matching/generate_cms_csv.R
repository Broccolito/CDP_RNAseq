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
  mutate(genes = as.character(genes)) %>%
  mutate(genes = sapply(genes,function(x){gsub(","," & ",x)}))
cms = cms[cms$cms_score >= 6,]

write.csv(cms, file = "CDP_CompMultSign_list.csv", quote = FALSE, row.names = FALSE)
