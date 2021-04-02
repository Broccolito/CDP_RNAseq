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
dwd = paste0(wd, "/raw_data")
simwd = paste0(wd, "/simulation")
reswd = paste0(wd, "/result")

setwd(dwd)

for(f in list.files(all.files = TRUE,pattern="csv",recursive = TRUE)){
  d = read.csv(f) %>%
    select(ENSEMBL,logFC,P.Value) %>%
    mutate(pi_value = abs(logFC)*(-log10(P.Value)))
  assign(paste0(unlist(strsplit(f,"_"))[1],"_data"),d)
}

setwd(simwd)
for(f in list.files(all.files = TRUE,pattern="csv",recursive = TRUE)){
  d = read.csv(f)
  assign(paste0(toupper(unlist(strsplit(f,"_"))[1]),"_sim"),d)
}

fl = cbind(paste0("C",1:19,"_data"),paste0("C",1:19,"_sim"),paste0("C",1:19))
for(i in 1:19){
  assign(fl[i,3], left_join(get(fl[i,1]),get(fl[i,2]),by = "ENSEMBL"))
} 

rm(list = ls(pattern = "_data"))
rm(list = ls(pattern = "_sim"))
rm(d,fl)
gc()

generate_ipa_input = function(d,fn){
  d = arrange(d,desc(pi_value))
  d_95 = d %>%
    mutate(pi_value = ifelse(pi_value>=sig95,pi_value,0)) %>%
    arrange(desc(pi_value)) %>%
    filter(pi_value!=0) %>%
    select(ENSEMBL:P.Value)
  d_90 = d %>%
    mutate(pi_value = ifelse(pi_value>=sig90,pi_value,0)) %>%
    arrange(desc(pi_value)) %>%
    filter(pi_value!=0) %>%
    select(ENSEMBL:P.Value)
  d_85 = d %>%
    mutate(pi_value = ifelse(pi_value>=sig85,pi_value,0)) %>%
    arrange(desc(pi_value)) %>%
    filter(pi_value!=0) %>%
    select(ENSEMBL:P.Value)
  d_80 = d %>%
    mutate(pi_value = ifelse(pi_value>=sig80,pi_value,0)) %>%
    arrange(desc(pi_value)) %>%
    filter(pi_value!=0) %>%
    select(ENSEMBL:P.Value)
  write.csv(d,file = paste0(fn,"_overall.csv"),quote = FALSE,row.names = FALSE)
  write.csv(d_95,file = paste0(fn,"_95.csv"),quote = FALSE,row.names = FALSE)
  write.csv(d_90,file = paste0(fn,"_90.csv"),quote = FALSE,row.names = FALSE)
  write.csv(d_85,file = paste0(fn,"_85.csv"),quote = FALSE,row.names = FALSE)
  write.csv(d_80,file = paste0(fn,"_80.csv"),quote = FALSE,row.names = FALSE)
}

setwd(reswd)
for(i in paste0("C",1:19)){
  generate_ipa_input(get(i),i)
}

