rm(list = ls())
graphics.off()

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

seq_data = read.csv(file = "RNAseq_Raw.csv")
seq_data_new = as.data.frame(t(seq_data[,-1]))
colnames(seq_data_new) = seq_data$X
seq_data_new = cbind.data.frame(data.frame(id = colnames(seq_data)[-1]), 
                                seq_data_new)

meta = read.csv(file = "RNAseq_Meta.csv", fileEncoding="UTF-8-BOM")
suppressWarnings({
  seq_data_new = na.omit(full_join(y = seq_data_new, x = meta, by = "id"))
  names(seq_data_new)[1] = "RNAid"
})

# Add more physiology data

write.csv(seq_data_new, file = "RNAseq_Processed.csv", quote = FALSE,
            row.names = FALSE)

# Generate RData file for shinyapp
cdp = seq_data_new
save(cdp, file = "shiny_data.RData")


