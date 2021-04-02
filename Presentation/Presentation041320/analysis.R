library(ggplot2)
library(dplyr)
library(ggpubr)
library(latex2exp)
library(ggpmisc)

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

cdp = read.csv(file = "CDP_Data_Spreadsheet.csv", fileEncoding = "UTF-8-BOM")

cdp = mutate(cdp, CONDITION = ifelse(CMS == 0, "CONTROL", ifelse(HMD == 1, "CMS-HEMO", "CMS")))
cdp$CMS = as.factor(cdp$CMS)
cdp$HMD = as.factor(cdp$HMD)

comparison1 = list(c("CMS", "CMS-HEMO"))
comparison2 = list(c("CMS-HEMO", "CONTROL"), c("CMS", "CONTROL"))

source("Boxplot_Figures.R")
source("Regression_Figures.R")