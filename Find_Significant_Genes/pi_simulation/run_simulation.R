rm(list = ls())

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

source("simulate_critical_pi.R")
set.seed(123) # Seed used for the simulation
# Iteration number of the simulation is 1000

cat("\n Simulation starting time: \n")
print(Sys.time())

#--------------------------------------------

critical_c1 = get_critical_pi(condition1="CMS",condition2="CMS",
                              treatment1="PREHEM_POSTEX",treatment2="PREHEM_PREEX",
                              entire_compare1 ="",entire_compare2 = "",sig_level=0.9)
write.csv(critical_c1, file = "critical_c1.csv", quote = FALSE, row.names = FALSE)

critical_c5 = get_critical_pi(condition1="CMS",condition2="CMS",
                              treatment1="POSTHEM_POSTEX",treatment2="POSTHEM_PREEX",
                              entire_compare1 ="",entire_compare2 = "",sig_level=0.9)
write.csv(critical_c5, file = "critical_c5.csv", quote = FALSE, row.names = FALSE)

critical_c7 = get_critical_pi(condition1="CON",condition2="CON",
                              treatment1="PREHEM_POSTEX",treatment2="PREHEM_PREEX",
                              entire_compare1 ="",entire_compare2 = "",sig_level=0.9)
write.csv(critical_c7, file = "critical_c7.csv", quote = FALSE, row.names = FALSE)

#--------------------------------------------

critical_c10 = get_critical_pi(condition1="",condition2="",
                               treatment1="",treatment2="",
                               entire_compare1='C5_CMS_POSTHEM_POSTEXvsCMS_POSTHEM_PREEX',
                               entire_compare2 = 'C1_CMS_PREHEM_POSTEXvsCMS_PREHEM_PREEX', sig_level=0.9)
write.csv(critical_c10, file = "critical_c10.csv", quote = FALSE, row.names = FALSE)

critical_c18 = get_critical_pi(condition1="",condition2="",
                               treatment1="",treatment2="",
                               entire_compare1='C1_CMS_PREHEM_POSTEXvsCMS_PREHEM_PREEX',
                               entire_compare2 = 'C7_CON_PREHEM_POSTEXvsCON_PREHEM_PREEX', sig_level=0.9)
write.csv(critical_c18, file = "critical_c18.csv", quote = FALSE, row.names = FALSE)

critical_c19 = get_critical_pi(condition1="",condition2="",
                               treatment1="",treatment2="",
                               entire_compare1='C5_CMS_POSTHEM_POSTEXvsCMS_POSTHEM_PREEX',
                               entire_compare2 = 'C7_CON_PREHEM_POSTEXvsCON_PREHEM_PREEX', sig_level=0.9)
write.csv(critical_c19, file = "critical_c19.csv", quote = FALSE, row.names = FALSE)

#--------------------------------------------

critical_c10_0.8 = get_critical_pi(condition1="",condition2="",
                               treatment1="",treatment2="",
                               entire_compare1='C5_CMS_POSTHEM_POSTEXvsCMS_POSTHEM_PREEX',
                               entire_compare2 = 'C1_CMS_PREHEM_POSTEXvsCMS_PREHEM_PREEX', sig_level=0.8)
write.csv(critical_c10_0.8, file = "critical_c10_0.8.csv", quote = FALSE, row.names = FALSE)

critical_c18_0.8 = get_critical_pi(condition1="",condition2="",
                               treatment1="",treatment2="",
                               entire_compare1='C1_CMS_PREHEM_POSTEXvsCMS_PREHEM_PREEX',
                               entire_compare2 = 'C7_CON_PREHEM_POSTEXvsCON_PREHEM_PREEX', sig_level=0.8)
write.csv(critical_c18_0.8, file = "critical_c18_0.8.csv", quote = FALSE, row.names = FALSE)

critical_c19_0.8 = get_critical_pi(condition1="",condition2="",
                               treatment1="",treatment2="",
                               entire_compare1='C5_CMS_POSTHEM_POSTEXvsCMS_POSTHEM_PREEX',
                               entire_compare2 = 'C7_CON_PREHEM_POSTEXvsCON_PREHEM_PREEX', sig_level=0.8)
write.csv(critical_c19_0.8, file = "critical_c19_0.8.csv", quote = FALSE, row.names = FALSE)

#--------------------------------------------

critical_c9 = get_critical_pi(condition1="CMS",condition2="CMS",
                              treatment1="POSTHEM_FASTING",treatment2="PREHEM_FASTING",
                              entire_compare1 ="",entire_compare2 = "",sig_level=0.9)
write.csv(critical_c9, file = "critical_c9.csv", quote = FALSE, row.names = FALSE)

critical_c12 = get_critical_pi(condition1="CMS",condition2="CON",
                              treatment1="PREHEM_FASTING",treatment2="PREHEM_FASTING",
                              entire_compare1 ="",entire_compare2 = "",sig_level=0.9)
write.csv(critical_c12, file = "critical_c12.csv", quote = FALSE, row.names = FALSE)

critical_c15 = get_critical_pi(condition1="CMS",condition2="CON",
                              treatment1="POSTHEM_FASTING",treatment2="PREHEM_FASTING",
                              entire_compare1 ="",entire_compare2 = "",sig_level=0.9)
write.csv(critical_c15, file = "critical_c15.csv", quote = FALSE, row.names = FALSE)

#--------------------------------------------

cat("\n Simulation ending time: \n")
print(Sys.time())