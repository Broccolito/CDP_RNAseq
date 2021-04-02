library("dplyr")

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

set.seed(123)
df = read.csv(file='cdp.csv',colClasses=c('character'))

#get all patients ids like "EL1_S1_L004_R1_001", etc.
df$File_Name<-gsub("\\..*","",df$File_Name)

#data frame contains patients ids and treatments conditions. etc
patients<- select(df,-c(File_Name.1,Iteration,libType,Tissue,Date_Collected))
# data with DE data for all patients
rawdata=read.csv(file="20191108_simonson_rnaseq_data_exploration_20200218092250_filtered_pc_genes_tmm_normed_cpms.csv")
#change the raw data to log2 expression level 
rawdata[rawdata == 0] <- NA
lg_rawdata= apply(rawdata[,-1], MARGIN = c(1,2), FUN = log2)
lg_rawdata<-as.data.frame(lg_rawdata)
lg_rawdata$X<-rawdata$X
#select patients according to conditions
#for example select_patients(patients,condition="CMS",treatment="PREHEM_POSTEX")
#return the file name of the pacients get selected
select_pat<-function(patients_menu,condition,treatment){
  df<-patients_menu %>%  filter(Condition == condition & 
                                  (Treatment == treatment ))
  return(df$File_Name)
}

select_patients<-function(patients_menu,condition,treatment,entire_comparison){
  if(nchar(entire_comparison)>1){
    subtext<-sub(".*[1-9]_",'',entire_comparison)
    group2<-sub('.*vs', '', entire_comparison)
    group1<-gsub('\\s*vs.*$', '', subtext)
    a<-unlist(strsplit(group1, '_'))
    b<-unlist(strsplit(group2, '_'))
    condition1<-a[1]
    condition2<-b[1]
    treatment1<-paste(a[2],a[3],sep="_")
    treatment2<-paste(b[2],b[3],sep="_")
    patient1<-select_pat(patients_menu,condition1,treatment1)
    patient2<-select_pat(patients_menu,condition2,treatment2)
    r<-c(patient1,patient2)
    return(r)
  }
  else{return(select_pat(patients_menu,condition,treatment))
    
  }
}

#simulate data under null hypothesis, which is there is no differential expression genes.
list_of_patients<-names(lg_rawdata)
simulate_null_dist<-function(condition1,condition2,treatment1,treatment2,entire_compare1,entire_compare2){
  #get patients file name
  g1_patients<-select_patients(patients,condition1,treatment1,entire_compare1)
  g2_patients<-select_patients(patients,condition2,treatment2,entire_compare2)
  #get patients that in the list of patients in lg_rawdata(since some do not appear in lg_rawdata)
  g1_patients<-g1_patients[g1_patients %in% list_of_patients]
  g2_patients<-g2_patients[g2_patients %in% list_of_patients]
  #get data for the patients of two treatment groups
  g<-lg_rawdata[,c(g1_patients,g2_patients)]
  num_m1<-length(g1_patients)
  num_m2<-length(g2_patients)
  num_m<-num_m1+num_m2
  m1<-sample(1:num_m,num_m1,TRUE)
  m2<-sample(1:num_m,num_m2,TRUE)
  g1<-g[,m1]
  g2<-g[,m2]
  simulated_data <- list(g1,g2)
  return (simulated_data)
}
#using data under null hypothesis, simulate the p value and obtain empirical distirbution
simulate_pi<-function(condition1,condition2,treatment1,treatment2,entire_compare1,entire_compare2){
  pi_values<-data.frame(gene=rawdata$X)
  # Using a large it number can increase the accuracy of the simulation
  for(i in seq(1,1000)){
    g<-simulate_null_dist(condition1,condition2,treatment1,treatment2,entire_compare1,entire_compare2)
    g1<-as.data.frame(g[1])
    g2<-as.data.frame(g[2])
    g1_mean<-rowMeans(g1,na.rm = TRUE)
    g2_mean<-rowMeans(g2,na.rm = TRUE)
    log_ratio<-as.vector(g1_mean-g2_mean)
    logfc<-abs(log_ratio)
    #calculate pvalue
    std<-as.vector(apply(g1,1, sd, na.rm = TRUE))
    t<-(g1_mean-g2_mean)/(std/sqrt(ncol(g1)))
    calc_p <- function(t) {2*pt(-abs(t),df=ncol(g1)-1) }
    pval<-sapply(t, calc_p)
    #calculate pi values
    pi_val<-as.vector(logfc* (-log10(pval)))
    #add results to list
    pi_values<-cbind(pi_values,pi_val)
  }
  return (pi_values)
  
}
#calculate the critical pi value of each gene
get_critical_pi<-function(condition1,condition2,
                          treatment1,treatment2,
                          entire_compare1,entire_compare2,
                          sig_level){
  l<-simulate_pi(condition1,condition2,treatment1,treatment2,entire_compare1,entire_compare2)
  critical_pi<-c()
  for(row in 1:nrow(l)){
    pi <- l[row,-1 ]
    sig_pi<-quantile(pi, sig_level,na.rm=TRUE)
    critical_pi<-c(critical_pi,sig_pi)
  }
  critical_val<-data.frame(ENSEMBL=rawdata$X,critical_Val=critical_pi)
  return(critical_val)
}
#example C18 (double comparison C1 VS C7)
# result<-get_critical_pi(condition1="",condition2="",
#                         treatment1="",treatment2="",
#                         entire_compare1='C1_CMS_PREHEM_POSTEXvsCMS_PREHEM_PREEX',
#                         entire_compare2 = 'C7_CON_PREHEM_POSTEXvsCON_PREHEM_PREEX', sig_level=0.8)
# head(result)




#example for ones like 

#C1_CMS_PREHEM_POSTEXvsCMS_PREHEM_PREEX
# microbenchmark::microbenchmark({
  # result<-get_critical_pi(condition1="CMS",condition2="CMS",
  #                         treatment1="PREHEM_POSTEX",treatment2="PREHEM_PREEX",
  #                         entire_compare1 ="",entire_compare2 = "",sig_level=0.9)
#   head(result)
# }, times = 1L)





