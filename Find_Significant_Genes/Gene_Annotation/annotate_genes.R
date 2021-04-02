library(dplyr)
library(GeneBook)
library(xlsx)

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

exercise_paired = read.csv(file = "exercise_paired.csv")
exercise_second_order = read.csv(file = "exercise_second_order.csv")
fasting_rest = read.csv(file = "fasting_rest.csv")


from_gene_to_annotation = function(gene_list){
  annotation_mat = vector()
  for(gene in gene_list){
    if(dim(GeneCard_Symbol_Details(gene))[1] == 0){
      NULL_result = data.frame(gene = NA,
                               type = NA,
                               description = NA,
                               summary_entrez = NA,
                               summary_genecard = NA,
                               summary_uniport = NA,
                               summary_Tocris = NA,
                               summary_CIViC = NA)
      annotation_mat = rbind.data.frame(annotation_mat,
                                        NULL_result)
    }else{
      annotation_mat = rbind.data.frame(annotation_mat,
                                        GeneCard_Symbol_Details(gene))
    }
  }
  names(annotation_mat) = names(GeneCard_Symbol_Details(NA))
  annotation_mat = select(annotation_mat, type:summary_CIViC)
  return(annotation_mat)
}

exercise_paired_annotate = select(cbind.data.frame(exercise_paired, 
                                                   from_gene_to_annotation(exercise_paired$SYMBOL)),
                                  ENSEMBL:SYMBOL, type, starts_with("pi_value"), description:summary_CIViC)

exercise_second_order_annotate = select(cbind.data.frame(exercise_second_order, 
                                                         from_gene_to_annotation(exercise_second_order$SYMBOL)),
                                        ENSEMBL:SYMBOL, type, starts_with("pi_value"), description:summary_CIViC)

fasting_rest_annotate = select(cbind.data.frame(fasting_rest, 
                                                from_gene_to_annotation(exercise_paired$SYMBOL)),
                               ENSEMBL:SYMBOL, type, starts_with("pi_value"), description:summary_CIViC)

write.xlsx(exercise_paired_annotate, file = "exercise_paired_annotate.xlsx", 
           sheetName = "Gene_Annotation", row.names = FALSE)
write.xlsx(exercise_second_order_annotate, file = "exercise_second_order_annotate.xlsx",
           sheetName = "Gene_Annotation", row.names = FALSE)
write.xlsx(fasting_rest_annotate, file = "fasting_rest_annotate.xlsx", 
           sheetName = "Gene_Annotation", row.names = FALSE)