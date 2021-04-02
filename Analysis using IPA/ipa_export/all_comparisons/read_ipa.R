library(readxl)

rm(list = ls())

read_ipa = function(filename = "c1_90_all.xls"){
  
  suppressWarnings({
    comp = read_xls(filename)
  })
  
  comp = cbind.data.frame(0,comp)
  n = 1
  for(s in comp[,2]){
    n = n + 1
    if(grepl("Canonical Pathways for My Projects", 
             s, fixed = TRUE)){
      comp[n,1] = 1
    }
    if(grepl("Upstream Regulators", 
             s, fixed = TRUE)){
      comp[n,1] = 1
    }
    if(grepl("Causal Networks for My Projects", 
             s, fixed = TRUE)){
      comp[n,1] = 1
    }
    if(grepl("Diseases and Bio Functions for My Projects", 
             s, fixed = TRUE)){
      comp[n,1] = 1
    }
    if(grepl("Tox Functions for My Projects", 
             s, fixed = TRUE)){
      comp[n,1] = 1
    }
    if(grepl("Regulator Effects for My Projects", 
             s, fixed = TRUE)){
      comp[n,1] = 1
    }
    if(grepl("Networks for My Projects", 
             s, fixed = TRUE)){
      comp[n,1] = 1
    }
    if(grepl("My Lists for My Projects", 
             s, fixed = TRUE)){
      comp[n,1] = 1
    }
    if(grepl("Tox Lists for My Projects", 
             s, fixed = TRUE)){
      comp[n,1] = 1
    }
    if(grepl("My Pathways for My Projects", 
             s, fixed = TRUE)){
      comp[n,1] = 1
    }
    if(grepl("Analysis Ready Molecules for My Projects", 
             s, fixed = TRUE)){
      comp[n,1] = 1
    }
  }
  
  sep = which(comp[,1]==1)
  df1 = comp[1:sep[1],]
  df2 = comp[(sep[1]):(sep[2]-2),]
  df3 = comp[(sep[2]):(sep[3]-2),]
  df4 = comp[(sep[3]):(sep[4]-2),]
  df5 = comp[(sep[4]):(sep[5]-2),]
  df6 = comp[(sep[5]):(sep[6]-2),]
  df7 = comp[(sep[6]):(sep[7]-2),]
  df8 = comp[(sep[7]):(sep[8]-2),]
  df9 = comp[(sep[8]):(sep[9]-2),]
  df10 = comp[(sep[9]):(sep[10]-2),]
  df11 = comp[(sep[10]):(sep[11]-2),]
  df12 = comp[(sep[11]):dim(comp)[1],]
  
  process_df = function(dff){
    dff = dff[,which(!is.na(dff[1,]))]
    names(dff) = dff[1,]
    dff = dff[-1,]
    dff = dff[,-1]
    return(dff)
  }
  
  df2 = process_df(df2)
  df3 = process_df(df3)
  df4 = process_df(df4)
  df5 = process_df(df5)
  df6 = process_df(df6)
  df7 = process_df(df7)
  df8 = process_df(df8)
  df9 = process_df(df9)
  df10 = process_df(df10)
  df11 = process_df(df11)
  df12 = process_df(df12)
  
  list_from_ipa = list(canonical_pathway = df2,
                       upstream_regulator = df3,
                       causal_network = df4,
                       diseases_bio_function = df5,
                       tox_function = df6,
                       regulator_effect = df7,
                       my_network = df8,
                       my_list = df9,
                       tox_list = df10,
                       my_pathway = df11,
                       analysis_ready_molecule = df12)
  
  return(list_from_ipa)
  
}

all_comps = list(
  c1 = read_ipa("c1_90_all.xls"),
  c2 = read_ipa("c2_90_all.xls"),
  c3 = read_ipa("c3_90_all.xls"),
  c4 = read_ipa("c4_90_all.xls"),
  c5 = read_ipa("c5_90_all.xls"),
  c6 = read_ipa("c6_90_all.xls"),
  c7 = read_ipa("c7_90_all.xls"),
  c8 = read_ipa("c8_90_all.xls"),
  c9 = read_ipa("c9_90_all.xls"),
  c10 = read_ipa("c10_90_all.xls"),
  c11 = read_ipa("c11_90_all.xls"),
  c12 = read_ipa("c12_90_all.xls"),
  c13 = read_ipa("c13_90_all.xls"),
  c14 = read_ipa("c14_90_all.xls"),
  c15 = read_ipa("c15_90_all.xls"),
  c16 = read_ipa("c16_90_all.xls"),
  c17 = read_ipa("c17_90_all.xls"),
  c18 = read_ipa("c18_90_all.xls"),
  c19 = read_ipa("c19_90_all.xls")
)

save(all_comps,file = "all_comparisons_from_ipa.RData")
