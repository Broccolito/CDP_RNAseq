library(edgeR)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggpubr)
library(ggthemes)

counts = read.table("cdp_counts_matrix")
groups = as.factor(c("CMS_POSTHEM_FASTING",
                     "CMS_POSTHEM_FASTING",
                     "CMS_POSTHEM_FASTING",
                     "CMS_POSTHEM_FASTING",
                     "CMS_POSTHEM_FASTING",
                     "CMS_POSTHEM_FASTING",
                     "CMS_POSTHEM_FASTING",
                     "CMS_POSTHEM_FASTING",
                     "CMS_POSTHEM_FASTING",
                     "CMS_POSTHEM_FASTING",
                     "CMS_POSTHEM_POSTEX",
                     # "CMS_POSTHEM_POSTEX",
                     "CMS_POSTHEM_POSTEX",
                     "CMS_POSTHEM_POSTEX",
                     "CMS_POSTHEM_POSTEX",
                     "CMS_POSTHEM_POSTEX",
                     "CMS_POSTHEM_POSTEX",
                     # "CMS_POSTHEM_PREEX",
                     "CMS_POSTHEM_PREEX",
                     "CMS_POSTHEM_PREEX",
                     "CMS_POSTHEM_PREEX",
                     "CMS_POSTHEM_PREEX",
                     "CMS_POSTHEM_PREEX",
                     "CMS_POSTHEM_PREEX",
                     "CMS_PREHEM_FASTING",
                     "CMS_PREHEM_FASTING",
                     "CMS_PREHEM_FASTING",
                     # "CMS_PREHEM_FASTING",
                     "CMS_PREHEM_FASTING",
                     "CMS_PREHEM_FASTING",
                     "CMS_PREHEM_POSTEX",
                     "CMS_PREHEM_POSTEX",
                     "CMS_PREHEM_POSTEX",
                     "CMS_PREHEM_POSTEX",
                     "CMS_PREHEM_POSTEX",
                     "CMS_PREHEM_PREEX",
                     "CMS_PREHEM_PREEX",
                     "CMS_PREHEM_PREEX",
                     "CMS_PREHEM_PREEX",
                     "CMS_PREHEM_PREEX",
                     "CMS_PREHEM_PREEX",
                     "CON_PREHEM_FASTING",
                     "CON_PREHEM_FASTING",
                     "CON_PREHEM_FASTING",
                     "CON_PREHEM_FASTING",
                     "CON_PREHEM_FASTING",
                     "CON_PREHEM_FASTING",
                     "CON_PREHEM_FASTING",
                     "CON_PREHEM_FASTING",
                     "CON_PREHEM_POSTEX",
                     "CON_PREHEM_POSTEX",
                     "CON_PREHEM_POSTEX",
                     "CON_PREHEM_POSTEX",
                     "CON_PREHEM_POSTEX",
                     "CON_PREHEM_POSTEX",
                     # "CON_PREHEM_POSTEX",
                     "CON_PREHEM_POSTEX",
                     "CON_PREHEM_PREEX",
                     "CON_PREHEM_PREEX",
                     "CON_PREHEM_PREEX",
                     "CON_PREHEM_PREEX",
                     "CON_PREHEM_PREEX",
                     "CON_PREHEM_PREEX",
                     "CON_PREHEM_PREEX",
                     "CON_PREHEM_PREEX"))

groups = factor(groups, levels = c(
  "CMS_POSTHEM_FASTING", # 1
  "CMS_POSTHEM_POSTEX",  # 2
  "CMS_POSTHEM_PREEX",   # 3
  "CMS_PREHEM_FASTING",  # 4
  "CMS_PREHEM_POSTEX",   # 5
  "CMS_PREHEM_PREEX",    # 6
  "CON_PREHEM_FASTING",  # 7
  "CON_PREHEM_POSTEX",   # 8
  "CON_PREHEM_PREEX"     # 9
))

counts = counts %>%
  select(-genes.EL8_aligned.genes.results) %>%
  select(-genes.EL6_aligned.genes.results) %>%
  select(-genes.EL42_aligned.genes.results) %>%
  select(-genes.EL61_aligned.genes.results)

d = DGEList(counts = counts, group = groups)

library_size = apply(d$counts, 2, sum)
library_size = tibble(sample = names(library_size),library_size)

keep = rowSums(cpm(d)>100) >= 2
d = d[keep,]
d$samples$lib.size = colSums(d$counts)

library_size = library_size %>%
  mutate(filtered_library_size = colSums(d$counts))

d = calcNormFactors(d)

# plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
# legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

d = estimateCommonDisp(d, verbose=F)
d = estimateTagwiseDisp(d)

get_dge = function(condition1, condition2){
  et = exactTest(d, pair=c(condition1,condition2))
  dge_list = topTags(et,1e6)$table
  dge_list = dge_list %>% 
    mutate(gene_name = rownames(dge_list)) %>%
    select(gene_name, everything()) %>%
    arrange(FDR)
  rownames(dge_list) = NULL
  return(dge_list)
}

# "CMS_POSTHEM_FASTING", # 1
# "CMS_POSTHEM_POSTEX",  # 2
# "CMS_POSTHEM_PREEX",   # 3
# "CMS_PREHEM_FASTING",  # 4
# "CMS_PREHEM_POSTEX",   # 5
# "CMS_PREHEM_PREEX",    # 6
# "CON_PREHEM_FASTING",  # 7
# "CON_PREHEM_POSTEX",   # 8
# "CON_PREHEM_PREEX"     # 9

# CMS PREHEM Exercise Response
cms_prehem = get_dge(5,6) %>%
  mutate(Condition = "CMS Pre-Hemoglobin")
write.csv(cms_prehem, file = "cms_prehem.csv",
          quote = FALSE, row.names = FALSE)

# CMS POSTHEM Exercise Response
cms_posthem = get_dge(2,3) %>%
  mutate(Condition = "CMS Post-Hemoglobin")
write.csv(cms_posthem, file = "cms_posthem.csv",
          quote = FALSE, row.names = FALSE)

# CON PREHEM Exercise Response
con_prehem = get_dge(8,9) %>%
  mutate(Condition = "Control")
write.csv(con_prehem, file = "con_prehem.csv",
          quote = FALSE, row.names = FALSE)


# Fasting
cmspre_con_fasting = get_dge(4,7) %>%
  mutate(Condition = "CMS PreHem Vs. Control")
write.csv(cmspre_con_fasting, file = "cmspre_con_fasting.csv",
          quote = FALSE, row.names = FALSE)

cmspost_con_fasting = get_dge(1,7) %>%
  mutate(Condition = "CMS PostHem Vs. Control")
write.csv(cmspost_con_fasting, file = "cmspost_con_fasting.csv",
          quote = FALSE, row.names = FALSE)

cmspost_cmspre_fasting = get_dge(1,4) %>%
  mutate(Condition = "CMS PostHem Vs. CMS PreHem")
write.csv(cmspost_cmspre_fasting, file = "cmspost_cmspre_fasting.csv",
          quote = FALSE, row.names = FALSE)


volcano_matrix = rbind.data.frame(cms_prehem,
                                  cms_posthem,
                                  con_prehem)

# Make Volcano Plot
ggplot(data = volcano_matrix, aes(x = logFC, y = -log(PValue,base = 10))) + 
  geom_point() + 
  facet_grid(.~Condition) + 
  xlab("Log Fold Change") + 
  ylab("-Log(P)") + 
  xlim(c(-2.6,2.6)) +
  ylim(c(0, 2.5)) + 
  theme_gray() + 
  theme(text = element_text(size = 15)) + 
  ggsave(filename = "CDP Exercise Response Volcano Plot.png",
         device = "png", dpi = 1200,
         width = 11, height = 4)




