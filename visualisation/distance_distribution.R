library(tidyverse)
library(data.table)
library(stringr)
library(ggplot2)
source("utils2.R")

loadings_file = "../data2/gtex/mfactorization/mapping/sn_spMF_K30_a1900_l11100_Loadings_beta_alpha0.05_corrected.txt"
loadings = read.delim(loadings_file)
loadings = loadings_to_tibble(loadings)

factor_names = tibble(factor = paste0("Factor", c(1:16)),
                      name = c("Universal", "iPSC", "Skin", "Blood",
                               "LCL", "T-cell", "Monocyte & Macrophage", "Fat",
                               "T-cell immune response", "Brain", "Macrophage", "Monocyte naive",
                               "BLUEPRINT T-cell", "ROSMAP Brain", "Muscle", "Neutrophil"))
loadings = loadings %>% rename_at(vars(starts_with("Factor")), ~ factor_names$name)

metadata = read_tsv("../data2/additional/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv")
metadata = distinct(metadata, gene_id, phenotype_pos, chromosome)

n_distinct(metadata$gene_id)

factors = tidyr::pivot_longer(loadings, cols = factor_names$name, names_to = "Factors", values_to = "Loadings")
factors = factors %>% filter(Loadings != 0)
factors = factors %>% group_by(eqtl_id) %>% filter(abs(Loadings) == max(abs(Loadings)))
# factors %>% filter()
# factors %>% group_by(eqtl_id) %>% summarise(Factors)
loadings 
factors$position = sapply(factors$variant, str_extract, "(?<=_)(.[0-9]+)(?=_)")
factors$position = as.double(factors$position)

factors = factors %>% left_join(metadata, by=c("gene"="gene_id"))

factors = factors %>% mutate(dist = abs(phenotype_pos-position))
factors = factors %>% mutate(Factors = if_else(Factors == "Universal", "Universal", "Other"))
factors = factors %>% filter(dist<50000)
ggplot(factors, aes(dist, fill=Factors)) +
  geom_density(alpha=0.4)


