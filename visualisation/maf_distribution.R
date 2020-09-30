library(tidyverse)
library(reshape2)

betas = read.table("../data2/gtex/mfactorization/mapping_sn_spMF_K50_a11000_l11000_Loadings_beta_alpha0.05.txt")
colnames(betas) = 1: 
eqtls = strsplit(rownames(betas), split=" ")
variants = sapply(eqtls, "[[", 2)
betas$variant = variants

factors = reshape2::melt(betas)
factors = as_tibble(factors)

factors %>% filter(value != 0, variable == 11) %>% select(variant) %>%  
  write_tsv("../data2/gtex/factor11_variants.txt")

factors %>% filter(value != 0, variable == 13) %>% select(variant) %>%  
  write_tsv("../data2/gtex/factor13_variants.txt")

afr11 = read_tsv("../data2/gtex/visualisation/afr.f11.tsv") %>% mutate(factor = 11, population="AFR")
afr13 = read_tsv("../data2/gtex/visualisation/afr.f13.tsv") %>% mutate(factor = 13, population="AFR")
eur11 = read_tsv("../data2/gtex/visualisation/eur.f11.tsv") %>% mutate(factor = 11, population="EUR")
eur13 = read_tsv("../data2/gtex/visualisation/eur.f13.tsv") %>% mutate(factor = 13, population="EUR")

distribution = bind_rows(afr13, eur13, eur11, afr11)
ggplot(distribution, aes(population, maf)) +
  ggplot2::geom_jitter(height = 0, alpha=0.1) + 
  ggplot2::geom_violin() + 
  theme_bw() +
  facet_grid(. ~ factor)
