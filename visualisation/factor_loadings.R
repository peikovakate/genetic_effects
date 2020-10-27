library(tidyverse)
source("utils2.R")
factor_names = tibble(factor = paste0("Factor", c(1:16)),
                      name = c("Universal", "iPSC", "Skin", "Blood",
                               "LCL", "T-cell", "Monocyte & Macrophage", "Fat",
                               "T-cell immune response", "Brain", "Macrophage", "Monocyte naive",
                               "BLUEPRINT T-cell", "ROSMAP Brain", "Muscle", "Neutrophil" ))

loadings_file = "../data2/gtex/mfactorization/mapping/sn_spMF_K30_a1900_l11100_Loadings_beta_alpha0.05_corrected.txt"
loadings = read.delim(loadings_file)

loadings = loadings_to_tibble(loadings)

n_eqtls = nrow(loadings)

loadings = loadings %>% select(starts_with("Factor"))
colnames(loadings) = factor_names$name


factors = tidyr::pivot_longer(loadings, cols = colnames(loadings), names_to = "Factors", values_to = "Loadings")
fractions = factors %>% group_by(Factors) %>% summarise(count = sum(Loadings != 0), Fraction = count / n_eqtls)


ggplot(fractions, aes(Factors, Fraction)) +
  geom_col() + 
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), text = element_text(size=16))
