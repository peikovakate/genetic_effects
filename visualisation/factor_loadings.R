library(tidyverse)
source("utils2.R")

loadings_file = "../data2/gtex/mfactorization/mapping/sn_spMF_K50_a11050_l11500_Loadings_beta_alpha0.05_corrected.txt"
loadings = read.delim(loadings_file)

loadings = loadings_to_tibble(loadings)

n_eqtls = nrow(loadings)

loadings = loadings %>% select(starts_with("Factor"))
colnames(loadings) = factor_names$name
factors = tidyr::pivot_longer(loadings, cols = colnames(loadings), names_to = "Factors", values_to = "Loadings")
fractions = factors %>% group_by(Factors) %>% summarise(count = sum(Loadings != 0), Fraction = count / n_eqtls)


ggplot(fractions, aes(Factors, Fraction)) +
  geom_col() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), text = element_text(size=16))
