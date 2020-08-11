library(tidyverse)
library(ggplot2)

gwas = "bmi"
folder = "../data/rnaseq/effects_with_na/cc_variants_with_na_mapped/heritability_maxscale/"
# folder = "../data/rnaseq/cc_average/heritability/"
herit = read_tsv(file.path(folder, sprintf("%s.results", gwas)))
nrow(herit)

herit = arrange(herit, Category)

categories = strsplit(herit$Category, split = "L2_0")
categories = unlist(lapply(categories, "[[", 1))

herit = mutate(herit, Category=categories)
herit = mutate(herit, p_value=if_else(Enrichment_p <= 0.05, "< 0.05", "> 0.05"))
significant = herit %>% filter(Enrichment_p <= 0.05) 
# herit = herit %>% filter(Enrichment_p <= 0.05)

ggplot(herit, aes(x=Category, y = Enrichment, fill=p_value))+
  ggplot2::geom_col() +
  ggplot2::geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2) + 
  # ggplot2::ggtitle(paste("GWAS", toupper(gwas))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values=c("#E69F00", "#999999"))+
  labs(fill="P-value")
  

# ggsave(file.path(folder, paste0(gwas, ".png")))
