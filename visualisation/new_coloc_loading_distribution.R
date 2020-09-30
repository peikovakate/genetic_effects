library(tidyverse)
source("utils2.R")

mapped_data = "../data2/gtex/mfactorization/mapping/mapping_sn_spMF_K50_a11050_l11500"
metadata = read_tsv("../data2/additional/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv")

colocs = read_tsv("../data2/gtex/novel_eqtlcat_colocs_in_novel_ldblocks_relative_to_GTExV8.tsv")
betas_tbl = read.table(paste0(mapped_data, "_Loadings_beta_alpha0.05.txt"))

colocs = mutate(colocs, eqtl_id = paste(variant, molecular_trait_id, sep="."))
gwas = "LC-ebi-a-GCST004627"
gwas_colocs = filter(colocs, gwas_id == gwas)

gwas_colocs = gwas_colocs %>% distinct(eqtl_id, .keep_all = T)

x = t(apply(betas_tbl, 1, function(row){
  abs(row)/max(abs(row))
  }))
x = data.frame(x)

betas = loadings_to_tibble(x)
new_colocs = left_join(gwas_colocs, betas, by="eqtl_id")
new_colocs = left_join(new_colocs, metadata[c("gene_name", "gene_id")], by=c("molecular_trait_id"="gene_id"))

loadings = new_colocs %>% select(starts_with("Factor"))
colnames(loadings) = factor_names$name

factors = tidyr::pivot_longer(loadings, cols = colnames(loadings), names_to = "Factors", values_to = "Loadings")
# factors = filter(factors, Loadings < 0.3)
factors = filter(factors, Loadings != 0) 
ggplot(factors, aes(Factors, Loadings)) +
  geom_jitter(height=0.02, alpha=0.8) +     
  geom_violin() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))



