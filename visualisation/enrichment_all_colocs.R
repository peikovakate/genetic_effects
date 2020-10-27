library(tidyverse)
source("utils2.R")

metadata = read_tsv("../data2/additional/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv")
colocs = read_tsv("../data2/additional/eqtlCat_all_sign_colocs.tsv")

colocs = mutate(colocs, eqtl_id = paste(variant, molecular_trait_id, sep="."))
colocs = colocs %>% select(gwas_id, eqtl_id) %>% distinct()


n_distinct(colocs$eqtl_id)

# loadings_file = "../data2/gtex/mfactorization/mapping/sn_spMF_K50_a11050_l11500_Loadings_beta_alpha0.05_corrected.txt"
loadings_file = "../data2/gtex/mfactorization/mapping/mapping_sn_spMF_K50_a11050_l11500_Loadings_beta_alpha0.05.txt"
loadings = read.delim(loadings_file)
loadings = loadings_to_tibble(loadings) 

pvalues_file = "../data2/gtex/mfactorization/mapping/mapping_sn_spMF_K50_a11050_l11500_Loadings_pvalue_BH.txt"
pvalues = read.delim(pvalues_file)
pvalues = loadings_to_tibble(pvalues) 
# pvalues = left_join(pvalues, metadata[c("gene_name", "gene_id")], by=c("gene"="gene_id"))

pvalues = pvalues %>% 
  rowwise() %>% 
  mutate_at(vars(starts_with("Factor")), function(x){if (x < 0) 1 else x}) %>% ungroup()

subset_eqtl_id = 
  pvalues %>% 
  pivot_longer(c(-gene, -variant, -eqtl_id), names_to="Factor") %>% 
  group_by(gene) %>% 
  mutate(min_p = min(value)) %>% 
  pivot_wider(names_from=Factor) %>% 
  filter_at(vars(starts_with("Factor")), any_vars(. == min_p)) %>% 
  sample_n(1) %>% 
  ungroup() %>% 
  arrange(gene)
# %>% select(-variant, -gene)

n_distinct(pvalues$gene)
nrow(subset_eqtl_id)

loadings = loadings %>% filter(eqtl_id %in% subset_eqtl_id$eqtl_id)

loading_coloc = left_join(loadings, colocs, by="eqtl_id")
# loading_coloc = left_join(loading_coloc, metadata[c("gene_name", "gene_id")], by=c("gene"="gene_id"))
# loading_coloc = loading_coloc %>% group_by(gene_name) %>% sample_n(1)

factors = tidyr::pivot_longer(loading_coloc, cols = colnames(loadings)[startsWith(colnames(loadings), "Factor")], 
                              names_to = "Factors", values_to = "Loadings")

factors = factors %>% mutate(abs_loading = abs(Loadings))

max_factor = factors %>% 
  group_by(eqtl_id) %>% 
  mutate( max_abs_loading = max(abs_loading)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(Loadings = if (abs_loading == max_abs_loading) Loadings else 0) %>% 
  ungroup()

max_factor

# factors %>% group_by(gwas_id) %>% summarise(c = n_distinct(eqtl_id))

fisher_test = function(factors){
  counts = factors %>% 
    filter(!is.na(gwas_id)) %>% 
    group_by(gwas_id, Factors) %>% 
    summarise(load = sum(Loadings != 0), dont_load = sum(Loadings == 0)) %>% 
    ungroup()
  
  counts_na = factors %>% 
    filter(is.na(gwas_id)) %>% 
    group_by(Factors) %>% 
    summarise(no_gwas_load = sum(Loadings != 0), no_gwas_dont_load = sum(Loadings == 0)) %>% 
    ungroup()
  
  test = left_join(counts, counts_na, by="Factors")
  
  run_test = function(a, b, c, d){
    data = rbind(c(a,b), c(c,d))
    t = fisher.test(data)
    return(list(p=t$p.value, or=t$estimate))
  }
  # log10(t$p.value) * (if (t$estimate < 0) -1 else 1)
  test = test %>%  
    rowwise() %>% 
    mutate(p = run_test(no_gwas_dont_load, no_gwas_load, dont_load, load)$p, 
           or = run_test(no_gwas_dont_load, no_gwas_load, dont_load, load)$or)
  
  test = test %>% 
    mutate(signed_log10p = log10(p) * (if (or < 0) -1 else 1))
  
  return(test)
  
}
test = fisher_test(factors)
test = fisher_test(max_factor)
test = left_join(test, factor_names, by=c("Factors"="factor"))
significant = test %>% filter(p < 0.05)

ggplot(test, aes(name, gwas_id, fill=signed_log10p)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", high="white") +
  geom_point(data=significant) + 
  xlab("Factors") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))

# gwas_colocs = filttate(p = run_er(loading_coloc, gwas_id == gwas) %>% distinct(eqtl_id, .keep_all = T)

# gwas_colocs %>% summarise(Factor)
# gwas_colocs %>% distinct(eqtl_id) %>% nrow()




