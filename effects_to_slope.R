library(tidyverse)
library(lazyeval)
library(rlang)

output_dir = "../data/mfactorization/20k/"
effects_file = "../thesis_experiments/pipeline_data/effects_microarr.tsv"
gene_names = "../thesis_experiments/additional_data/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv"
one_per_gene = F

dir.create(output_dir)
effects <- read_tsv(effects_file)
gene_names_tbl = read_tsv(gene_names) 
gene_names_tbl = mutate(gene_names_tbl, molecular_trait_id = phenotype_id)

effects = left_join(effects, gene_names_tbl[c("molecular_trait_id", "gene_name")], by="molecular_trait_id")

pvalues <- select(effects, starts_with('pvalue.')) %>% rename_all(function(x){sub("pvalue.", "", x)})
min_cols = max.col(-pvalues)

effects = mutate(effects, min.pval = as.data.frame(pvalues)[cbind(seq_along(min_cols), min_cols)])
nrow(effects)

if(one_per_gene){
  effects = effects %>% 
    group_by(gene_name) %>% 
    top_n(1, min.pval) %>% 
    ungroup()  
}

nrow(effects)

# beta_sums = effects %>% 
#   select(starts_with('beta.')) %>% 
#   rowSums()
# 
# effects = effects %>% 
#   mutate(beta_sum = beta_sums) %>% 
#   arrange(desc(beta_sum)) %>% 
#   top_n(3000)

nrow(effects)
betas <- select(effects, starts_with('beta.')) %>% rename_all(function(x){sub("beta.", "", x)})
pvalues <- select(effects, starts_with('pvalue.')) %>% rename_all(function(x){sub("pvalue.", "", x)})
ses <- select(effects, starts_with('se.')) %>% rename_all(function(x){sub("se.", "", x)})

write_tsv(effects, paste0(output_dir, "effects.tsv"))
write_tsv(betas,  paste0(output_dir, "betas.txt"))
write_tsv(pvalues, paste0(output_dir, "pvalues.txt"))
write_tsv(ses,  paste0(output_dir, "ses.txt"))

write_lines(colnames(betas), paste0(output_dir, "tissues.txt"))

select(effects, variant, molecular_trait_id, starts_with('beta.')) %>%
  rename_all(function(x){sub("beta.", "", x)}) %>%
  rename(SNP = variant, Gene = molecular_trait_id) %>%
  write_tsv(paste0(output_dir, "slope.txt"))

select(effects, variant, molecular_trait_id, starts_with('se.')) %>%
  rename_all(function(x){sub("se.", "", x)}) %>%
  rename(SNP = variant, Gene = molecular_trait_id) %>%
  write_tsv(paste0(output_dir, "se.txt"))

