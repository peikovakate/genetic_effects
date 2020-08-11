library(tidyverse)
library(lazyeval)
library(rlang)

effects_file = "../data/rnaseq/rnaseq_effects.tsv"
output_library = "../data/rnaseq/unique_genes/"
genes_metadata = "../thesis_experiments/additional_data/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv"
one_eqtl_per_gene = T

dir.create(output_library)

# read effects file across all tissues
# eqtl per row, contains beta, se and pval across set of tissues and cell types
effects <- read_tsv(effects_file)
print(sprintf("Reading effects file with %i records", nrow(effects)))

# file with gene names
gene_names_tbl = read_tsv(genes_metadata) 
# gene_names_tbl = mutate(gene_names_tbl, molecular_trait_id = phenotype_id)
gene_names_tbl = mutate(gene_names_tbl, molecular_trait_id = gene_id)

effects = left_join(effects, gene_names_tbl[c("molecular_trait_id", "gene_name")], by="molecular_trait_id")
effects %>% distinct(gene_name) %>% nrow()

effects = effects %>% distinct(.keep_all = T)

if(one_eqtl_per_gene){
  # find cell type with smallest pvalue across eqtls
  pvalues <- select(effects, starts_with('pvalue.')) %>% rename_all(function(x){sub("pvalue.", "", x)})
  # replace na with inf, so max.col does not return index with na 
  min_cols = max.col(-replace(pvalues, is.na(pvalues), Inf))
  effects = mutate(effects, min.pval = as.data.frame(pvalues)[cbind(seq_along(min_cols), min_cols)])
  nrow(effects)
  
  # select one eqtl per gene (eqtl with smallest pval across all cell types)
  effects = effects %>% 
    group_by(gene_name) %>% 
    # not sure that it returns exactyly one row
    # also min.pval soring is correct?
    top_n(1, min.pval) %>% 
    ungroup()
  
  nrow(effects)
}


betas <- select(effects, starts_with('beta.')) %>% rename_all(function(x){sub("beta.", "", x)})
pvalues <- select(effects, starts_with('pvalue.')) %>% rename_all(function(x){sub("pvalue.", "", x)})
ses <- select(effects, starts_with('se.')) %>% rename_all(function(x){sub("se.", "", x)})

print(sprintf("writing files to %s", output_library))

write_tsv(effects, file.path(output_library, "effects.tsv"))
write_tsv(betas,  file.path(output_library, "betas.txt"))
write_tsv(pvalues, file.path(output_library, "pvalues.txt"))
write_tsv(ses,  file.path(output_library, "ses.txt"))
write_lines(colnames(betas), file.path(output_library, "tissues.txt"))

# this is input for matrix factorizaiton: sn_spMF folder
select(effects, variant, molecular_trait_id, starts_with('beta.')) %>%
  rename_all(function(x){sub("beta.", "", x)}) %>%
  rename(SNP = variant, Gene = molecular_trait_id) %>%
  write_tsv(file.path(output_library, "slope.txt"))

select(effects, variant, molecular_trait_id, starts_with('se.')) %>%
  rename_all(function(x){sub("se.", "", x)}) %>%
  rename(SNP = variant, Gene = molecular_trait_id) %>%
  write_tsv(file.path(output_library, "se.txt"))

