library(tidyverse)
source("utils2.R")

mapped_data = "../data2/gtex/mfactorization/mapping/mapping_sn_spMF_K50_a11050_l11500"
metadata = read_tsv("../data2/additional/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv")

colocs = read_tsv("../data2/gtex/novel_eqtlcat_colocs_in_novel_ldblocks_relative_to_GTExV8.tsv")
betas_tbl = read.table(paste0(mapped_data, "_Loadings_beta_alpha0.05.txt"))
effects = read_tsv("../data2/gtex/lead_effects_na.tsv")

colocs = mutate(colocs, eqtl_id = paste(variant, molecular_trait_id, sep="."))
gwas = "LC-ebi-a-GCST004627"
gwas_colocs = filter(colocs, gwas_id == gwas)
data = filter(gwas_colocs, qtl_subset == "BLUEPRINT_PE.T-cell_ge")

betas = loadings_to_tibble(betas_tbl)

new_colocs = left_join(data, betas)
new_colocs = left_join(new_colocs, effects)
new_colocs = left_join(new_colocs, metadata[c("gene_name", "gene_id")], by=c("molecular_trait_id"="gene_id"))
matrs = effects_to_matricies(new_colocs)

matrs$beta
loadings = new_colocs %>% select(starts_with("Factor"))
colnames(loadings) = factor_names$name

plot_forrest_loadings = function(beta, se, variant_id, gene_name, loadings){
  rmeta::metaplot(unlist(beta), unlist(se), labels = names(beta),
           xlab = sprintf("%s effect size on %s gene", variant_id, gene_name), ylab="qtl groups")
  
  p1 <- recordPlot()
  
  factors = tidyr::pivot_longer(loadings, cols = colnames(loadings), names_to = "Factors", values_to = "Loadings")
  p2 <-ggplot(factors, aes(x=factor(Factors, levels = colnames(loadings)), y=Loadings)) + 
    geom_col() + 
    xlab("Factors") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size=16))
  
  plt = cowplot::plot_grid(p1, p2, nrow=2, rel_heights = c(3,1))
  return(plt)
}

# ind = 1
# dir.create(file.path("figures", gwas))
for(ind in c(1:nrow(new_colocs))){
  plt = plot_forrest_loadings(matrs$beta[ind,], matrs$se[ind,], new_colocs$variant[ind], new_colocs$gene_name[ind], loadings[ind,])
  ggsave(filename=sprintf("figures/%s/%s.png", gwas, new_colocs$eqtl_id[ind]), plot=plt, width=8, height=20)
}
  
