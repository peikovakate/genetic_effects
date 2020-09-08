`%>%` <- magrittr::`%>%`
library(ggplot2)
library(rmeta)
library(gridGraphics)
library(cowplot)

eqtls_file = "../data2/additional/novel_sign_eqtlcat_colocs.tsv"
variants_file = "../data/rnaseq/variants_across_tissues_rnaseq_no_effect_to_na_fract0.tsv"
loadings_file = "../data/rnaseq/effects_with_na/cc_variants_with_na_mapped/mapping_sn_spMF_K39_a1700_l1700_Loadings_beta.txt"

variants = readr::read_tsv(variants_file)
variants = dplyr::rename(variants, qtlGroup = cell_type)
eqtl_ids = paste(variants, genes, sep=".")

factor_loadings = read.table(loadings_file)
genes = sapply(strsplit(rownames(factor_loadings), split = " "), "[[", 1)
variants = sapply(strsplit(rownames(factor_loadings), split = " "), "[[", 2)
eqtl_ids = paste(variants, genes, sep=".")
rownames(factor_loadings) = eqtl_ids
colnames(factor_loadings) = as.character(1:ncol(factor_loadings))


###########################
gwas = "IBD-ebi-a-GCST004131"
dir.create(sprintf("figures/%s", gwas))

eqtls = readr::read_tsv(eqtls_file)
eqtls = dplyr::filter(eqtls, gwas_id == gwas)

print("Records")
eqtls %>% nrow()

print("Unique pairs")
eqtls %>% dplyr::distinct(variant, molecular_trait_id) %>% nrow()

eqtls = dplyr::distinct(eqtls, variant, molecular_trait_id, .keep_all = T)
eqtls = dplyr::mutate(eqtls, eqtl_id = paste(variant, molecular_trait_id, sep="."))

effects = sumstat_to_effects(variants)
column_names = c("eqtl_id", "variant", "molecular_trait_id")

betas = dplyr::inner_join(eqtls[column_names], effects$beta, by="eqtl_id")
ses =  dplyr::inner_join(eqtls[column_names], effects$beta, by="eqtl_id")
loadings = factor_loadings[betas$eqtl_id,]

to_vector = function(data){
  data %>% dplyr::select(!column_names) %>% as.vector()
}

plot_eqtl = function(beta, se, loadings){
  eqtl_id = beta["eqtl_id"]
  variant_id = beta["variant"]
  gene_name = se["molecular_trait_id"]
  beta = to_vector(beta)
  se = to_vector(se)
  
  metaplot(unlist(beta), unlist(se), labels = colnames(bs), cex.axis = 3, 
           xlab = sprintf("%s effect size on %s gene", variant_id, gene_name), ylab="qtl groups")
  p1 <- recordPlot()
  
  factors = tidyr::pivot_longer(loadings, cols = colnames(loadings), names_to = "Factors", values_to = "Loadings")
  p2 <-ggplot(factors, aes(x=factor(Factors, levels = colnames(loadings)), y=Loadings)) + 
    geom_col() + 
    xlab("Factors") + 
    theme_bw()
  
  plt = plot_grid(p1, p2, nrow=2, rel_heights = c(4,1))
  
  ggsave(filename=sprintf("figures/%s/%s.png", gwas, eqtl_id), plot=plt, width=9, height = 12)
}

for (i in 1:nrow(betas)) {
  plot_eqtl(betas[i, ], ses[i,], loadings[i, ])
}

mapply(plot_eqtl, betas, ses, loadings)

plot_eqtl(betas[2,], ses[2,], loadings[2,])



