library(dplyr)

beta_file = "../data/mfactorization/credible_set/mapping_sn_spMF_K15_a1800_l1300_Loadings_beta_alpha0.05.txt"
variants_metadata_file = "../data/found_variants_in_sumstat/variants_metadata.tsv"
betas = read.table(beta_file)

nrow(betas)

genes = sapply(strsplit(rownames(betas), split = " "), "[[", 1)
variants = sapply(strsplit(rownames(betas), split = " "), "[[", 2)
eqtl_ids = paste(variants, genes, sep=".")

membership_matr = matrix(0, dim(betas)[1], dim(betas)[2])
membership_matr[betas > 0] = 1

eqtl_factors <- apply(membership_matr, 1, function(x)(which(x==1)))
names(eqtl_factors) = eqtl_ids
eqtl_factors = mapply(function(l, n){
  tibble(factor=l, eqtl_id=n)
}, eqtl_factors, names(eqtl_factors), SIMPLIFY = F)
eqtl_factors = bind_rows(eqtl_factors)
nrow(eqtl_factors)

variants_metadata = read_tsv(variants_metadata_file)
eqtl_factors = left_join(eqtl_factors, variants_metadata, by="eqtl_id")
eqtl_factors %>% distinct(eqtl_id) %>% nrow()
write_tsv(eqtl_factors, "../data/found_variants_in_sumstat/found_variants_factors.tsv")


