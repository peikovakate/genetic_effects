library(tidyverse)
library(reshape2)
library(optparse)
library(ggplot2)

model_path = "../data/rnaseq/effects_with_na/cc_variants_with_na_mapped/mapping_sn_spMF_K39_a1700_l1700_Loadings"
variants_metadata_file = "../data/rnaseq/effects_with_na/cc_variants_with_na_mapped/rnaseq_variants_metadata.tsv"
factors_tbl_output =  "../data/rnaseq/effects_with_na/cc_variants_with_na_mapped/factor_binary_one_per_eqtl.tsv"

betas = read.table(paste0(model_path, "_beta.txt"))
N_factors = ncol(betas)
pvalues = read.table(paste0(model_path, "_pvalue_BH.txt"))
corrected_betas = read.table(paste0(model_path, "_beta_alpha0.05_corrected.txt"), sep="\t")

nrow(corrected_betas)
nrow(betas)
nrow(pvalues)

colnames(betas) = sapply(strsplit(colnames(betas), split = "V"), "[[", 2)
colnames(corrected_betas) = sapply(strsplit(colnames(corrected_betas), split = "V"), "[[", 2)
colnames(pvalues) = sapply(strsplit(colnames(pvalues), split = "X"), "[[", 2)

genes = sapply(strsplit(rownames(betas), split = " "), "[[", 1)
variants = sapply(strsplit(rownames(betas), split = " "), "[[", 2)
eqtl_ids = paste(variants, genes, sep=".")

######################### binary ###########################
stand_betas = betas
stand_betas[,] = 1
stand_betas$eqtl_id = eqtl_ids
stand_betas_factors = reshape2::melt(stand_betas, id="eqtl_id", variable.name = "factor", value.name="stand_beta")
############################################################

################# stand beta or max scale ##################
# stand_betas = as.data.frame(t(scale(t(betas), center = F))) # divide by sd
#stand_betas = as.data.frame(t(apply(betas, 1, function(x){(abs(x)-min(abs(x)))/max(abs(x))}))) # min-max scale
stand_betas = as.data.frame(t(apply(betas, 1, function(x){x/max(abs(x))}))) # max scale

stand_betas$eqtl_id = eqtl_ids
stand_betas_factors = reshape2::melt(stand_betas, id="eqtl_id", variable.name = "factor", value.name="stand_beta")
############################################################

########################### pip ############################
# cc_file = "../data/connected_components/microarr_cc.tsv"
# ccs = read_tsv(cc_file)
# ccs = ccs %>% mutate(eqtl_id = paste(variant_id, phenotype_id, sep="."))
# 
# pip_df = 
#   ccs %>% 
#   distinct(eqtl_id, cell_type, pip) %>% 
#   group_by(eqtl_id) %>% 
#   summarise(max_pip = max(pip)) 
# 
# stand_betas_factors = expand_grid(pip_df, factor = as.character(seq(1:N_factors)))
# stand_betas_factors = rename(stand_betas_factors, stand_beta = max_pip)
###########################################################

#################### product ##############################
# stand_betas = as.data.frame(t(apply(betas, 1, function(x){abs(x)/max(abs(x))})))
# stand_betas$eqtl_id = eqtl_ids
# stand_betas_factors = reshape2::melt(stand_betas, id="eqtl_id", variable.name = "factor", value.name="stand_beta")
# 
# cc_file = "../data/connected_components/microarr_cc.tsv"
# ccs = read_tsv(cc_file)
# ccs = ccs %>% mutate(eqtl_id = paste(variant_id, phenotype_id, sep="."))
# 
# pip_df = 
#   ccs %>% 
#   distinct(eqtl_id, cell_type, pip) %>% 
#   group_by(eqtl_id) %>% 
#   summarise(max_pip = max(pip)) 
# 
# stand_betas_factors = 
#   stand_betas_factors %>% left_join(pip_df, by="eqtl_id") %>% mutate(stand_beta = stand_beta*max_pip)

#########################################################

corrected_betas$count = apply(corrected_betas, 1, function(row){
  return(sum(row != 0))
})

pvalues$eqtl_id = eqtl_ids
betas$eqtl_id = eqtl_ids
corrected_betas$eqtl_id = eqtl_ids


eqtl_factors = reshape2::melt(pvalues, id="eqtl_id", variable.name = "factor", value.name="pvalue")
betas_factors = reshape2::melt(betas, id="eqtl_id", variable.name = "factor", value.name="beta")
corrected_factors = reshape2::melt(corrected_betas, id.vars=c("eqtl_id", "count"), variable.name = "factor", value.name="corrected_beta")

corrected_betas %>% filter(count == 1) %>% nrow()

eqtl_factors = left_join(eqtl_factors, betas_factors, by=c("eqtl_id","factor")) %>% 
  left_join(stand_betas_factors, by=c("eqtl_id","factor")) %>% 
  left_join(corrected_factors, by=c("eqtl_id", "factor"))

nrow(eqtl_factors)
# eqtl_factors = filter(eqtl_factors, pvalue <= 0.05, pvalue > 0)

eqtl_factors = filter(eqtl_factors, corrected_beta != 0, count==1)
nrow(eqtl_factors)

eqtl_factors = mutate(eqtl_factors, abs_beta = abs(beta))

# pick one factor per eqtl
# eqtl_factors = eqtl_factors %>% 
#   group_by(factor) %>% 
#   top_n(15000, abs_beta) %>% 
#   ungroup()

ggplot(eqtl_factors, aes(x=stand_beta, fill=factor)) + geom_density(alpha=0.3)

variants_metadata = read_tsv(variants_metadata_file)
eqtl_factors = left_join(eqtl_factors, variants_metadata, by="eqtl_id")
eqtl_factors %>% nrow()
eqtl_factors %>% distinct(eqtl_id) %>% nrow()

write_tsv(eqtl_factors, factors_tbl_output)

