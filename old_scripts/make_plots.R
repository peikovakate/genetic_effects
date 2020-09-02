library(tidyverse)
library(ggplot2)
library(rmeta)
library(gridGraphics)
library(cowplot)

ccs = read_tsv("../data/connected_components/microarr_cc.tsv")
cc = ccs %>% filter(phenotype_id == "ILMN_1651438")

p1 = ggplot(cc, aes(x=pos,  y = cell_type, colour = cell_type)) + 
  geom_point(aes(fill=cell_type), size=4, colour="black", pch=21) +
  theme_bw()

p2 = ggplot(cc, aes(x=pos, fill = cell_type)) + 
  geom_dotplot(method="histodot", stackgroups = T, binwidth = 1.4, dotsize = 500, alpha=0.8, stackdir="down") +
  theme_bw() +
  theme(panel.grid=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

plt = plot_grid(p1+theme(legend.position = "none"), 
                p2+theme(legend.position = "none"), 
                nrow=2, align = "v", rel_heights=c(3,1))

legend <- get_legend(p1)
plot_grid(plt, legend, rel_widths = c(3, 1))

variants = read_tsv("../data/found_variants_in_sumstat/found_in_sumstatv2_microarr.tsv")

variants = variants %>% mutate(cell_type = sapply(strsplit(cell_type, "/"), "[[", 2))

variants = variants %>% mutate(eqtl = paste(molecular_trait_id, variant, sep="."))
x = variants %>% group_by(eqtl) %>% summarise(count=n_distinct(cell_type))
xx = x %>% filter(count==17)

summary = variants %>% filter(molecular_trait_id == "ILMN_1651278")
ggplot(summary, aes(x=position, y=beta, colour=cell_type)) +
  geom_point(alpha=0.7) +
  xlab("bp") +
  ylab("Effect size") +
  labs(colour="Cell type")

################################

# mapped_factors = read.table("../data/mfactorization/all_variants_from_credible_sets/mapping_sn_spMF_K15_a1800_l1300_Loadings_beta.txt")
# mapped_factors_pval = read.table("../data/mfactorization/all_variants_from_credible_sets/mapping_sn_spMF_K15_a1800_l1300_Loadings_pvalue_BH.txt")

mapped_factors = read.table("../data/rnaseq/effects_with_na/cc_variants_with_na_mapped/mapping_sn_spMF_K39_a1700_l1700_Loadings_beta.txt")
mapped_factors_pval = read.table("../data/rnaseq/effects_with_na/cc_variants_with_na_mapped/mapping_sn_spMF_K39_a1700_l1700_Loadings_pvalue_BH.txt")

gene_metadata <- read_tsv("additional_data/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv")
colnames(mapped_factors) = 1:ncol(mapped_factors)
colnames(mapped_factors_pval) = 1:ncol(mapped_factors_pval)

pairs = strsplit(row.names(mapped_factors), " ")
genes = sapply(pairs, "[[", 1)
variants = sapply(pairs, "[[", 2)
eqtl_ids = paste(variants, genes, sep=".")
rownames(mapped_factors) = eqtl_ids
rownames(mapped_factors_pval) = eqtl_ids

# effects = read_tsv("pipeline_data/effects_microarr.tsv")
effects = read_tsv("../data/rnaseq/rnaseq_effects.tsv")
variants_metadata = read_tsv("../data/found_variants_in_sumstat/found_variants_microarr.tsv")

betas = effects %>% select(eqtl_id, starts_with("beta."))
ses = effects %>% select(eqtl_id, starts_with("se."))

# betas %>% filter(beta.CEDAR.monocyte_CD14 > 2) %>% select(eqtl_id, beta.CEDAR.monocyte_CD14)
betas %>% filter(beta.GENCORD.LCL > 2) %>% select(eqtl_id, beta.TwinsUK.blood )

# gene = "ILMN_2353161"
# variant_name = "chr16_760040_C_T"

variant_name = "chr9_98193600_A_C"
gene = "ENSG00000106789"

eqtl = paste0(variant_name, ".", gene)
gene_name = filter(gene_metadata, gene_id == gene)[1,]$gene_name
variant_id = filter(variants_metadata, variant == variant_name)[1, ]$rsid


sample_beta = betas %>% filter(eqtl_id==eqtl) %>% select(-eqtl_id)
tissues = sapply(strsplit(names(sample_beta), "beta."), "[[", 2)

colnames(sample_beta) = tissues
sample_se = ses %>% filter(eqtl_id==eqtl) %>% select(-eqtl_id)
colnames(sample_se) = tissues

data = bind_rows(sample_beta, sample_se) 
data = t(data)
colnames(data) = c("beta", "se")
data = as_tibble(data)
data = mutate(data, cell_type = tissues, lower=beta-se**2, upper=beta+sqrt(se))

par(mfrow=c(1,1))

# fp <- ggplot(data, aes(x=cell_type, y=beta, ymin=lower, ymax=upper)) +
#   geom_pointrange() + 
#   geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
#   coord_flip() +  # flip coordinates (puts labels on y axis)
#   xlab("Cell type") + ylab("Effect size") +
#   theme_bw()  # use a white background
# 
# print(fp)

metaplot(unlist(sample_beta), unlist(sample_se), labels = colnames(sample_beta), cex.axis = 3, 
         xlab = sprintf("%s effect size on %s gene", variant_id, gene_name), ylab="Cell types and tissues")
p1 <- recordPlot()

factors = as_tibble(mapped_factors[eqtl,])
factor_pvals = as_tibble(mapped_factors_pval[eqtl, ])
factors = pivot_longer(factors, cols = colnames(factors), names_to = "Factors", values_to = "Loadings")
factor_pvals = pivot_longer(factor_pvals, cols = colnames(factor_pvals), names_to = "Factors", values_to = "p_value")

factors = inner_join(factors, factor_pvals) %>% mutate(p_value = if_else(p_value < 0.05, "< 0.05", "> 0.05"))

p2 <- ggplot(factors, aes(x=Factors, y=Loadings, fill=p_value)) +
  geom_col() +
  scale_fill_manual(values=c("#E69F00", "#999999"))+
  labs(fill="P-value")+
  theme_bw()

plot_grid(p1, p2, nrow=2, rel_heights = c(3,1), labels = "AUTO")

effects




