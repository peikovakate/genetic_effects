colours = read_tsv("../data2/gtex/qtl_groups_colours.tsv")
metadata = read_tsv("../data2/gtex/qtl_groups_data.tsv")
metadata = left_join(metadata, colours)

load("../data2/gtex/mfactorization/results/sn_spMF_K50_a11050_l11500_Run1.RData")
pheatmap::pheatmap(FactorM, fontsize=8)

factor_names = tibble(factor = paste0("Factor", c(1:12)),
                      name = c("Universal", "T-cell", "LCL",  "Blood",
                               "iPSC", "Neutrophil", "Fat", "Muscle", 
                               "T-cell immune response","Brain", "Monocyte naive","Monocyte & Macrophage"))
FactorM
factror_loadings = tibble(reshape2::melt(FactorM, value.name="loading"))
factror_loadings = dplyr::rename(factror_loadings, qtl_group=Var1, factor=Var2)
factror_loadings = dplyr::left_join(factror_loadings, metadata[c("qtl_group", "colour", "group")])


df = dplyr::filter(factror_loadings, factor == 1)
ggplot2::ggplot(df, aes(qtl_group, loading, fill=group)) +
  ggplot2::geom_col() +
  scale_fill_manual(values=colours$colour) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
