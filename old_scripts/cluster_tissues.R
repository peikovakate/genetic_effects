library(tidyverse)
library(mashr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
source("visualisation_utils.R")

# effects_file = "../data/microarray/effects_microarr.tsv"
effects_file <- "../data/rnaseq/rnaseq_effects.tsv"

effects <- read_tsv(effects_file)
nrow(effects)

effects_matrix <-
  select(effects, starts_with('beta.')) %>% rename_all(function(x) {
    sub("beta.", "", x)
  })

errors_matrix <- select(effects, starts_with('se.')) %>% rename_all(function(x){sub("se.", "", x)})
pval_matrix <- select(effects, starts_with('pvalue.')) %>% rename_all(function(x){sub("pvalue.", "", x)})

effects_matrix <- as.matrix(effects_matrix)
# effects_matrix <- effects_matrix*sign(effects_matrix[max.col(abs(effects_matrix))])

# multiply by the sign of the strongest effect across tissues
# effects_matrix <- t(apply(effects_matrix, 1, function(x){
#   abs_x <- abs(x)
#   x*sign(x[which.max(abs_x)])
# }))

errors_matrix <- as.matrix(errors_matrix)
pval_matrix <- as.matrix(pval_matrix)

alpha_value = 1

data   = mash_set_data(effects_matrix, errors_matrix, alpha=alpha_value)
m.1by1 = mash_1by1(data, alpha=alpha_value)
strong = get_significant_results(m.1by1)
U.c    = cov_canonical(data)
U.pca  = cov_pca(data, 5, strong)
U.ed   = cov_ed(data, U.pca)

# m = mash(data, c(U.c, U.ed))
m = mash(data, U.ed)
# save(m, file = "../data/mash/microarr_mash_1_ed_only.R")
# load("../data/mash/rnaseq_mash_1_ed_only.R")

sharing = get_pairwise_sharing(m)
coords = MASS::isoMDS(1-sharing)

cols.cor <- cor(effects_matrix, method = m)
coords = MASS::isoMDS(1-cols.cor)

coords = tibble(x = coords$points[,1], y=coords$points[,2], qtlGroup=rownames(coords$points))
ggplot2::ggplot(coords, aes(x, y)) +
  ggplot2::geom_point() +
  geom_label_repel(aes(label = qtlGroup), size=3) +
  theme_classic() 


  # geom_label_repel(aes(label = qtlGroup))
  # ggplot2::geom_text(aes(label=qtlGroup),hjust=0, vjust=0)

# pheatmap(get_pairwise_sharing(m.ed))
hmap(
  as.data.frame(sharing),
  fontSize = 6.5,
  colorPalette = 8,
  dendroLineSize = 0.3,
  scaleName="eQTL sharing"
  # main="MashR pairwise tissue sharing"
)

pheatmap(sharing, cutree_cols = 6, 
         filename = "final_figures/microarr/mash_a1_raw_effects_only_ed_cuttree.png",
         colorRampPalette(brewer.pal(n = 7, name ="OrRd"))(100), width=7, height=6)

ggsave("final_figures/rnaseq/mash_a1_raw_effects_only_ed.png", scale=1, width = 6, height = 5)

# pm <- get_pm(m)
# clusters <- kmeans(pm, 10, nstart = 10, iter.max = 30)
# clusters$size
# for (i in 1:10) {
#   c <- pm[clusters$cluster == i,]
#   tissue_names <- names(c[1, ])
#   c <- aperm(c)
#   # c <- aperm(sign(c) * c)
#   matplot(
#     c,
#     type = "l",
#     lwd = 0.1,
#     lty = 1,
#     main = paste("cluster", i),
#     axes = F,
#     las = 2,
#     ylab = "Posterior means",
#     xlab = "cell types"
#   )
#   axis(2)
#   # axis(side=1,at=1:nrow(c),labels=tissue_names)
# }
# 
# sub_pm <- pm[sample(nrow(pm), 300), ]
# apm <- aperm(sub_pm)
# cors <- cor(apm)
# pheatmap(cors)
