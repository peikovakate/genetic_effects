library(tidyverse)
library(mashr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
source("utils2.R")


effects_file <- "../data2/gtex/lead_effects_na.tsv"

effects <- read_tsv(effects_file)
nrow(effects)

eqtls = effects_to_matricies(effects, replace_na_with="zero")
cor_method = "spearman"
cols.cor <- cor(eqtls$beta, method = cor_method)
pheatmap(cols.cor, cutree_rows = 7, colorRampPalette(brewer.pal(n = 7, name ="OrRd"))(100),
         fontsize=7,border_color=NA,
# )
         filename=sprintf("figures/gtex/%s_55k_zero_pheatmap.png", cor_method), width=12, height=10)

coords = MASS::isoMDS(1-cols.cor)

plot_coords = function(coords){
  coords = tibble(x = coords$points[,1], y=coords$points[,2], qtlGroup=rownames(coords$points))
  plt = ggplot2::ggplot(coords, aes(x, y)) +
    ggplot2::geom_point() +
    geom_label_repel(aes(label = qtlGroup), size=3) +
    theme_classic()
  return(plt)
}

plot_coords(coords)



pheatmap(sharing, cutree_rows = 7, colorRampPalette(brewer.pal(n = 7, name ="OrRd"))(100),
         fontsize=7,border_color=NA,
         filename=sprintf("figures/gtex/mash_ed_only_heatmap.png", cor_method), width=12, height=10)

coords = MASS::isoMDS(1-sharing)
plot_coords(coords)



