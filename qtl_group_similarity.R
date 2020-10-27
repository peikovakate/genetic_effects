library(tidyverse)
library(mashr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(plotly)
source("utils2.R")

colours = read_tsv("../data2/gtex/qtl_groups_colours.tsv")
metadata = read_tsv("../data2/gtex/qtl_groups_data.tsv")
leads = read_tsv("../data2/gtex/qtl_groups_leads.tsv")
metadata = left_join(metadata, colours)

plot_coords = function(coords){
  coords = tibble(x = coords$points[,1], y=coords$points[,2], qtl_group=rownames(coords$points))
  coords = left_join(coords, metadata[c("qtl_group", "colour", "group")])
  coords
  plt = ggplot2::ggplot(coords, aes(x, y, colour=group, grp=qtl_group)) +
    ggplot2::geom_point() +
    scale_colour_manual(values=colours$colour)  +
    # geom_label_repel(aes(label = qtlGroup), size=3) +
    theme_classic()
  fig <- ggplotly(plt, tooltip = c("qtl_group"))
  return(fig)
}


effects_file <- "../data2/microarr/subset_lead_effects_na.tsv"
effects <- readr::read_tsv(effects_file)
nrow(effects)

eqtls = effects_to_matricies(effects, replace_na_with="zero")
cor_method = "spearman"
cols.cor <- cor(eqtls$beta, method = cor_method)
pheatmap::pheatmap(cols.cor, cutree_rows = 5, colorRampPalette(brewer.pal(n = 7, name ="OrRd"))(100),
         fontsize=12,border_color=NA,
         filename=sprintf("figures/microarr/%s_12k_zero_pheatmap.png", cor_method), width=12, height=10)

coords = MASS::isoMDS(1-cols.cor)
plot_coords(coords)

########## sharing ##########

load("../data2/microarr/mash_13k_na_to_mean_mash_ed_sharing.R")

pheatmap(sharing, cutree_rows = 5, colorRampPalette(brewer.pal(n = 7, name ="OrRd"))(100),
         fontsize=12,border_color=NA,
         filename="figures/microarr/mash_13k_mean_ed_heatmap.png", width=12, height=10)

coords = MASS::isoMDS(1-sharing)
plot_coords(coords)



