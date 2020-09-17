library(tidyverse)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(plotly)

colours = read_tsv("../data2/gtex/qtl_groups_colours.tsv")
metadata = read_tsv("../data2/gtex/qtl_groups_data.tsv")
leads = read_tsv("../data2/gtex/qtl_groups_leads.tsv")

load("../data2/gtex/mash_55k_mash_ed_sharing.R")

metadata = left_join(metadata, colours)

plot_distribution <- function(similarity_matrix, y_axis="Mash sharing", measure="sharing"){
  df_sharing = tibble(reshape2::melt(similarity_matrix))
  
  df_sharing = df_sharing %>% filter(value < 1)
  df_sharing = inner_join(df_sharing, leads, by=c("Var1"="lead"))
  df_sharing = left_join(df_sharing, metadata[c("qtl_group", "colour", "group")], by=c("Var2"="qtl_group"))
  df_sharing = dplyr::rename(df_sharing, lead=Var1, qtl_group=Var2, !!measure:=value)
  
  plt <- ggplot(df_sharing, aes(lead, get(measure), colour=group, grp=qtl_group)) +
    # geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.008)+
    geom_jitter(width = 0.2) +
    xlab("QTL group") +
    ylab(y_axis) +
    scale_colour_manual(values=colours$colour)  +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  fig <- ggplotly(plt, tooltip = c("qtl_group", measure))
  return(fig)
}


plot_distribution(sharing)


