library(tidyverse)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

colours = read_tsv("../data2/gtex/qtl_groups_colours.tsv")
metadata = read_tsv("../data2/gtex/qtl_groups_data.tsv")
leads = read_tsv("../data2/gtex/qtl_groups_leads.tsv")


metadata = left_join(metadata, colours)

sharing
df_sharing = tibble(reshape2::melt(sharing))
df_sharing = df_sharing %>% filter(value < 1)

df_sharing = inner_join(df_sharing, leads, by=c("Var1"="lead"))
df_sharing = left_join(df_sharing, metadata[c("qtl_group", "colour", "group")], by=c("Var2"="qtl_group"))

ggplot(df_sharing, aes(Var1, value, colour=group)) +
# geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.008)+
  geom_jitter(width = 0.2) +
  xlab("QTL group") +
  ylab("Sharing") +
  scale_colour_manual(values=colours$colour)  +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
