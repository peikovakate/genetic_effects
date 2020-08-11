library(igraph)
library(tidyverse)
library(ggplot2)
library(gridGraphics)
library(cowplot)

# susie_path = "../data/susie/"
susie_path = "../data/susie_rnaseq/"
# suffix = ".HumanHT12V4.txt.gz"
suffix = ".ge.txt.gz"
# output_file = "../data/connected_components/microarr_cc.tsv"
output_file = "../data/connected_components/rnaseq_cc.tsv"

files <- list.files(susie_path, pattern=suffix)

# load susie files
tbls <- lapply(files, function(x){
  file = file.path(susie_path, x)
  print(file)
  cell_type_name = strsplit(x, split=suffix)[[1]]
  print(cell_type_name)
  tbl <- read.table(file, header = T, sep='\t', stringsAsFactors = F) %>% 
    mutate(cell_type = cell_type_name)
  return(tbl)
})

# unique_phenotypes <- lapply(tbls, function(x){unique(x$phenotype_id)})

tbls <- bind_rows(tbls)
# make unique id for credible sets
tbls <- mutate(tbls, cs_uid = paste(phenotype_id, cs_id, cell_type, sep="_"))
nrow(tbls)

tissue_names = unique(tbls$cell_type)


rnaseq_tbls = tbls
microarr_tbls = tbls

# cs summary per tissue
tissue_counts_microarr = tbls %>% 
  group_by(cell_type) %>% 
  summarise(genes = n_distinct(phenotype_id), credible_sets = n_distinct(cs_uid)) %>% 
  pivot_longer(cols=c(genes, credible_sets), values_to="Count", names_to="Term")

mp = ggplot(tissue_counts_microarr, aes(x=cell_type, y=Count, fill=Term)) +
  geom_col(position = position_dodge()) + 
  scale_fill_brewer(palette = "Set2") +
  xlab("Conditions") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14))

rp = ggplot(tissue_counts_rnaseq, aes(x=cell_type, y=Count, fill=Term)) +
  geom_col(position = position_dodge()) + 
  # scale_fill_manual(values = c("darkgreen", "darkblue"))+
  scale_fill_brewer(palette = "Set2") +
  xlab("Conditions") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14))

plt = plot_grid(mp+theme(legend.position = "none"), 
                rp+theme(legend.position = "none"), 
                ncol=2, align = "h", rel_widths=c(2,1), labels="AUTO")

legend <- get_legend(mp)
plot_grid(plt, legend, rel_widths = c(6, 1))


plot_grid(mp, rp, ncol=2, align = "h", rel_widths = c(2,1), labels="AUTO")
plt

counts <- tbls %>% 
  group_by(phenotype_id) %>% 
  dplyr::distinct(phenotype_id, cell_type, cs_id) %>% 
  summarise(count = n())

# ggplot(counts, aes(x=count)) +
#   ggplot2::geom_histogram() +
#   geom_density(alpha=0.6)

z_threshold <- 3
cs_size_threshold <- 50

# group into unique credible sets
# filter by max z score in the credible sets and credible set size
# for each phenotype count distinct credible sets in all tissues together
credible_sets <- tbls %>% 
  group_by(cs_uid) %>% 
  mutate(max_abs_z = max(abs(z))) %>% 
  filter(max_abs_z > z_threshold, cs_size < cs_size_threshold) %>% 
  ungroup()

phenotypes <- credible_sets %>% 
  group_by(phenotype_id) %>% 
  distinct(cs_id, cell_type) %>% 
  summarise(count = n()) 

ggplot(phenotypes, aes(count)) +
  geom_histogram(binwidth = 1) +
  labs(title = "N credible sets per phenotype in all tissues and cell types together") +
  stat_bin(binwidth=1, geom="text", colour="black", size=3.5, aes(label=..count..), vjust=-0.5)

nrow(credible_sets)
max(credible_sets$cs_size)



################################
### Find connected components
################################

get_cc <- function(phenotype_group){
  # phenotype_group <-  ungroup(phenotype_group)
  variants <- csDfToList(phenotype_group)
  if(length(variants) < 2){
    combinations <- matrix(c(1,1), nrow=2)
  }else{
    combinations <- combn(length(variants), 2, simplify = T)
  }
  
  # combinations[1,]
  intersections <- map2(combinations[1,], combinations[2,], function(x, y){
    intr <- intersect(variants[[x]], variants[[y]])
    return(length(intr)>0)
  })
  
  n_verts <-  length(variants)
  self_edges <- matrix(c(1:n_verts, 1:n_verts), nrow=n_verts)
  combs <- t(combinations[, unlist(intersections)])
  edge_list <- rbind(self_edges, combs)
  
  g <- graph_from_edgelist(edge_list)
  # g <- as.undirected(g)
  
  components_sets <- 1:components(g)$no %>% 
    lapply(function(x){return(names(variants)[components(g)$membership == x])})
  
  res <- lapply(components_sets, function(x){
    phenotype_group %>% 
      dplyr::filter(cs_uid %in% x) %>% 
      dplyr::mutate(cc_id = paste(x, collapse=""))
  }) %>% bind_rows()
  
  return(res)
  
}

csDfToList <- function(cs_df){
  grouped_df = dplyr::group_by(cs_df, cs_uid)
  cs_list = setNames(dplyr::group_split(grouped_df), dplyr::group_keys(grouped_df)$cs_uid) %>%
    purrr::map(~.$variant_id)
  return(cs_list)
}

groups <- credible_sets %>% 
  group_by(phenotype_id) %>% 
  group_split()

connected_components <- lapply(groups, get_cc) %>% 
  bind_rows()

# saving connected components to file
connected_components %>% 
  write_tsv(output_file, quote_escape=F)

# check that function for cc search works
h <- credible_sets %>% filter(phenotype_id == "ILMN_1655990")
get_cc(h)

# this threshold is not applied to the final file
# used only for plotting
pip_threshold <- 0.1

# just plotting
genes <- connected_components %>% 
  group_by(phenotype_id)%>% 
  summarise(distinct = length(unique(cc_id)))

ggplot(genes, aes(distinct)) +
  geom_histogram(binwidth = 1) + 
  labs(title = "Count connected components per phenotype ") + 
  xlab("N of conncted components in phenotype") +
  ylab("Phenotype counts") +
  stat_bin(binwidth=1, geom="text", colour="black", size=3.5, aes(label=..count..), vjust=-0.5)

genes <- connected_components %>% 
  filter(pip < pip_threshold) %>% 
  group_by(phenotype_id)%>% 
  summarise(distinct = length(unique(cc_id)))

ggplot(genes, aes(distinct)) +
  geom_histogram(binwidth = 1) + 
  labs(title = sprintf("Count connected components (pip < %.2f) per phenotype ", pip_threshold)) + 
  xlab("N of conncted components in phenotype") +
  ylab("Phenotype counts") +
  stat_bin(binwidth=1, geom="text", colour="black", size=3.5, aes(label=..count..), vjust=-0.5)
