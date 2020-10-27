library(tidyverse)
library(data.table)
source("utils2.R")
# ccs = read_tsv("../data2/gtex/cc.tsv")
# ccs %>% select(phenotype_id, variant_id, qtl_group, cs_index) %>% write_tsv("../data2/gtex/credible_sets_across_datasets.tsv")

ontology_map = readr::read_tsv("../eQTL-Catalogue-resources/ontology_mappings/tissue_ontology_mapping.tsv")
friendly_names = readr::read_tsv("../eQTL-Catalogue-resources/ontology_mappings/friendly_names.tsv") %>%
  dplyr::select(ontology_term, ontology_tissue)
# rename blueprint dataset for consistency with effect dataset names
ontology_map$study[c(1,2,3)] = c("BLUEPRINT_SE", "BLUEPRINT_SE", "BLUEPRINT_PE")
ontology_map <- ontology_map %>% dplyr::mutate(study_qtlgroup = paste0(study, ".", qtl_group)) %>%
  dplyr::left_join(friendly_names) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "brain", "brain", "other")) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "LCL", "LCL", sample_class)) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "monocyte", "monocyte", sample_class)) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "macrophage", "macrophage", sample_class)) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "blood", "blood", sample_class)) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "neutrophil", "neutrophil", sample_class)) %>%
  dplyr::mutate(sample_class = ifelse(ontology_term %in% c("CL_0000236","CL_0002677","CL_0002678","CL_0000624","CL_0000625","CL_0000623","CL_0000899","CL_0000546","CL_0000545","CL_0000899","CL_0002038","CL_0000084"), "lymphocyte", sample_class)) %>%
  dplyr::mutate(sample_class = ifelse(ontology_tissue %like% "iPSC", "iPSC", sample_class))

fct_levels = c("blood","lymphocyte","LCL","neutrophil","monocyte","macrophage","brain","iPSC","other")
ontology_map = dplyr::mutate(ontology_map, tissue_fct = factor(sample_class, levels = fct_levels))
# place immune t-cells together
ontology_map = ontology_map %>% dplyr::mutate(dummy = ifelse(qtl_group %like% "anti", "ANTI", ontology_tissue), 
                                              ontology_tissue = ifelse(qtl_group %like% "anti", paste(ontology_tissue, "(anti-CD3-CD28)"), ontology_tissue))
# place quach monocyte and nedelec macrophage together because they load on one factor
ontology_map = ontology_map %>% dplyr::mutate(dummy = ifelse(study_qtlgroup == "Quach_2016.monocyte_naive", "z-monocyte", dummy))
ontology_map = ontology_map %>% dplyr::mutate(dummy = ifelse(study_qtlgroup == "Nedelec_2016.macrophage_naive", "a-macrophage", dummy))
# sort studies
ontology_map = ontology_map %>% dplyr::arrange(tissue_fct, dummy, study)
ontology_map = ontology_map %>% dplyr::mutate(heatmap_label = paste(study, ontology_tissue, sep=" "))

credible_sets = read_tsv("../data2/gtex/credible_sets_across_datasets.tsv")
loadings_file = "../data2/gtex/mfactorization/mapping/sn_spMF_K30_a1900_l11100_Loadings_beta_alpha0.05_corrected.txt"
loadings = read.delim(loadings_file)
loadings = loadings_to_tibble(loadings)

factor_names = tibble(factor = paste0("Factor", c(1:16)),
                      name = c("Universal", "iPSC", "Skin", "Blood",
                               "LCL", "T-cell", "Monocyte & Macrophage", "Fat",
                               "T-cell immune response", "Brain", "Macrophage", "Monocyte naive",
                               "BLUEPRINT T-cell", "ROSMAP Brain", "Muscle", "Neutrophil"))

# loadings = loadings %>% select(starts_with("Factor"))
# colnames(loadings) = factor_names$name
loadings = loadings %>% rename_at(vars(starts_with("Factor")), ~ factor_names$name)

tissue_eqtl_counts = loadings %>% left_join(credible_sets, by=c("variant"="variant_id", "gene"="phenotype_id")) %>% 
  group_by(qtl_group) %>% summarise(total_tissue_eqtls = n()) %>% ungroup()

factors = tidyr::pivot_longer(loadings, cols = factor_names$name, names_to = "Factors", values_to = "Loadings")
factors = factors %>% filter(Loadings != 0)

loadings %>% filter(Universal != 0) %>% nrow()

factors = factors %>% left_join(credible_sets, by=c("variant"="variant_id", "gene"="phenotype_id"))
# factors = factors %>% 
factor_counts = factors %>% group_by(qtl_group, Factors) %>% summarise(count = n()) %>% ungroup()
factor_counts = factor_counts %>% left_join(tissue_eqtl_counts)

datasets = group_split(factor_counts, qtl_group)

lapply(datasets, function(data){
  data = data %>% mutate(fraction = count/total_tissue_eqtls)
  plt = ggplot(data, aes(Factors, fraction)) +
    geom_col() + 
    ggtitle(data$qtl_group[1]) +
    ylab("eQTLs assigned %") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), 
          text = element_text(size=16), panel.grid.major = element_blank())
  ggsave(file.path("figures/dataset_factor/", paste0(data$qtl_group[1], ".png")), plt, width=8, height = 6)
})

universal_factor_counts = factor_counts %>% filter(Factors=="Universal")
universal_factor_counts = universal_factor_counts %>% left_join(ontology_map[c("study_qtlgroup", "tissue_fct")], by=c("qtl_group"="study_qtlgroup"))

# total_tissue_counts = factors %>% group_by(qtl_group) %>% summarise(total_count = n()) %>% ungroup()

# universal_factor_counts = universal_factor %>% group_by(tissue_fct, qtl_group) %>% summarise(count=n()) %>% ungroup()
# universal_factor_counts = universal_factor_counts %>% left_join(total_tissue_counts)
universal_factor_counts = universal_factor_counts %>% mutate(proportion = count / total_tissue_eqtls)
plt = ggplot(universal_factor_counts, aes(tissue_fct, proportion, qtl_group=qtl_group, count=count)) +
  geom_jitter(width = 0.2) +
  theme_light() 

plotly::ggplotly(plt)

# factors %>% group_by(qtl_group, Factors) %>% summarise(count = n())
# universal_variants = factors %>% filter(Factors == "Factor1", Loadings != 0)
# universal_variants = universal_variants %>% 
#   left_join(credible_sets, by=c("variant"="variant_id", "gene"="phenotype_id"))
# universal_variants = universal_variants %>% group_by(qtl_group) %>% summarise(count=n())
# universal_variants = inner_join(ontology_map, universal_variants, by=c("study_qtlgroup"="qtl_group"))
# 
# specific_variants = factors %>% filter(Factors!="Factor1", Loadings != 0)
# specific_variants = specific_variants %>% left_join(credible_sets, by=c("variant"="variant_id", "gene"="phenotype_id"))
# specific_variants = specific_variants %>% group_by(qtl_group) %>% summarise(count=n())
# specific_variants = inner_join(ontology_map, specific_variants, by=c("study_qtlgroup"="qtl_group"))
# 
# ggplot(universal_variants, aes(qtl_group, count)) +
#   geom_col() + 
#   theme_light() +
#   xlab("Dataset") +
#   ylab("eQTLs assigned") +
#   theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), 
#         text = element_text(size=16), panel.grid.major = element_blank())
# 
# ggplot(specific_variants, aes(qtl_group, count)) +
#   geom_col() + 
#   xlab("Dataset") +
#   ylab("eQTLs assigned") +
#   theme_light() +
#   theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), 
#         text = element_text(size=16), panel.grid.major = element_blank())
