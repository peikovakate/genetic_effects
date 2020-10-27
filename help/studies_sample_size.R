library(tidyverse)
library(ggplot2)

dir_path = "../SampleArcheology/studies/cleaned"

files = list.files(dir_path, full.names = F)
files = files[!files %in% c("GEUVADIS_EUR.tsv", "BLUEPRINT_PE.tsv", "BLUEPRINT_SE.tsv")]

tbls = lapply(files, function(file){
  tbl = read_tsv(file.path(dir_path, file), col_types = cols_only(rna_qc_passed="l", 
                                                                  genotype_qc_passed="l", 
                                                                  sample_id="c", 
                                                                  qtl_group="c", 
                                                                  cell_type="c", 
                                                                  type="c", 
                                                                  study="c"))
  return(tbl)
})

studies = bind_rows(tbls)
studies = studies %>% filter(rna_qc_passed == TRUE, genotype_qc_passed == TRUE)
studies[is.na(studies$type), "type"] = "RNA-seq"

condition = read_tsv("../eQTL-Catalogue-resources/tabix/tabix_ftp_paths.tsv")
condition = select(condition, study, qtl_group, condition_label) %>% distinct()
studies = left_join(studies, condition)

naive_studies = studies %>% mutate(condition_label = if_else(study=="GTEx", "naive", condition_label))
naive_studies = naive_studies %>% filter(condition_label == "naive")

sample_size = naive_studies %>% group_by(cell_type, qtl_group, study, type) %>%  summarise(dataset_sample_size=n())
sample_size %>% write_tsv("dataset_sample_sizes.tsv")

# plotting
sample_size = readr::read_tsv("dataset_sample_sizes.tsv")

ontology =  readr::read_tsv("../eQTL-Catalogue-resources/ontology_mappings/tissue_ontology_mapping.tsv")
friendly_names =  readr::read_tsv("../eQTL-Catalogue-resources/ontology_mappings/friendly_names.tsv")

ontology = ontology %>% dplyr::left_join(friendly_names[c("ontology_term", "ontology_tissue")], by="ontology_term")
ontology = ontology %>% dplyr::arrange(ontology_tissue, study)

studies = dplyr::left_join(sample_size, ontology)

rnaseq_studies = dplyr::filter(studies, type == "RNA-seq")
microarray_studies = dplyr::filter(studies, type == "microarray", study != "Kolberg_2020")

draw_plot = function(studies, legend_rows = 4, legen_y_position=0.8){
  sample_sizes = studies %>% 
    dplyr::group_by(ontology_tissue) %>% 
    dplyr::mutate(total_sample_size=sum(dataset_sample_size))
  
  ggplot(sample_sizes, aes(x = reorder(ontology_tissue, -total_sample_size), dataset_sample_size, fill=study)) +
    geom_col() + 
    guides(fill=guide_legend(nrow=legend_rows,byrow=TRUE))+
    xlab("Tissue (naive)") + 
    ylab("Sample size") +
    theme_light() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 10),
          panel.grid = element_blank(),
          legend.position=c(0.5, legen_y_position),
          legend.background = element_rect(colour="lightgrey", 
                                           size=0.5, linetype="solid"))
}
#6x12inch
rnaseq_plt = draw_plot(rnaseq_studies)
ggsave("rnaseq_sample_size.pdf", rnaseq_plt, width = 12, height = 6)
#6x9
microarr_plt = draw_plot(microarray_studies, legend_rows = 1, legen_y_position=0.9)
ggsave("microarr_sample_size.pdf", microarr_plt, width = 9, height = 6)



