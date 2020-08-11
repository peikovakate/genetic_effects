library(forestplot)
library(rmeta)
genes <- read_tsv("additional_data/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv")
cluster_file <- read_tsv("additional_data/Figure_2a_response_eQTLs.txt")

cluster_genes = filter(cluster_file, new_cluster_id == 3)$gene_name

eqtls <- lapply(cluster_genes, function(gene){
  genes_info <- genes %>%
    filter(gene_name == gene) 
  if(nrow(genes_info) > 0){
    eqtls <- effects %>%
      filter(molecular_trait_id %in% genes_info$phenotype_id) %>% 
      mutate(gene_name = gene)
    return(eqtls)
  }
})

eqtls <- bind_rows(eqtls)


ind = 11
# ind = 53

eqtl = eqtls[ind, ]
# eqtl = effects[ind,]

betas <-
  select(eqtl, starts_with('beta.')) %>% rename_all(function(x) {
    sub("beta.", "", x)
  })

ses <-
  select(eqtl, starts_with('se.')) %>% rename_all(function(x) {
    sub("se.", "", x)
  })

zs <- 
  select(eqtl, starts_with('z.')) %>% rename_all(function(x) {
    sub("z.", "", x)
  })

pair_ind = which(
  effects$variant == eqtl$variant &
    effects$molecular_trait_id == eqtl$molecular_trait_id
)

# xx <-
#   tibble(mean = as.numeric(zs)) %>% mutate(lower = mean, upper = mean)
# 
# forestplot::forestplot(names(r),
#                        xx$mean,
#                        xx$lower - 1,
#                        xx$upper + 1,
#                        title = paste("prior", gene, eqtl$variant))
par(mar=c(5,4,4,2)+0.1)
metaplot(
  as.numeric(betas),
  se =  m$result$PosteriorSD[pair_ind,],
  labels = names(betas),
  ylab = "Cell type",
  xlab = paste("Prior Beta:", eqtl$gene_name, eqtl$variant)
)

# mash_plot_meta(m, pair_ind, labels = colnames(m$result$PosteriorMean))
metaplot(
  m$result$PosteriorMean[pair_ind,],
  se = m$result$PosteriorSD[pair_ind,],
  labels = colnames(m$result$PosteriorMean),
  ylab = "Cell type",
  xlab = paste("Posterior Beta:", eqtl$gene_name, eqtl$variant)
)
# MGAT3 - B-cell specific
# TRAF1 - LPS2 hours (rs10985070 )
# SLFN5 - stimulated monocytes
# ARHGEF3 - platelet specific