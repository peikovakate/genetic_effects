library(pheatmap)
library(tidyverse)
source("../thesis_experiments/pipeline/visualisation_utils.R")

colors = read_tsv("../data/mfactorization/unique_genes/tissues_colors.txt")

k = 40
a = 800
l = 600
r = 17

matrix_file = "../data/mfactorization/output/sn_spMF_K15_a1800_l1300/sn_spMF_K15_a1800_l1300_Run28.RData"
matrix_file = sprintf("../data/mfactorization/output_rnaseq/sn_spMF_K%i_a1%i_l1%i/sn_spMF_K%i_a1%i_l1%i_Run%i.RData", k,a,l,k,a,l,r)
matrix_file = "../data/mfactorization/output_rnaseq_with_na/sn_spMF_K39_a1700_l1700/sn_spMF_K39_a1700_l1700_Run4.RData"
# output_rnaseq_with_na/sn_spMF_K39_a1700_l1700/sn_spMF_FactorMatrix_K39_a1700_l1700_Run4.txt

load(matrix_file)
p = pheatmap(FactorM, 
         clustering_distance_rows = "correlation",
         file="../thesis_experiments/final_figures/rnaseq/factor_heatmap.png")

  my_gtable = p$gtable

my_gtable$grobs[[3]]$gp=gpar(col="#ffffff", fontsize=20)

hmap(
  as.data.frame(FactorM),
  fontSize = 6.5,
  tileColor="white",
  # colorPalette = 7,
  dendroLineSize = 0.2,
  scaleName="Loading",
  revColors = T
  # main="MashR pairwise tissue sharing"
)

# ggsave("../thesis_experiments/final_figures/rnaseq/factors_with_na.png", scale=1, width = 6, height = 5)

