library(tidyverse)

variants_folder = "../data/sumstat_cs_variants/"

files = list.files(variants_folder)
tbls = lapply(file.path(variants_folder, files), read_tsv)
tbls = bind_rows(tbls)
head(tbls)

write_tsv(tbls, "../data/sumstat_cs_variants/tissues_combined.tsv")
