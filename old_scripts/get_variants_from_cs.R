library(tidyverse)

susie_path = "../data/susie/"
suffix = ".HumanHT12V4.txt.gz"
output_file = "../data/credible_sets/credible_set_variants.tsv"

files <- list.files(susie_path)

# load susie files
tbls <- lapply(files, function(x){
  file = paste0(susie_path, x)
  print(file)
  cell_type_name = strsplit(x, split=suffix)[[1]]
  print(cell_type_name)
  tbl <- read.table(file, header = T, sep='\t', stringsAsFactors = F) %>% 
    mutate(cell_type = cell_type_name)
  return(tbl)
})

tbls <- bind_rows(tbls)

# assign to each pair a cc_id, 
# then it's easy to use the same code for querying eQTL Catalogue sumstat, 
# that was written for connected components
tbls = tbls %>% 
  distinct(variant_id, phenotype_id, .keep_all = T) %>% 
  mutate(cc_id = paste(variant_id, phenotype_id, sep="."))

tbls %>% nrow()

write_tsv(tbls, output_file)
