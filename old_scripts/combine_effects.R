library(tidyverse)

sumstat_across_tissues_file = "../data/rnaseq/variants_across_tissues_rnaseq.tsv"
output_file = "../data/rnaseq/rnaseq_effects.tsv"

init_tbl <- read_tsv(sumstat_across_tissues_file)

tbl <- init_tbl %>%
  select(variant,
         molecular_trait_id,
         cc_id,
         pvalue,
         beta,
         se,
         cell_type) %>% 
  mutate(eqtl_id = paste(variant, molecular_trait_id, sep="."))

tbls <- tbl %>%
  group_by(cell_type) %>%
  group_split()

suffixs <- sapply(tbls, function(x) {
  cell_type_name = x[1,]$cell_type
  sprintf(".%s", cell_type_name)
})

effects <- tbls[[1]]
for (i in 2:length(tbls)) {
  subtbl = tbls[[i]][c("eqtl_id","pvalue", "beta", "se")]
  if (i %% 2 == 0) {
    effects <-
      inner_join(
        effects,
        subtbl,
        by = "eqtl_id",
        suffix = suffixs[c(i - 1, i)]
      )
  } else{
    effects <-
      inner_join(effects, subtbl, by = "eqtl_id")
  }
  print(tbls[[i]][1, ]$cell_type)
}

if (length(tbls) %% 2 == 1) {
  l = length(tbls)
  b_name = paste0("beta", suffixs[[l]])
  se_name = paste0("se", suffixs[[l]])
  pvalue_name = paste0("pvalue", suffixs[[l]])
  effects <-
    dplyr::rename(effects,
                  !!b_name := beta,
                  !!se_name := se,
                  !!pvalue_name := pvalue)
}
nrow(effects)
effects = mutate(effects, eqtl_id = paste(variant, molecular_trait_id, sep="."))

write_tsv(effects, output_file)
