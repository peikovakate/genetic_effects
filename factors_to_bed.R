library(rtracklayer)
library(tidyverse)
library(GenomicRanges)
library(dplyr)

chain <- import.chain("../data/hg38ToHg19.over.chain")
eqtls_one_factor_per_line = read_tsv("../data/found_variants_in_sumstat/found_variants_factors.tsv")

regions <-  GRanges(
  seqnames = paste0("chr", eqtls_one_factor_per_line$chromosome),
  ranges = IRanges(start = eqtls_one_factor_per_line$position,
                   end = eqtls_one_factor_per_line$position),
  mcols = eqtls_one_factor_per_line$factor
)
regions_hg19 <- liftOver(regions, chain)
regions_hg19_df = as_tibble(regions_hg19)

factors <- regions_hg19_df %>% 
  dplyr::rename(factor = mcols) %>% 
  select(seqnames, start, end, factor) %>% 
  group_split(factor)

lapply(factors, function(df){
  factor = df[1,]$factor
  df = select(df, seqnames, start, end)
  write_tsv(df, sprintf("../data/found_variants_in_sumstat/factor_annotations/factor%i.bed", factor), col_names = F)
})

factors


