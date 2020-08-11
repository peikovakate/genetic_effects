suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRanges))

chain_file = "../data/hg38ToHg19.over.chain"
# factors_tbl_input =  "../data/found_variants_in_sumstat/eqtl_factors_for_cont_annot.tsv"

parser <- OptionParser()
parser <- add_option(parser, c("-b", "--bim"), type="character")
parser <- add_option(parser, c('-o', '--out'), type="character")
parser <- add_option(parser, c('-f', '--factors'), type="character")
args = parse_args(parser)

bim_file_path = args$bim
output_dir = args$out
factors_tbl_input = args$factors

# bim_file_path = "../data/bim_files/eur_chr10_biall_maf_cm.bim"
# output_dir = "../data/found_variants_in_sumstat/factor_annotations/10/"
# factors_tbl_input = "../data/found_variants_in_sumstat/eqtl_factors_for_cont_annot.tsv"

dir.create(output_dir, recursive = T)

eqtl_factors = read_tsv(factors_tbl_input)
chain <- import.chain(chain_file)

bim_file = read_tsv(bim_file_path, 
                    col_names = c("chr", "variant", "CM", "pos", "ref", "alt"), 
                    col_types = cols(chr = col_character()))

bim_file = bim_file %>% mutate(len_ref = sapply(ref, nchar), len_alt = sapply(alt, nchar))
bim_file = bim_file %>% mutate(end_pos = pos + len_ref)

bim_ranges = GRanges(
  seqnames = paste0("chr", bim_file$chr),
  ranges = IRanges(start = bim_file$pos, end = bim_file$end_pos)
)

annotate_factor <- function(one_factor_eqtls){
  f = as.numeric(one_factor_eqtls[1, "factor"])
  one_factor_eqtls = eqtl_factors %>% filter(factor == f)
  one_factor_eqtls %>% distinct(eqtl_id) %>% nrow()
  print(paste("factor", f, "has", nrow(one_factor_eqtls), "eqtls"))
  
  regions <-  GRanges(
    seqnames = paste0("chr", one_factor_eqtls$chromosome),
    ranges = IRanges(start = one_factor_eqtls$position,
                     end = one_factor_eqtls$position),
    mcols = select(one_factor_eqtls, factor, stand_beta, eqtl_id)
  )
  
  regions_hg19 <- liftOver(regions, chain)
  regions_hg19_df = as.data.frame(regions_hg19)
  regions_hg19 <- GRanges(regions_hg19_df)
  
  overlaps = GenomicRanges::findOverlaps(bim_ranges, regions_hg19)
  overlaps = as.data.frame(overlaps)
  
  bim_file = mutate(bim_file, stand_beta = 0)
  bim_file$stand_beta[overlaps$queryHits] <- regions_hg19[overlaps$subjectHits, ]$mcols.stand_beta
  
  category_name = paste0("factor", f)
  
  annot = bim_file %>% 
    dplyr::rename(CHR = chr, SNP = variant, BP = pos, !!(category_name) := stand_beta) %>% 
    dplyr::select(CHR, SNP, CM, BP, !!category_name)
  
  write_tsv(annot, file.path(output_dir, paste0(category_name, ".bed.annot")))
  return()
}

factors_dfs <- dplyr::group_split(eqtl_factors, factor)
lapply(factors_dfs, annotate_factor)
