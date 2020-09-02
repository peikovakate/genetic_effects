suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
source("scripts/utils.R")

parser <- OptionParser()
parser <- add_option(parser, c("-v", "--variants"), type = "character", help="file with variants")
parser <- add_option(parser, c('-s', "--sumstat_file"), 
                     type = "character", help = "eqtl catalogue sumstat file")
parser <- add_option(parser, c('-p', '--pattern'),
                     default = "_microarray.nominal.sorted.txt.gz",
                     type="character", 
                     help = "sumstat file pattern ending similar for all tissues")
parser <- add_option(parser, c('-o', '--output_dir'),
                     type="character", 
                     help = "output dir")
args = parse_args(parser)

sumstat <- args$sumstat_file
pattern <- args$pattern
variants_file <- args$variants
output <- args$output_dir

# sumstat <- "../data/sumstat_2/CEDAR.B-cell_CD19_microarray.nominal.sorted.txt.gz"
# pattern <- args$pattern
# variants_file <- "../data/credible_sets/credible_set_variants.tsv"

print(paste("Sumstats file", sumstat))
print(paste("Sumstats file pattern", pattern))
print(paste("Read file with variants", variants_file))

variants = read_tsv(variants_file)
head(variants)

cell_type = str_split(basename(sumstat), pattern, simplify = T)[1]

region <-  GRanges(
  seqnames = variants$chr,
  ranges = IRanges(start = variants$pos,
                   end = variants$pos)
)

print("Start tabix scan")
start_time <- Sys.time()

regions <- scanTabixDataFrame(sumstat,
  region,
  col_names = eqtl_colnames,
  col_types = cols(alt = "c", ref = "c")
)

print(Sys.time() - start_time)
print("Join and postporcess regions")
start_time <- Sys.time()

variant_regions = bind_rows(regions) %>% distinct(.keep_all = T)

extracted_sumstat = 
  inner_join(variant_regions, 
             variants[c("variant_id", "phenotype_id")],
             by=c("variant"="variant_id", "molecular_trait_id"="phenotype_id")) %>% 
  distinct(.keep_all = T) %>% 
  mutate(cell_type = cell_type)

print(Sys.time() - start_time)
write_tsv(extracted_sumstat, file.path(output, paste0(cell_type, ".tsv")))
print("Done")

