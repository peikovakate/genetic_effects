`%>%` <- magrittr::`%>%`
source("utils.R")

parser <- optparse::OptionParser()
parser <- optparse::add_option(parser, c("-v", "--variants"), 
                               type = "character", 
                               help="file with variants, should have phenotype_id, variant_id, cc_id columns")
parser <- optparse::add_option(parser, c('-s', "--sumstat_file"), 
                               type = "character", 
                               help = "eqtl catalogue sumstat file, with molecular_id column")
parser <- optparse::add_option(parser, c('-p', '--pattern'),
                               default = "_ge.nominal.sorted.tsv.gz",
                               type="character", 
                               help = "sumstat file pattern ending, similar for all tissues")
parser <- optparse::add_option(parser, c('-o', '--output_dir'),
                               type="character", 
                               help = "output dir")
args = optparse::parse_args(parser)
# 
sumstat <- args$sumstat_file
pattern <- args$pattern
variants_file <- args$variants
output_dir <- args$output_dir

# sumstat <- "../data/sumstat_2/TwinsUK.skin_ge.nominal.sorted.tsv.gz"
# pattern <- args$pattern
# variants_file <- "../data2/connected_components/rnaseq_cc.tsv"
# output_dir <- "../data2/sumstat/"

print(paste("Sumstats file", sumstat))
print(paste("Sumstats file pattern", pattern))
print(paste("Read file with variants", variants_file))

variants = readr::read_tsv(variants_file, col_types = readr::cols())
# some of the variants can be duplicated, we don't want to query them multiple times
unique_variants = dplyr::distinct(variants, variant_id, pos, chr)
# unique_variants = unique_variants[1:5000, ]
# head(variants)

# extract the name of QTL group
qtlGroup = stringr::str_split(basename(sumstat), pattern, simplify = T)[1]

region <- GenomicRanges::GRanges(
  seqnames = unique_variants$chr,
  ranges = IRanges::IRanges(start = unique_variants$pos,
                   end = unique_variants$pos)
)

print("Start tabix scan")
start_time <- Sys.time()

regions <- scanTabixDataFrame(sumstat,
  region,
  col_names = eqtl_colnames,
  col_types = readr::cols(alt = "c", ref = "c")
)

print(Sys.time() - start_time)

print("Join and postporcess regions")
start_time <- Sys.time()

variant_regions = dplyr::bind_rows(regions) %>% 
  dplyr::distinct(.keep_all = T)

extracted_sumstat = 
  dplyr::inner_join(variant_regions, 
             variants[c("variant_id", "phenotype_id", "cc_id")],
             by=c("variant"="variant_id", "molecular_trait_id"="phenotype_id")) %>% 
  dplyr::distinct(.keep_all = T) %>% 
  dplyr::mutate(qtlGroup = qtlGroup)

print(Sys.time() - start_time)

if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = T)
}

readr::write_tsv(extracted_sumstat, file.path(output_dir, paste0(qtlGroup, ".tsv")))
print("Done")

