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
# unique_variants = dplyr::distinct(variants, variant_id, pos, chr)

# grouping by connected component to query regions rather than each variant separately
# to reduce number of queries
cc_coords = variants %>% 
  dplyr::group_by(cc_id) %>% 
  dplyr::summarise(start = min(pos), end = max(pos), chr=chr[1])
# cc_coords = cc_coords[1:1000, ]

# head(cc_coords)

# extract the name of QTL group
qtlGroup = stringr::str_split(basename(sumstat), pattern, simplify = T)[1]

region <- GenomicRanges::GRanges(
  seqnames = cc_coords$chr,
  ranges = IRanges::IRanges(start = cc_coords$start,
                   end = cc_coords$end)
)

print("Start tabix scan")
start_time <- Sys.time()

regions <- scanTabixDataFrame(sumstat,
  region,
  col_names = eqtl_colnames,
  col_types = readr::cols(alt = "c", ref = "c", type="c", rsid="c", r2="c")
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

