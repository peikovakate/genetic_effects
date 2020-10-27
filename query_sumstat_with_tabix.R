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
parser <- optparse::add_option(parser, 
                               c('-t', '--tabix'), 
                               type="character", 
                               default = "tabix",
                               help="path to tabix cmd tool")
args = optparse::parse_args(parser)

sumstat_file <- args$sumstat_file
pattern <- args$pattern
pairs_file <- args$variants
output_dir <- args$output_dir
tabix_path = args$tabix

if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = T)
}

qtl_group = strsplit(basename(sumstat_file), split=pattern)[[1]]
# getting the name of file with coords (CHR, POS)
coords_file = sprintf("%s_coords.tsv", sub('\\.tsv$', '', pairs_file))
temp_sumstat_file = file.path(output_dir, paste0(qtl_group, "_temp.tsv"))

variants = readr::read_tsv(pairs_file, col_types = readr::cols())

# in case of file with variants has inconsitent name (microarray vs rna-seq)
if("phenotype_id" %in% names(variants)){
  variants = dplyr::rename(variants, molecular_trait_id = phenotype_id)  
}
if("variant_id" %in% names(variants)){
  variants = dplyr::rename(variants, variant = variant_id)
}

# set file with variant coordinates as input and redirect output of tabix to temporary file
input = paste(sumstat_file, "-R", coords_file, ">", temp_sumstat_file)
sumstat_out = system2(tabix_path, args=input, stdout = T)

# read results in temp file
sumstat = readr::read_tsv(temp_sumstat_file, col_names = eqtl_colnames, col_types = eqtl_col_types)

sumstat = dplyr::inner_join(sumstat, 
                            variants[c("molecular_trait_id", "variant", "cc_id")], 
                            by=c("variant", "molecular_trait_id"))
sumstat = dplyr::mutate(sumstat, qtl_group = qtl_group)
# remove duplicates with rsid
sumstat = dplyr::select(sumstat, -rsid) %>%  dplyr::distinct(.keep_all = T)
# delete temp file
file.remove(temp_sumstat_file)

print(qtl_group)
print(paste("Found pairs:", nrow(sumstat)))

readr::write_tsv(sumstat, file.path(output_dir, paste0(qtl_group, ".tsv")))

