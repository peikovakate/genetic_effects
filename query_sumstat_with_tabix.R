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

qtl_group = strsplit(basename(sumstat_file), split=pattern)[[1]]
coords_file = sprintf("%s_coords.tsv", sub('\\.tsv$', '', pairs_file))

variants = readr::read_tsv(pairs_file)
if("phenotype_id" %in% names(variants)){
  variants = dplyr::rename(variants, molecular_trait_id = phenotype_id)  
}
if("variant_id" %in% names(variants)){
  variants = dplyr::rename(variants, variant = variant_id)
}

coords = variants %>% dplyr::select(chr, pos) %>% 
  dplyr::rename(CHROM = chr, POS=pos) %>% 
  dplyr::distinct(.keep_all = T) %>% 
  dplyr::arrange(CHROM, POS)

readr::write_tsv(coords, coords_file, col_names = F)

input = paste(sumstat_file, "-R", coords_file)
sumstat_out = system2(tabix_path, args=input, stdout = T)

handle_sumstat_output = function(output){
  if (length(output) > 0) {
    if (length(output) == 1) {
      #Hack to make sure that it also works for data frames with only one row
      #Adds an empty row and then removes it
      result = paste(paste(sumstat_out, collapse = "\n"), "\n", sep = "")
      result = readr::read_tsv(sumstat_out, col_names = eqtl_colnames, col_types = eqtl_col_types)[1,]
    } else{
      result = paste(sumstat_out, collapse = "\n")
      result = readr::read_tsv(sumstat_out, col_names = eqtl_colnames, col_types = eqtl_col_types)
    }
  } else{
    # Return NULL if the nothing is returned from tabix file
    result = NULL
  }
  return(result)
}


sumstat = handle_sumstat_output(sumstat_out)
sumstat = dplyr::inner_join(sumstat, 
                            variants[c("molecular_trait_id", "variant", "cc_id")], 
                            by=c("variant", "molecular_trait_id"))
sumstat = dplyr::mutate(sumstat, qtl_group = qtl_group)
# remove duplicates with rsid
sumstat = dplyr::select(sumstat, -rsid) %>%  dplyr::distinct(.keep_all = T)

print(qtl_group)
print(paste("Found pairs:", nrow(sumstat)))

if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = T)
}

readr::write_tsv(sumstat, paste0(output_dir, qtl_group, ".tsv"))

