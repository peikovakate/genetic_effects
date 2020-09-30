parser <- optparse::OptionParser()
parser <- optparse::add_option(parser, c('-f', "--freq_file"), 
                               type = "character")
parser <- optparse::add_option(parser, c('-o', '--output'),
                               type="character", 
                               help = "output file")
args = optparse::parse_args(parser)

variants_file = args$freq_file
output_file = args$output
  
variants = readr::read_tsv(variants_file)
variants = dplyr::mutate(variants, variant=paste0('chr', chr, '_', pos, '_', ref, '_', alt))

readr::write_tsv(variants, output_file)
