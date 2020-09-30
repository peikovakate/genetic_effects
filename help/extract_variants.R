parser <- optparse::OptionParser()
parser <- optparse::add_option(parser, c("-v", "--variants"), 
                               type = "character")
parser <- optparse::add_option(parser, c('-f', "--freq_file"), 
                               type = "character")
parser <- optparse::add_option(parser, c('-o', '--output'),
                               type="character", 
                               help = "output file")
args = optparse::parse_args(parser)

variants = readr::read_tsv(args$variants)
freq = readr::read_tsv(args$freq_file)

res = dplyr::inner_join(variants, freq, by="variant")

readr::write_tsv(res, args$output)

