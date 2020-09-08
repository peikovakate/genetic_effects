`%>%` <- magrittr::`%>%`

parser <- optparse::OptionParser()
parser <- optparse::add_option(parser, c("-f", "--folder"), 
                               type = "character", 
                               help="folder with variants files")
parser <- optparse::add_option(parser, c('-o', '--output'),
                               type="character", 
                               help = "output filename tsv")
args = optparse::parse_args(parser)

variants_folder = args$folder
output_file = args$output

files = list.files(variants_folder)
tbls = lapply(file.path(variants_folder, files), readr::read_tsv, col_types = readr::cols())
tbls = dplyr::bind_rows(tbls)
head(tbls)

readr::write_tsv(tbls, output_file)
