variants_folder = "../data2/sumstat/"
output_file = "../data2/sumstat_from_cc/cc_eqtls.tsv"

files = list.files(variants_folder)
tbls = lapply(file.path(variants_folder, files), readr::read_tsv)
tbls = dplyr::bind_rows(tbls)
head(tbls)

readr::write_tsv(tbls, output_file)
