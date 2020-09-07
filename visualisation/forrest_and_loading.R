`%>%` <- magrittr::`%>%`

gwas = "RA-ebi-a-GCST002318"

eqtls = readr::read_tsv("../data2/additional/novel_sign_eqtlcat_colocs.tsv")
eqtls = dplyr::filter(eqtls, gwas_id == gwas)

print("Records")
eqtls %>% nrow()

print("Unique pairs")
eqtls %>% dplyr::distinct(variant, molecular_trait_id) %>% nrow()

