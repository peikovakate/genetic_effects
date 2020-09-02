x = read_tsv("../data/rnaseq/unique_genes/slope.txt")
names = paste(x$Gene, x$SNP)
length(unique(names))
x %>% distinct() %>% nrow()
