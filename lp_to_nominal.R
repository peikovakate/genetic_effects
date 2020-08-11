suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(T)

gwas_file = args[1]
out_file = args[2]
n = args[3]

print(args)

gwas = read_tsv(gwas_file, col_names = c("SNP", "A1", "A2", "ES", "SE", "LP"))

gwas = gwas %>%
  rename(B=ES) %>% 
  mutate(P=10^(-LP), N=n)

gwas$P[gwas$P == 0] <- 10^(-3)

gwas %>% select(SNP, A1, A2, P, N, B) %>% write_tsv(out_file)
