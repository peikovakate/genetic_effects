################
## Gather the results
################
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rslurm))
suppressPackageStartupMessages(library(optparse))

parser <- OptionParser()
parser <- add_option(parser, c("-o", "--output"), type = "character", help="output filename for found variants in sumstats")
parser <- add_option(parser, c('-j', "--jobname"), type = "character", help="slurm job name")

args = parse_args(parser)
print(args)

sjob = slurm_job(args$jobname, 25)
print("Loading batch results")
res = get_slurm_out(sjob, outtype="raw")
found_variants = bind_rows(res)

print("Writing to file")
write_tsv(found_variants, args$output, quote_escape=F)