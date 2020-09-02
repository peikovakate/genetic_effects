# script for querying (tabix) eQTL catalogue sumstat files
# for connected components variants
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rslurm))
suppressPackageStartupMessages(library(optparse))

parser <- OptionParser()
parser <- add_option(parser, c("-c", "--cc_file"), type = "character", help="file with variants and cc_id")
parser <- add_option(parser, c('-j', "--jobname"), type = "character", help="slurm job name")
parser <- add_option(parser, c('-s', "--sumstat_folder"), 
                     default = "/gpfs/hpc/projects/eQTLCatalogue/summary_stats/v0.2/final/",
                     type = "character", help = "folder with eqtl catalogue sumstat")
parser <- add_option(parser, c('-p', '--pattern'),
                     default = "_microarray.nominal.sorted.txt.gz",
                     type="character", 
                     help = "sumstat file pattern ending similar for all tissues")
parser <- add_option(parser, c('-t', '--time'), 
                     default = "4:00:00", 
                     type="character", 
                     help="max job time, ex. 4:00:00")
parser <- add_option(parser, c('-m', '--mem'), 
                     default = "25G", 
                     type="character", 
                     help="max job memory usage, ex. 10G")

args = parse_args(parser)

sumstat <- args$sumstat_folder
pattern <- args$pattern
cc_file <- args$cc_file
jobname <- args$jobname

print(paste("Sumstats folder", sumstat))
print(paste("Sumstats file pattern", pattern))
print(paste("Read file with variants and cc_id from", cc_file))
print(paste("Run slurm jobs with name", jobname))

N_chunks = 10

# getting sumstat file names
files = list.files(sumstat, recursive=T, pattern=paste0(pattern, "$"))

tissues <- sapply(files, strsplit, split = pattern , simplify = T)
tissues <- unlist(tissues)

print("Cell types and tissues:")
print(tissues)

process_job <- function(cell_type, chunk_number, N_chunks, cc_file, sumstat_folder, pattern){
  print(getwd())
  source("../scripts/utils.R")
  get_chunk <- function(chunk_number, connected_components, cc_variants){
    indices <- splitIntoChunks(chunk_number, N_chunks, nrow(connected_components))
    cc_chunk <- connected_components[indices, ]
    cc_variants_chunk <- cc_variants[indices]
    return(list(connected_components = cc_chunk, variants = cc_variants_chunk))
  }
  
  # read file with eQTL-gene pairs and connetected components ids assigned to pairs
  all_pairs <-
    read_tsv(cc_file)
  
  # define minimum and maximum position
  all_pairs <- all_pairs %>%
    group_by(cc_id) %>%
    mutate(
      min_pos = min(pos),
      max_pos = max(pos),
      min_pos_chr = chr[which.min(pos)],
      max_pos_chr = chr[which.max(pos)]
    )
  
  connected_components <- all_pairs %>%
    select(cc_id, phenotype_id, min_pos, max_pos, min_pos_chr) %>%
    distinct(.keep_all = T)
  
  # find variants in each connected component
  cc_groups <- split(all_pairs, all_pairs$cc_id)
  # vector of variants per connected component
  cc_variants <- lapply(cc_groups, "[[", "variant_id")
  
  # 
  chunk_data <- get_chunk(chunk_number, connected_components, cc_variants)
  
  start_time <- Sys.time()
  res <- analyse_chunk(chunk_data$connected_components, chunk_data$variants, cell_type, sumstat_folder, pattern)
  print(Sys.time() - start_time)
  print("Chunk done")
  
  return(res)
}

chunks <- 1:N_chunks
# generate parameters for every job
params = tidyr::crossing(tissues, chunks)
params <- params %>% 
  dplyr::rename(cell_type = tissues, chunk_number = chunks) %>% 
  dplyr::mutate(N_chunks = N_chunks, cc_file = cc_file, sumstat_folder = sumstat, pattern=pattern)

# shuffle rows of parameters to avoid one node processing heavy files
rows <- sample(nrow(params))
params <- params[rows,]

# run a slurm jobs
sjob <- slurm_apply(
  process_job,
  params,
  jobname = jobname,
  submit = T,
  nodes = 25,
  cpus_per_node = 1,
  slurm_options = list(time = args$time, mem=args$mem)
)


