# script for querying (tabix) eQTL catalogue sumstat files
# for connected components variants
`%>%` <- magrittr::`%>%`


parser <- optparse::OptionParser()
parser <- optparse::add_option(parser, 
                               c("-c", "--cc_file"), 
                               type = "character", 
                               help="file with variants and cc_id (one directory upper)")
parser <- optparse::add_option(parser, 
                               c('-j', "--jobname"), 
                               type = "character", 
                               help="slurm job name")
parser <- optparse::add_option(parser, 
                               c('-s', "--sumstat_folder"), 
                               default = "/gpfs/hpc/projects/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out/sumstats/",
                               type = "character", help = "folder with eqtl catalogue sumstat")
parser <- optparse::add_option(parser, c('-p', '--pattern'),
                               default = "_ge.nominal.sorted.tsv.gz",
                               type="character", 
                               help = "sumstat file pattern ending similar for all tissues")
parser <- optparse::add_option(parser, c('-t', '--time'), 
                               default = "4:00:00", 
                               type="character", 
                               help="max job time, ex. 4:00:00")
parser <- optparse::add_option(parser, c('-m', '--mem'), 
                               default = "25G", 
                               type="character", 
                               help="max job memory usage, ex. 10G")

args = optparse::parse_args(parser)

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

qtlGroups <- sapply(files, strsplit, split = pattern , simplify = T)
qtlGroups <- unlist(qtlGroups)

print("Cell types and tissues:")
print(qtlGroups)

process_job <- function(qtlGroup, chunk_number, N_chunks, cc_file, sumstat_folder, pattern){
  print(getwd())
  # the script is executed from directory _rslurm_[jobname]
  # so we need to source utils from one level higher
  source("../utils.R")
  `%>%` <- magrittr::`%>%`
  
  get_chunk <- function(chunk_number, connected_components, cc_variants){
    indices <- splitIntoChunks(chunk_number, N_chunks, nrow(connected_components))
    cc_chunk <- connected_components[indices, ]
    cc_variants_chunk <- cc_variants[indices]
    return(list(connected_components = cc_chunk, variants = cc_variants_chunk))
  }
  
  # read file with eQTL-gene pairs and connetected components ids assigned to pairs
  all_pairs <- readr::read_tsv(cc_file, col_types = readr::cols(alt = "c", ref = "c"))
  
  # define minimum and maximum position
  all_pairs <- all_pairs %>%
    dplyr::group_by(cc_id) %>%
    dplyr::mutate(
      min_pos = min(pos),
      max_pos = max(pos),
      min_pos_chr = chr[which.min(pos)],
      max_pos_chr = chr[which.max(pos)]
    )
  
  connected_components <- all_pairs %>%
    dplyr::select(cc_id, phenotype_id, min_pos, max_pos, min_pos_chr) %>%
    dplyr::distinct(.keep_all = T)
  
  # find variants in each connected component
  cc_groups <- split(all_pairs, all_pairs$cc_id)
  # vector of variants per connected component
  cc_variants <- lapply(cc_groups, "[[", "variant_id")
  
  # 
  chunk_data <- get_chunk(chunk_number, connected_components, cc_variants)
  
  start_time <- Sys.time()
  res <- analyse_chunk(chunk_data$connected_components, chunk_data$variants, qtlGroup, sumstat_folder, pattern)
  print(Sys.time() - start_time)
  print("Chunk done")
  
  return(res)
}

chunks <- 1:N_chunks
# generate parameters for every job
params = tidyr::crossing(qtlGroups, chunks)
params <- params %>% 
  dplyr::rename(qtlGroup = qtlGroups, chunk_number = chunks) %>% 
  dplyr::mutate(N_chunks = N_chunks, cc_file = cc_file, sumstat_folder = sumstat, pattern=pattern)

# shuffle rows of parameters to avoid one node processing heavy files
rows <- sample(nrow(params))
params <- params[rows,]

# run a slurm jobs
sjob <- rslurm::slurm_apply(
  process_job,
  params,
  jobname = jobname,
  submit = T,
  nodes = 25,
  cpus_per_node = 1,
  slurm_options = list(time = args$time, mem=args$mem)
)


