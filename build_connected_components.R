`%>%` <- magrittr::`%>%`

parser <- optparse::OptionParser()
parser <- optparse::add_option(parser, c("-m", "--fine_mapping"), 
                               type = "character", 
                               default = "../data/susie_rnaseq/",
                               help="folder with fine-mapping results")
parser <- optparse::add_option(parser, c('-p', '--pattern'),
                               default = ".ge.txt.gz",
                               type ="character", 
                               help = "fine-mapping file pattern ending, same for all qtl groups")
parser <- optparse::add_option(parser, c('-o', '--output'), 
                               default = "../data2/connected_components/rnaseq_cc.tsv",
                               type="character", 
                               help="output filename")
parser <- optparse::add_option(parser, c('-z', '--z_score'), 
                               default = 3, 
                               type="double", 
                               help="z score threshold, (absolute z values are compared)")
parser <- optparse::add_option(parser, c('-c', "--cs_size"), 
                               default = 50, 
                               type="integer", 
                               help="max size of credible set to include into connected components")

args = optparse::parse_args(parser)
print(args)

susie_path = args$fine_mapping
suffix = args$pattern
output_file = args$output
z_threshold <- args$z_score
cs_size_threshold <- args$cs_size

# print("Files with given pattern")
# print(files)

#' read_fine_mapping
#'
#' @param file_names vector of file names in the folder to read,
#' files should have cs_id unique across genes
#' @param susie_folder_path folder path with files
#'
#' @return tibble combining all files with additional two columns:
#' qtl_group (name of the cell type + study) and cs_uid (credibe set id that includes qtl_group)
read_fine_mapping <- function(file_names, susie_folder_path){
  print("Loading fine-mapped files")
  tbls <- lapply(file_names, function(x){
    file = file.path(susie_folder_path, x)
    qtl_group_name = strsplit(x, split=suffix)[[1]]
    print(sprintf("%s: %s", qtl_group_name, file))
    tbl <- readr::read_tsv(file, progress=F, col_types = readr::cols()) %>% dplyr::mutate(qtl_group = qtl_group_name)
    return(tbl)
  })
  tbl_fine_mapped <- dplyr::bind_rows(tbls)
  # make unique id for credible sets
  # if cs_id is unique, cs_id and cell type is enough
  tbl_fine_mapped <- dplyr::mutate(tbl_fine_mapped, cs_uid = paste(cs_id, qtl_group, sep="_"))
  qtl_groups = unique(tbl_fine_mapped$qtl_group)
  return(list(credible_sets = tbl_fine_mapped, qtl_groups = qtl_groups))
}

#' cs_counts 
#'
#' @param fine_mapped tibble with qtl group, cs_uid (credible set unique id)
#'
#' @return list with counts of credible sets and genes
cs_counts <-  function(fine_mapped){
  cs = dplyr::group_by(tbl_fine_mapped, qtl_group) %>% 
    dplyr::summarise(cs = dplyr::n_distinct(cs_id))
  genes = dplyr::group_by(tbl_fine_mapped, qtl_group) %>% 
    dplyr::summarise(genes = dplyr::n_distinct(phenotype_id))
  return(list(credible_sets = cs, genes = genes))
}

csDfToList <- function(cs_df){
  grouped_df = dplyr::group_by(cs_df, cs_uid)
  cs_list = setNames(dplyr::group_split(grouped_df), dplyr::group_keys(grouped_df)$cs_uid) %>%
    purrr::map(~.$variant_id)
  return(cs_list)
}

get_cc <- function(phenotype_group){
  # phenotype_group <-  ungroup(phenotype_group)
  variants <- csDfToList(phenotype_group)
  if(length(variants) < 2){
    combinations <- matrix(c(1,1), nrow=2)
  }else{
    combinations <- combn(length(variants), 2, simplify = T)
  }
  
  # combinations[1,]
  intersections <- purrr::map2(combinations[1,], combinations[2,], function(x, y){
    intr <- intersect(variants[[x]], variants[[y]])
    return(length(intr)>0)
  })
  
  n_verts <-  length(variants)
  self_edges <- matrix(c(1:n_verts, 1:n_verts), nrow=n_verts)
  combs <- t(combinations[, unlist(intersections)])
  edge_list <- rbind(self_edges, combs)
  
  g <- igraph::graph_from_edgelist(edge_list)
  # g <- as.undirected(g)
  
  components_sets <- 1:igraph::components(g)$no %>% 
    lapply(function(x){return(names(variants)[igraph::components(g)$membership == x])})
  
  res <- lapply(components_sets, function(x){
    phenotype_group %>% 
      dplyr::filter(cs_uid %in% x) %>% 
      dplyr::mutate(cc_id = paste(x, collapse="_"))
  }) %>% dplyr::bind_rows()
  
  return(res)
  
}

#' build_connected_components
#'
#' @param credible_sets tibble with credible set variants, cs_uid
#'
#' @return same tibble with additional column cc_id (id of connected component)
build_connected_components <- function(credible_sets){
  gene_groups <- credible_sets %>% 
    dplyr::group_by(phenotype_id) %>% 
    dplyr::group_split()
  
  print(paste("Unique genes", length(gene_groups)))
  
  variants_connected_components <- lapply(gene_groups, get_cc)
  variants_connected_components <- dplyr::bind_rows(variants_connected_components)
  
  print("Total connected components:")
  print(nrow(dplyr::distinct(variants_connected_components, cc_id)))
  return(variants_connected_components)
}


paths = unlist(strsplit(susie_path, " "))
tbl_fine_mapped = lapply(paths, function(fine_mapping_path){
  files <- list.files(fine_mapping_path, pattern=suffix)
  # reading the files with fine-mapped credible sets from susie
  fine_mapped = read_fine_mapping(files, susie_folder_path = fine_mapping_path)
  # fine_mapped$qtl_groups
  return(fine_mapped$credible_sets)
})


tbl_fine_mapped = dplyr::bind_rows(tbl_fine_mapped)

# find number of unique genes and credible sets per qtl group
cs_counts(tbl_fine_mapped) 

# filter the credible sets by its size and z-score
filtered_cs <- tbl_fine_mapped %>% 
  dplyr::group_by(cs_uid) %>% 
  dplyr::mutate(max_abs_z = max(abs(z))) %>% 
  dplyr::filter(max_abs_z > z_threshold, cs_size < cs_size_threshold) %>% 
  dplyr::ungroup()

# find connected components among credible sets
variants_connected_components <- build_connected_components(filtered_cs)

# save fine-mapped variants with connected component id to file
directory_name = dirname(output_file)
if(!dir.exists(directory_name)){
  dir.create(directory_name, recursive = T)
}

print(paste("writing to file", output_file))
readr::write_tsv(variants_connected_components, output_file, quote_escape=F)

coords_file = sprintf("%s_coords.tsv", sub('\\.tsv$', '', output_file))

print(paste("writing variant coords to file", coords_file))

coords = variants_connected_components %>% dplyr::select(chr, pos) %>% 
  dplyr::rename(CHROM = chr, POS=pos) %>% 
  dplyr::distinct(.keep_all = T) %>% 
  dplyr::arrange(CHROM, POS)
readr::write_tsv(coords, coords_file, col_names = F)

