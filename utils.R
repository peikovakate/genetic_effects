#' A general function to quickly import tabix indexed tab-separated files into data_frame
#'
#' @param tabix_file Path to tabix-indexed text file
#' @param param An instance of GRanges, RangedData, or RangesList
#' provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()
#' @param ... Additional parameters to be passed on to readr::read_delim()
#'
#' @return List of data_frames, one for each entry in the param GRanges object.
#' @export
scanTabixDataFrame <- function(tabix_file, param, ...) {
  tabix_list = Rsamtools::scanTabix(tabix_file, param = param)
  df_list = lapply(tabix_list, function(x, ...) {
    if (length(x) > 0) {
      if (length(x) == 1) {
        #Hack to make sure that it also works for data frames with only one row
        #Adds an empty row and then removes it
        result = paste(paste(x, collapse = "\n"), "\n", sep = "")
        result = readr::read_delim(result, delim = "\t", ...)[1, ]
      } else{
        result = paste(x, collapse = "\n")
        result = readr::read_delim(result, delim = "\t", ...)
      }
    } else{
      # Return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  return(df_list)
}

# colomn names for old sumstat with no se
# eqtl_colnames = c("variant","type","rsid","ref","r2","pvalue","position","molecular_trait_object_id","molecular_trait_id","median_tpm","maf","gene_id","chromosome","beta","an","alt","ac")

# microarr columns
# eqtl_colnames = c(
#   "molecular_trait_id",
#   "chromosome",
#   "position",
#   "ref",
#   "alt",
#   "variant",
#   "ma_samples",
#   "ma_count",
#   "maf",
#   "pvalue",
#   "beta",
#   "se"
# )

# rnaseq columns
eqtl_colnames = c(
  "molecular_trait_id",
  "chromosome",
  "position",
  "ref",
  "alt",
  "variant",
  "ma_samples",
  "ac",
  "an",
  "maf",
  "pvalue",
  "beta",
  "se",
  "molecular_trait_object_id",
  "gene_id",
  "mediat_tpm",
  "r2",
  "type",
  "rsid"
)

splitIntoBatches <- function(n, batch_size) {
  n_batches = ceiling(n / batch_size)
  batch_ids = rep(seq(1:n_batches), each = batch_size)[1:n]
  return(batch_ids)
}

splitIntoChunks <- function(chunk_number, n_chunks, n_total) {
  chunk_size = floor(n_total / (n_chunks))
  batches = splitIntoBatches(n_total, chunk_size)
  batches[batches > n_chunks] = n_chunks
  selected_batch = batches == chunk_number
  return(selected_batch)
}

analyse_chunk <-
  function(chunk,
           chunk_variants,
           qtlGroup,
           sumstat_folder,
           pattern) {
    region <-  GenomicRanges::GRanges(
      seqnames = chunk$min_pos_chr,
      ranges = IRanges::IRanges(start = chunk$min_pos,
                       end = chunk$max_pos)
    )
    
    regions <- scanTabixDataFrame(
      sprintf(paste0(sumstat_folder, "%s", pattern),
              qtlGroup),
      region,
      col_names = eqtl_colnames,
      col_types = readr::cols(alt = "c", ref = "c", type="c", rsid="c")
    )
    
    # name regions after connected component id
    names(regions) <- chunk$cc_id
    
    # in all eQTLs found in region, 
    # filter by phenotype and the varaints of connected component
    eqtls <-
      mapply(function(r, target_variants, phenotype) {
        if (!is.null(r)) {
          dplyr::filter(r, variant %in% target_variants & molecular_trait_id == phenotype)
        }
      },
      regions,
      chunk_variants,
      chunk$phenotype_id,
      SIMPLIFY = F)
    
    eqtls_mapped <- mapply(function(eqtl_region, cc_id_name) {
      if (!is.null(eqtl_region)) {
        dplyr::mutate(eqtl_region, cc_id = cc_id_name)
      }
    }, eqtls, names(eqtls), SIMPLIFY = F)
    
    eqtls_mapped <- dplyr::bind_rows(eqtls_mapped)
    # eqtls_mapped <- dplyr::mutate(eqtls_mapped, qtlGroup = qtlGroup)
    return(eqtls_mapped)
  }
