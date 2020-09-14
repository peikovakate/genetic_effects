`%>%` <- magrittr::`%>%`

output_file = "variants_across_tissues_rnaseq_no_effect_to_zero_fract0.tsv"
sumstat_variants = "rnaseq_variants.tsv"
is_mean = F
fraction = 0.95

# file with sumstat of variants from connected components
# found_variants = read_tsv("data/found_in_sumstat2_microarr.tsv")
found_variants = readr::read_tsv(sumstat_variants)
unique_tissues = unlist(distinct(found_variants, cell_type))
n_unique_tissues = length(unique_tissues)

print(paste("File with sumstat data contains", nrow(found_variants), "records"))

found_variants = found_variants %>% distinct(.keep_all = T)
print(paste("There are", nrow(found_variants), "unique records"))

# for every variant in some cc count number of tissues and cell types where it's found
found_variants <- found_variants %>%
  dplyr::group_by(cc_id, variant) %>%
  dplyr::mutate(count = length(unique(cell_type))) %>% 
  dplyr::ungroup()

print("Variants per number of cell types or tissues (where it's found)")
print(table(found_variants$count))

connected_components <- found_variants %>%
  dplyr::group_by(cc_id) %>%
  dplyr::group_split()

print("Number of cc")
print(length(connected_components))

most_freq_across_component <-
  lapply(connected_components, function(cc) {
    cc <- ungroup(cc)
    # chose one variant from cc with smallest pvalue
    target_eqtl <- dplyr::filter(cc, count == max(count)) %>%
      dplyr::filter(pvalue == min(pvalue)) %>% 
      sample_n(1)
  
    if (target_eqtl$count == n_unique_tissues) {
      dplyr::filter(cc, variant == target_eqtl$variant) %>% distinct(cell_type, .keep_all = T)
    } else if (target_eqtl$count >= floor(n_unique_tissues * fraction)) {
      found = dplyr::filter(cc, variant == target_eqtl$variant) %>% distinct(cell_type, .keep_all = T)
      tissues_not_found = dplyr::setdiff(unique_tissues, unique((found$cell_type)))
      rest <- lapply(tissues_not_found, function(x) {
        found[1,] %>% mutate(cell_type = x)
      })
      rest <- bind_rows(rest)
      if(is_mean){
        rest <-
          mutate(
            rest,
            pvalue = mean(found$pvalue),
            beta = mean(found$beta),
            se = mean(found$se)
          )        
      }else{
        rest <-
          mutate(
            rest,
            pvalue = NA,
            beta = NA,
            se = NA
          )
      }
      bind_rows(list(found, rest))
    }
  })

all_variants <- dplyr::bind_rows(most_freq_across_component)
all_variants %>% dplyr::select(cc_id) %>% dplyr::distinct() %>% nrow()

# removing directory name from cell type names
# cell_types = strsplit(all_variants$cell_type, split = "/")
# cell_types = unlist(lapply(cell_types, "[[", 2))
# all_variants = mutate(all_variants, cell_type = cell_types)

readr::write_tsv(all_variants, output_file, quote_escape = F)
