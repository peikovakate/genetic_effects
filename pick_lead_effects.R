`%>%` <- magrittr::`%>%`
source("utils2.R")

parser <- optparse::OptionParser()
parser <- optparse::add_option(parser, c("-s", "--sumstat"), 
                               type = "character", 
                               default = "../data/gtex/sumstat_comb.tsv",
                               help="sumstat data for all qtl groups combined")
parser <- optparse::add_option(parser, c('-o', '--output'), 
                               default = "../data/gtex/",
                               type="character", 
                               help="output folder")
parser <- optparse::add_option(parser, c('-f', '--frac'), 
                               default = 0.95, 
                               type="double", 
                               help="fraction of qtl groups (present in at least fraction) for subset of effects")

args = optparse::parse_args(parser)
print(args)

sumstat_file = args$sumstat
output_folder = args$output
frac = args$frac

if(!dir.exists(output_folder)){
  dir.create(output_folder, recursive = T)
}

sumstat = readr::read_tsv(sumstat_file, 
                        col_types=readr::cols_only(variant='c', molecular_trait_id='c', chromosome='c', position='c',
                                                   beta='d', pvalue='d', se='d', cc_id='c', qtl_group='c'))

sumstat = dplyr::mutate(sumstat, eqtl_id = paste(variant, molecular_trait_id, sep="."))
sumstat = dplyr::distinct(sumstat, eqtl_id, qtl_group, .keep_all = T)

ccs = dplyr::distinct(sumstat, eqtl_id, cc_id, variant, molecular_trait_id, chromosome, position)

data = sumstat_to_effects(sumstat)

effects = dplyr::left_join(data$beta, data$pvalue, by="eqtl_id", suffix=c(".beta", ".pvalue"))
names(data$se) = paste0(names(data$se), ".se")
effects = dplyr::left_join(effects, data$se, by=c("eqtl_id"="eqtl_id.se"))

count_occurrence = function(data){
  counts = data %>% 
    dplyr::select(dplyr::ends_with('.beta')) %>% 
    is.na %>% 
    `!` %>% 
    rowSums()
  data = dplyr::mutate(data, qtl_group_count = counts)
  return(data)
}

effects = count_occurrence(effects)
min_pvalue = apply(dplyr::select(effects, dplyr::ends_with(".pvalue")), 1, min, na.rm=T) 
effects = dplyr::mutate(effects, min_pvalue = min_pvalue)

effects = dplyr::left_join(effects, ccs, by="eqtl_id")

lead_effects = effects %>% 
  dplyr::group_by(cc_id) %>% 
  dplyr::filter(qtl_group_count == max(qtl_group_count)) %>% 
  dplyr::filter(min_pvalue == min(min_pvalue, na.rm=T)) %>% 
  dplyr::sample_n(1) %>% 
  dplyr::ungroup()

nrow(lead_effects)

readr::write_tsv(lead_effects, file.path(output_folder, "lead_effects_na.tsv"))

lead_pairs = dplyr::distinct(lead_effects, variant, molecular_trait_id, position, chromosome)
readr::write_tsv(dplyr::select(lead_pairs, molecular_trait_id, variant, chromosome, position), 
                 file.path(output_folder, "lead_pairs.tsv"))

max_n = max(lead_effects$qtl_group_count)
subset_lead_effects = lead_effects %>% 
  dplyr::filter(qtl_group_count >= max_n*frac)

nrow(subset_lead_effects)

readr::write_tsv(subset_lead_effects, file.path(output_folder, "subset_lead_effects_na.tsv"))

                      