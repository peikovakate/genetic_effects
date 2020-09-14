source("utils2.R")
`%>%` <- magrittr::`%>%`

parser <- optparse::OptionParser()
parser <- optparse::add_option(parser, c("-e", "--effects"), 
                               type = "character", 
                               help="tsv with effects, output from pick_lead_effects.R script")
parser <- optparse::add_option(parser, c('-o', '--output'),
                               type="character", 
                               help = "output filename for model")
args = optparse::parse_args(parser)

effects_file = args$effects
output_model_path = args$output


effects <- readr::read_tsv(effects_file)
nrow(effects)
eqtls = effects_to_matricies(effects, replace_na_with="zero")

fit_mash = function(eqtls){
  alpha_value = 1
  
  data   = mashr::mash_set_data(eqtls$beta, eqtls$se, alpha=alpha_value)
  m.1by1 = mashr::mash_1by1(data, alpha=alpha_value)
  strong = mashr::get_significant_results(m.1by1)
  U.c    = mashr::cov_canonical(data)
  # n_strong = 1000
  print("Compute PCA covariance matrix")
  U.pca = mashr::cov_pca(data, 5, strong)
  print("Compute extreme deconcolution covariance matrix")
  U.ed = mashr::cov_ed(data, U.pca, strong)
  print("Fit the model")
  # m = mash(data, c(U.c, U.ed))
  m = mashr::mash(data, U.ed, algorithm.version = "Rcpp")
  # save(m, file = "../data/mash/microarr_mash_1_ed_only.R")
  # load("../data/mash/rnaseq_mash_1_ed_only.R")
  return(m)
}

m = fit_mash(eqtls)
sharing = mashr::get_pairwise_sharing(m)

save(m, file = paste0(output_model_path, "_mash_ed.R"))
save(sharing, file = paste0(output_model_path, "_mash_ed_sharing.R"))