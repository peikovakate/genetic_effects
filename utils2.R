`%>%` <- magrittr::`%>%`

effects_to_matricies = function(effects, replace_na_with = FALSE){
  effects_matrix <- dplyr::select(effects, ends_with('.beta')) %>% 
    dplyr::rename_all(function(x) {sub(".beta", "", x)})
  errors_matrix <- dplyr::select(effects, ends_with('.se')) %>% 
    dplyr::rename_all(function(x){sub(".se", "", x)})
  
  effects_matrix <- as.matrix(effects_matrix)
  errors_matrix <- as.matrix(errors_matrix)
  
  missing_values <- which(is.na(effects_matrix), arr.ind=TRUE)
  if(replace_na_with == "mean"){
    effects_matrix[missing_values] <- rowMeans(effects_matrix, na.rm=TRUE)[missing_values[,1]]
    errors_matrix[missing_values] <- rowMeans(errors_matrix, na.rm=TRUE)[missing_values[,1]]
  }else if(replace_na_with == "zero"){
    effects_matrix[missing_values] <-  0
    errors_matrix[missing_values] <- 1
  }
  
  # effects_matrix <- effects_matrix*sign(effects_matrix[max.col(abs(effects_matrix))])
  
  # multiply by the sign of the strongest effect across tissues
  # effects_matrix <- t(apply(effects_matrix, 1, function(x){
  #   abs_x <- abs(x)
  #   x*sign(x[which.max(abs_x)])
  # }))
  
  return(list(beta=effects_matrix, se=errors_matrix))
}


sumstat_to_effects <- function(sumstat){
  # sumstat = dplyr::arrange(sumstat, eqtl_id)
  pvalues = dplyr::as_tibble(reshape2::dcast(sumstat, eqtl_id ~ qtl_group, value.var = "pvalue", fill = NA))
  betas = dplyr::as_tibble(reshape2::dcast(sumstat, eqtl_id ~ qtl_group, value.var = "beta", fill = NA))
  ses = dplyr::as_tibble(reshape2::dcast(sumstat, eqtl_id ~ qtl_group, value.var = "se", fill = NA))
  return(list(pvalue = pvalues, beta = betas, se = ses))
}

loadings_to_tibble <- function(loadings, factor_names = NULL){
  colnames(loadings) = paste0("Factor", 1:ncol(loadings))
  genes = sapply(strsplit(rownames(loadings), split = " "), "[[", 1)
  variants = sapply(strsplit(rownames(loadings), split = " "), "[[", 2)
  eqtl_ids = paste(variants, genes, sep=".")
  loadings$eqtl_id = eqtl_ids
  loadings$variant = variants
  loadings$gene = genes
  return(dplyr::as_tibble(loadings))
}

