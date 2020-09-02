for chr in {1..22};
do echo ${chr}; 
Rscript factors_to_bed_continuous.R \
  -b  ../data/bim_files/eur_chr${chr}_biall_maf_cm.bim \
  -f ../data/rnaseq/effects_with_na/cc_variants_with_na_mapped/factor_binary_one_per_eqtl.tsv\
  -o ../data/rnaseq/effects_with_na/cc_variants_with_na_mapped/factor_binary_one_per_eqtl/${chr};
done;