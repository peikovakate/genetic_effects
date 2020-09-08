# genetic_effects

Workflow:

1. Build connected components from fine-mapped credible sets
2. Fetch for variants from connected components in eQTL sumstat files for each qtl group

```
s=/gpfs/hpc/projects/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out/sumstats/
p=_ge.nominal.sorted.tsv.gz
v=../data/cc_rnaseq.tsv
o=../data/sumstat_rnaseq/

for f in $(ls $s/*$p); 
do sbatch query_sumstat_with_tabix.sh $f $p $v $o;
done
```

3. Combine eqtl susmstats across all qtl groups (merge tsvs into one)

```
srun --mem=20G -t 30 Rscript combine_query_sumstat_for_variants.R -f ../data/sumstat_rnaseq/ -o ../data/sumstat_rnaseq_comb.tsv
```

4. Select eqtls that are most frequent (or at least in 95% of qtlGroups)
5. Extract effects (beta, se, p-value) into sep tibble
6. Find sharing and similarities