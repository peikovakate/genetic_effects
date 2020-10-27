# genetic_effects

Workflow:

1. Build connected components from fine-mapped credible sets

```
srun --mem=20G -t 40 Rscript build_connected_components.R -m "/gpfs/hpc/projects/eQTLCatalogue/susie-finemapping/GTEx_v7/susie /gpfs/hpc/projects/eQTLCatalogue/susie-finemapping/eQTL_Catalogue_r3/susie" -p _ge.purity_filtered.txt.gz -o ../data/gtex/cc.tsv
```

2. Fetch for variants from connected components in eQTL sumstat files for each qtl group

```
s=/gpfs/hpc/projects/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out/sumstats/
s=/gpfs/hpc/projects/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/GTEx_v7/sumstats/
p=_ge.nominal.sorted.tsv.gz
v=../data/gtex/cc.tsv
o=../data/gtex/sumstat

for f in $(ls $s/*$p); 
do sbatch query_sumstat_with_tabix.sh $f $p $v $o;
done
```

3. Combine eqtl susmstats across all qtl groups (merge tsvs into one)

```
srun --mem=40G -t 50 Rscript combine_query_sumstat_for_variants.R -f ../data/gtex/sumstat -o ../data/gtex/sumstat_comb.tsv
``` 

4. Select lead eqtls

```
srun --mem=20G -t 50 Rscript pick_lead_effects.R -s ../data/microarr/sumstat_comb.tsv -o ../data/microarr/
```

5. Extract effects (beta, se, p-value) into sep tibble
6. Find sharing and similarities