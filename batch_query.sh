s=/gpfs/hpc/projects/eQTLCatalogue/summary_stats/v0.2/final/

for f in $(ls $s); 
do sbatch query_sumstat_for_variants.sh $s${f}/${f}.nominal.sorted.txt.gz 
done

