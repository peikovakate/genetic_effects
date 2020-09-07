#!/bin/bash
#SBATCH -p main
#SBATCH -t 2:00:00
#SBATCH --mem=50G
#SBATCH --output=query_sumstat_for_variants.out

s=$1
p=$2
v=$3
o=$4

Rscript query_sumstat_for_variants.R -v $v -p $p -s $s -o $o
	
