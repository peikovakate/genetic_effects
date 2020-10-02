#!/bin/bash
#SBATCH -p main
#SBATCH -t 00:30:00
#SBATCH --mem=10G
#SBATCH --output=query_sumstat_with_tabix.out

module load bcftools-1.9

s=$1
p=$2
v=$3
o=$4

Rscript query_sumstat_with_tabix.R -v $v -p $p -s $s -o $o

