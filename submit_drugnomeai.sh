#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=48:0:0
##SBATCH -o drugnomeai.out

#module load Anaconda3/5.3.0
#source activate drugnome_ai

out=$1
seed_genes_file=$2
time drugnomeai -o $out -k $seed_genes_file -l -r boruta
