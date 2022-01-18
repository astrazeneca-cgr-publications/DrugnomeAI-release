#!/bin/bash
#SBATCH -o drugnome_puNB_agnostic.out
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:0:0


drugnomeai -c mantis_ml/conf/IPF_config.yaml -o out/IPF-custom-seeds-no_HLA -i 5 -n 20 -k custom-data/IPF-seed-genes/IPF.custom_seeds-no_HLA.txt --fast


drugnomeai -c mantis_ml/conf/IPF_config.yaml -o out/IPF-HPO-seeds-no_HLA -i 5 -n 20 -k custom-data/IPF-seed-genes/IPF.HPO_seeds-no_HLA.txt --fast


drugnomeai -c mantis_ml/conf/Covid19.yaml -o out/Covid_19 -i 5 -n 20 -k custom-data/Covid19/seed_genes.txt --fast


drugnomeai -c conf/CKD_config.yaml -o out/drugnome_agnosticT1 -i 10 -n 20 -s nb -r pu -t 1
