#!/bin/bash

#sbatch -o ./outputs/open_targets/antibody/log_ab.out submit_drugnomeai.sh ./outputs/open_targets/antibody ./ot_seed_genes/ab_genes.txt

sbatch -o ./outputs/open_targets/small_molecule/log_sm.out submit_drugnomeai.sh ./outputs/open_targets/small_molecule ./ot_seed_genes/sm_genes.txt

sbatch -o ./outputs/open_targets/oncology_sm/log_onco_sm.out submit_drugnomeai.sh ./outputs/open_targets/oncology_sm ./ot_seed_genes/oncology_sm_all.txt

sbatch -o ./outputs/open_targets/non_oncology_sm/log_non_onco_sm.out submit_drugnomeai.sh ./outputs/open_targets/non_oncology_sm ./ot_seed_genes/non_oncology_sm_all.txt



