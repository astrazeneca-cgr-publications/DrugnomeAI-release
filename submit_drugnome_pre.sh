#!/bin/bash
#SBATCH -o drugnome_pre_pure.out
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:0:0


drugnomeai -c conf/CKD_config.yaml -o out/drugnome_PUREIPTest -i 10 -n 20 -f -r pre -t 1 -x dom
