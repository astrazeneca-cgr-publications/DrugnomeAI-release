#!/bin/bash
#SBATCH -o ckd-XGBoost.out
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:0:0


# Various Test runs
iterations=10

python __main__.py -c ../../conf/CKD_config.yaml -o ../../../../../../../projects/cgr/users/cgr-ds/kvtp611/CKD-XGBoost-debug -n 10 -i $iterations -m xgb


#python __main__.py -c ../../conf/CKD_config.yaml -o CKD-test-data -n 10 -i $iterations -r post -f


#python __main__.py -c ../../conf/CKD_config.yaml -o CKD-et_rf -n 10 -i $iterations -m et,rf


#python __main__.py -c ../../conf/CKD_config.yaml -o CKD-stacking -n 10 -i $iterations -m stack

