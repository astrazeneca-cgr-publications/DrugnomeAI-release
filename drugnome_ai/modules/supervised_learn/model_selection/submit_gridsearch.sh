#!/bin/bash 
#SBATCH -J gridsearch 
#SBATCH -o out_gridsearch.txt 
#SBATCH -e err_gridsearch.txt 
#SBATCH --time=24:00:00 
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G


out_dir=$1
# small_search_space: 0 (default - extended parameter search space) or 1 (for debugging: shorter running time)
small_search_space=0 

python -u sklearn_grid_search_cv.py $out_dir $small_search_space
