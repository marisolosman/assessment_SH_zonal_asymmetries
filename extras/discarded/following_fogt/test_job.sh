#!/bin/bash
 
#SBATCH --ntasks=1 #(default, so here optional)
#SBATCH --array=1-1000
#SBATCH --cpus-per-task=1  #(default, so here optional)
#SBATCH --job-name=correlations_ninio_antarctica
#SBATCH --output=out_correlations.out 
#SBATCH --partition=cluster #(optional, default is cluster)
#SBATCH --time=120:00 #(optional, default is 24 hours)
#SBATCH --mem-per-cpu=1G #(optional, default is 1024 MB)
iter=`expr 9000 + ${SLURM_ARRAY_TASK_ID}`
echo $iter
python regress_ninio_antarctica.py ${iter}̣
