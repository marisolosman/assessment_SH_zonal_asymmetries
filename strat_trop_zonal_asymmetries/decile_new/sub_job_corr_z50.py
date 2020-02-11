#!/usr/bin/env python

import os
import numpy as np

def mkdir_p(dir):
    '''make a directory (dir) if it doesn't exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

job_directory = "%s/.job" %os.getcwd()

mkdir_p(job_directory)

for i in np.arange(8, 11):
	job_file = os.path.join(job_directory, "task_%s.job" %i)
	with open(job_file, 'w') as fh:
		fh.writelines("#!/bin/bash\n")
		fh.writelines("#SBATCH --ntasks=1 #(default, so here optional)\n")
		fh.writelines("#SBATCH --array=1-1000\n")
		fh.writelines("#SBATCH --cpus-per-task=1  #(default, so here optional)\n")
		fh.writelines("#SBATCH --job-name=correlations_asymmetries_%s\n" %str(i))
		fh.writelines("#SBATCH --output=out_correlations_%s.out\n" %str(i))
		fh.writelines("#SBATCH --partition=cluster #(optional, default is cluster\n")
		fh.writelines("#SBATCH --time=10:00 #(optional, default is 24 hours)\n")
		fh.writelines("#SBATCH --mem-per-cpu=2048MB #(optional, default is 1024 MB)\n")
		fh.writelines("iter=`expr %s + ${SLURM_ARRAY_TASK_ID}`\n" %str((i-1)*1000))
		fh.writelines("python correlations_asymmetries_z50_new.py ${iter}̣\n")

	os.system("sbatch %s" %job_file)

