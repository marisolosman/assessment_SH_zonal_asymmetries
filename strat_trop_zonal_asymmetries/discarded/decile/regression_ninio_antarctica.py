#!/usr/bin/python3

#SBATCH --job-name=correlations_ninio_antarctica
#SBATCH --output=out_correlations.txt 
#SBATCH --partition=cluster #(optional, default is cluster)
#SBATCH -c 10 #10 cores
#SBATCH --time=120:00 #(optional, default is 24 hours)
#SBATCH --mem-per-cpu=4GB #(optional, default is 1024 MB)

import os
import numpy as np
import pandas as pd
import xarray as xr
import time
from random import randint
import multiprocessing as mp
from pathos.multiprocessing import ProcessingPool as Pool
import sys
import os

sys.path.append(os.getcwd())

cores = 10
start=time.time()
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/data/'
# load ninio time series
ninio34 = xr.open_dataset(RUTA + 'ninio34_index.nc')
# load  stratospheric index
PV = xr.open_dataset(RUTA + 'PV_index.nc')
# load antarctic points time series
antarct = xr.open_dataset(RUTA + 'hgt_points.nc')
antarct = antarct.drop('number')
year = np.arange(1981, 2017)
number = np.arange(0, 51)
realiz = pd.MultiIndex.from_product((year, number), names=('year', 'number'))

antarct = antarct.assign(realiz=realiz).unstack('realiz')
PV = PV.assign(dim_0=realiz).unstack('dim_0')
ninio34 = ninio34.assign(dim_0=realiz).unstack('dim_0')

"""
Following Byrne et al (2019) select one realization per year and compute regression between ninio and antarctic point conditioned and unconditioned on PV. 
do 10000
"""
ntimes = 10000
npoints = 2
nlevels = 2
nmonths = 7
i = np.repeat(np.arange(ntimes, dtype=int), npoints * nlevels * nmonths)
l = np.tile(np.repeat(np.arange(npoints, dtype=int), nlevels * nmonths), ntimes)
j = np.tile(np.repeat(np.arange(nlevels, dtype=int), nmonths), ntimes * npoints)
k = np.tile(np.arange(nmonths, dtype=int), ntimes * npoints * nlevels)
# create the multiprocessing pool
pool = Pool(cores) 
pool.close()
pool.join()
pool.clear()
def ComputeCorrelation(i, l, j, k, PV=PV, ninio34=ninio34, antarct=antarct):
	#select randomly one realiz per year
	PV_x = xr.DataArray([PV['PV_mon'].sel(year=p, number=randint(0, 50)) for p in range(1981,2017)])
	ninio34_x = xr.DataArray([ninio34['ninio34_mon'].sel(year=p, number=randint(0, 50)) for p in range(1981,2017)])
	antarct_y = xr.DataArray([antarct.z.sel(year=p, number=randint(0, 50)) for p in range(1981,2017)])
	#correlation antarct against PV
	R_antarct_PV = np.corrcoef(PV_x.values, antarct_y.values[:, l, j, k])[0, 1]
	#correlation antarct against ninio34
	R_antarct_ninio34 = np.corrcoef(ninio34_x.values, antarct_y.values[:, l, j, k])[0, 1]
	#regression antarct against PV
	A = np.vstack([PV_x.values, np.ones(36)])
	beta_antarct_PV, alpha_antarct_PV = np.linalg.lstsq(A.T, antarct_y.values[:, l, j, k], rcond=None)[0]
	#regression antarct against ninio34 condition
	antarct_residual = antarct_y.values[:, l, j, k] - (beta_antarct_PV * PV_x.values + alpha_antarct_PV)
	R_antarct_ninio34_cond = np.corrcoef(ninio34_x.values, antarct_residual)[0, 1]
	return R_antarct_PV, R_antarct_ninio34, R_antarct_ninio34_cond

print("Start")
# process the correlation by mapping
res = pool.map(ComputeCorrelation, i.tolist(), l.tolist(), j.tolist(), k.tolist())
print("Termino Paralelizacion")
print(time.time()-start)
# close down the pool and join
pool.close()
pool.join()
pool.clear()
#compute median and 5 and 92 percentile
res = np.stack(res, axis=1)
R_antarct_PV = np.reshape(res[0], [ntimes, npoints, nlevels, nmonths])
R_antarct_ninio34 = np.reshape(res[1], [ntimes, npoints, nlevels, nmonths])
R_antarct_ninio34_cond = np.reshape(res[2], [ntimes, npoints, nlevels, nmonths])
correlations = np.array([R_antarct_PV, R_antarct_ninio34, R_antarct_ninio34_cond])
median_corr = np.median(correlations, axis=1)
interval_corr = np.percentile(correlations, [5, 95], axis=1 )
np.savez('correlaciones.npz', corr=correlations, med=median_corr, interv=interval_corr)
