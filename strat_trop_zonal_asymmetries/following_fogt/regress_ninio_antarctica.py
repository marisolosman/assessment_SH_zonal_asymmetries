#!/usr/bin/python3
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

input1 = sys.argv[1]

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
npoints = 2
nlevels = 2
nmonths = 7
R_antarct_PV = np.empty( [npoints, nlevels, nmonths])
R_antarct_ninio34 = np.empty( [npoints, nlevels, nmonths])
R_antarct_ninio34_cond = np.empty([npoints, nlevels, nmonths])
print("Start")

for l in np.arange(npoints):
	for j in np.arange(nlevels):
		for k in np.arange(nmonths):
			PV_x = xr.DataArray([PV['PV_mon'].sel(year=p, number=randint(0, 50)) for p in range(1981,2017)])
			ninio34_x = xr.DataArray([ninio34['ninio34_mon'].sel(year=p, number=randint(0, 50)) for p in range(1981,2017)])
			antarct_y = xr.DataArray([antarct.z.sel(year=p, number=randint(0, 50)) for p in range(1981,2017)])
#correlation antarct against PV
			R_antarct_PV[l, j, k] = np.corrcoef(PV_x.values, antarct_y.values[:, l, j, k])[0, 1]
#correlation antarct against ninio34
			R_antarct_ninio34[l, j, k] = np.corrcoef(ninio34_x.values, antarct_y.values[:, l, j, k])[0, 1]
#regression antarct against PV
			A = np.vstack([PV_x.values, np.ones(36)])
			beta_antarct_PV, alpha_antarct_PV = np.linalg.lstsq(A.T, antarct_y.values[:, l, j, k], rcond=None)[0]
			#regression antarct against ninio34 condition
			antarct_residual = antarct_y.values[:, l, j, k] - (beta_antarct_PV * PV_x.values + alpha_antarct_PV)
			R_antarct_ninio34_cond[l, j, k] = np.corrcoef(ninio34_x.values, antarct_residual)[0, 1]

# process the correlation by mapping
print("End")
print(time.time()-start)
#compute median and 5 and 92 percentile
correlations = np.array([R_antarct_PV, R_antarct_ninio34, R_antarct_ninio34_cond])
median_corr = np.median(correlations, axis=1)
interval_corr = np.percentile(correlations, [5, 95], axis=1 )
np.savez('./correlaciones/correlations_' + input1 + '.npz', corr=correlations, med=median_corr, interv=interval_corr)

