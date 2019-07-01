#sort el ninio events
import numpy as np
import xarray as xr
import pandas as pd
import plots
import os
from windspharm.xarray import VectorWind
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/data/'
S = xr.open_dataset(RUTA + 'monthly_RWS.nc', chunks={'latitude':10})
S = S.__xarray_dataarray_variable__.rename('S')
uchi = xr.open_dataset(RUTA + 'monthly_uchi.nc', chunks={'latitude':10})
vchi = xr.open_dataset(RUTA + 'monthly_vchi.nc', chunks={'latitude':10})
ninio34 = xr.open_dataset(RUTA + 'ninio34_index.nc')
month = ['Aug', 'Sep', 'Oct', 'Nov']
seas = ['ASO', 'SON']
index_monthly_normal = np.logical_and(ninio34.ninio34_mon < ninio34.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34.ninio34_mon > ninio34.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear'))

for i in np.arange(0,4):
	var = np.mean(S.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(S.values[i, index_monthly_normal.values, :, :], axis=0)
	uchi_1 = np.mean(uchi.u_chi.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(uchi.u_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	vchi_1 = np.mean(vchi.v_chi.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(vchi.v_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa EN+LN Years - ' + month[i]
	filename = './figures/200_RWS_chi_composites_sum_NINIO_' + month[i] +'.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

S = xr.open_dataset(RUTA + 'seasonal_RWS.nc', chunks={'latitude':10})
S = S.__xarray_dataarray_variable__.rename('S')
uchi = xr.open_dataset(RUTA + 'seasonal_uchi.nc', chunks={'latitude':10})
vchi = xr.open_dataset(RUTA + 'seasonal_vchi.nc', chunks={'latitude':10})

for i in np.arange(0, 2):
	var = np.mean(S.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(S.values[i, index_monthly_normal.values, :, :], axis=0)
	uchi_1 = np.mean(uchi.u_chi.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(uchi.u_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	vchi_1 = np.mean(vchi.v_chi.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(vchi.v_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa EN+LN Years - ' + seas[i]
	filename = './figures/200_RWS_chi_composites_sum_NINIO_' + seas[i] +'.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

