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
PV_index = xr.open_dataset(RUTA + 'PV_index.nc')
month = ['Aug', 'Sep', 'Oct', 'Nov']
seas = ['ASO', 'SON']
index_monthly_normal = np.logical_and(ninio34.ninio34_mon < ninio34.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34.ninio34_mon > ninio34.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear'))
# compute PV composites conditioned on non-EN anomalies
PV_index_EN = PV_index.sel(dim_0 = ~index_monthly_normal.values)
S_EN = S.sel(realiz=~index_monthly_normal.values)
uchi_EN = uchi.sel(realiz=~index_monthly_normal.values)
vchi_EN = vchi.sel(realiz=~index_monthly_normal.values)

index_monthly_upper = PV_index_EN.PV_mon >= PV_index_EN.PV_mon.quantile(0.75, dim='dim_0', interpolation='linear')
index_monthly_lower = PV_index_EN.PV_mon <= PV_index_EN.PV_mon.quantile(0.25, dim='dim_0', interpolation='linear')

for i in np.arange(0,4):
	var = np.mean(S_EN.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(S_EN.values[i, index_monthly_upper.values, :, :], axis=0)
	uchi_1 = np.mean(uchi_EN.u_chi.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(uchi_EN.u_chi.values[i, index_monthly_upper.values, :, :], axis=0)
	vchi_1 = np.mean(vchi_EN.v_chi.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(vchi_EN.v_chi.values[i, index_monthly_upper.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa SPV-WPV Years - ' + month[i] + ' - No ENSO'
	filename = './figures/200_RWS_chi_composites_diff_PV_' + month[i] +'_NoENSO.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

S = xr.open_dataset(RUTA + 'seasonal_RWS.nc', chunks={'latitude':10})
S = S.__xarray_dataarray_variable__.rename('S')
uchi = xr.open_dataset(RUTA + 'seasonal_uchi.nc', chunks={'latitude':10})
vchi = xr.open_dataset(RUTA + 'seasonal_vchi.nc', chunks={'latitude':10})
S_EN = S.sel(realiz=~index_monthly_normal.values)
uchi_EN = uchi.sel(realiz=~index_monthly_normal.values)
vchi_EN = vchi.sel(realiz=~index_monthly_normal.values)

for i in np.arange(0, 2):
	var = np.mean(S_EN.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(S_EN.values[i, index_monthly_upper.values, :, :], axis=0)
	uchi_1 = np.mean(uchi_EN.u_chi.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(uchi_EN.u_chi.values[i, index_monthly_upper.values, :, :], axis=0)
	vchi_1 = np.mean(vchi_EN.v_chi.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(vchi_EN.v_chi.values[i, index_monthly_upper.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa SPV-WPV Years - ' + seas[i] + ' - No ENSO'
	filename = './figures/200_RWS_chi_composites_diff_PV_' + seas[i] +'_NoENSO.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

