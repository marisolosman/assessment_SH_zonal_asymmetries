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

#search for years with weak PV
index_monthly_upper = PV_index.PV_mon >= PV_index.PV_mon.quantile(0.75, dim='dim_0', interpolation='linear')
ninio34_WPV = ninio34.sel(dim_0=index_monthly_upper.values)
S_WPV = S.sel(realiz=index_monthly_upper.values)
uchi_WPV = uchi.sel(realiz=index_monthly_upper.values)
vchi_WPV = vchi.sel(realiz=index_monthly_upper.values)
index_monthly_normal = np.logical_and(ninio34_WPV.ninio34_mon < ninio34_WPV.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34_WPV.ninio34_mon > ninio34_WPV.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear'))

for i in np.arange(0,4):
	var = np.mean(S_WPV.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(S_WPV.values[i, index_monthly_normal.values, :, :], axis=0)
	uchi_1 = np.mean(uchi_WPV.u_chi.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(uchi_WPV.u_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	vchi_1 = np.mean(vchi_WPV.v_chi.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(vchi_WPV.v_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa EN+LN Years - ' + month[i] + ' - Weak PV'
	filename = './figures/200_RWS_chi_composites_sum_NINIO_' + month[i] +'_WPV.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

#search for years with strong PV
index_monthly_lower = PV_index.PV_mon <= PV_index.PV_mon.quantile(0.25, dim='dim_0', interpolation='linear')
ninio34_SPV = ninio34.sel(dim_0=index_monthly_upper.values)
S_SPV = S.sel(realiz=index_monthly_upper.values)
uchi_SPV = uchi.sel(realiz=index_monthly_upper.values)
vchi_SPV = vchi.sel(realiz=index_monthly_upper.values)
index_monthly_normal = np.logical_and(ninio34_SPV.ninio34_mon < ninio34_SPV.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34_SPV.ninio34_mon > ninio34_SPV.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear'))

for i in np.arange(0,4):
	var = np.mean(S_SPV.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(S_SPV.values[i, index_monthly_normal.values, :, :], axis=0)
	uchi_1 = np.mean(uchi_SPV.u_chi.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(uchi_SPV.u_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	vchi_1 = np.mean(vchi_SPV.v_chi.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(vchi_SPV.v_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa EN+LN Years - ' + month[i] + ' - Strong PV'
	filename = './figures/200_RWS_chi_composites_sum_NINIO_' + month[i] +'_SPV.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

S = xr.open_dataset(RUTA + 'seasonal_RWS.nc', chunks={'latitude':10})
S = S.__xarray_dataarray_variable__.rename('S')
uchi = xr.open_dataset(RUTA + 'seasonal_uchi.nc', chunks={'latitude':10})
vchi = xr.open_dataset(RUTA + 'seasonal_vchi.nc', chunks={'latitude':10})

#search for years with weak PV
index_monthly_upper = PV_index.PV_mon >= PV_index.PV_mon.quantile(0.75, dim='dim_0', interpolation='linear')
ninio34_WPV = ninio34.sel(dim_0=index_monthly_upper.values)
S_WPV = S.sel(realiz=index_monthly_upper.values)
uchi_WPV = uchi.sel(realiz=index_monthly_upper.values)
vchi_WPV = vchi.sel(realiz=index_monthly_upper.values)
index_monthly_normal = np.logical_and(ninio34_WPV.ninio34_mon < ninio34_WPV.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34_WPV.ninio34_mon > ninio34_WPV.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear'))

for i in np.arange(0, 2):
	var = np.mean(S_WPV.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(S_WPV.values[i, index_monthly_normal.values, :, :], axis=0)
	uchi_1 = np.mean(uchi_WPV.u_chi.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(uchi_WPV.u_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	vchi_1 = np.mean(vchi_WPV.v_chi.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(vchi_WPV.v_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa EN+LN Years - ' + seas[i] + ' - Weak PV'
	filename = './figures/200_RWS_chi_composites_sum_NINIO_' + seas[i] +'_WPV.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

#search for years with strong PV
index_monthly_lower = PV_index.PV_mon <= PV_index.PV_mon.quantile(0.25, dim='dim_0', interpolation='linear')
ninio34_SPV = ninio34.sel(dim_0=index_monthly_upper.values)
S_SPV = S.sel(realiz=index_monthly_upper.values)
uchi_SPV = uchi.sel(realiz=index_monthly_upper.values)
vchi_SPV = vchi.sel(realiz=index_monthly_upper.values)
index_monthly_normal = np.logical_and(ninio34_SPV.ninio34_mon < ninio34_SPV.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34_SPV.ninio34_mon > ninio34_SPV.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear'))

for i in np.arange(0, 2):
	var = np.mean(S_SPV.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(S_SPV.values[i, index_monthly_normal.values, :, :], axis=0)
	uchi_1 = np.mean(uchi_SPV.u_chi.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(uchi_SPV.u_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	vchi_1 = np.mean(vchi_SPV.v_chi.values[i, ~index_monthly_normal.values, :, :], axis=0) - np.mean(vchi_SPV.v_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa EN+LN Years - ' + seas[i] + ' - Strong PV'
	filename = './figures/200_RWS_chi_composites_sum_NINIO_' + seas[i] +'_SPV.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

