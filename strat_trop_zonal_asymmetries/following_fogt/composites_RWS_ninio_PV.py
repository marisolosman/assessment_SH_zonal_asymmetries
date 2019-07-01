#sort el ninio events
import numpy as np
import xarray as xr
import pandas as pd
import plots
import os
from windspharm.xarray import VectorWind
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/data/'
#abro el archivo de geopotencial y junto la coordenada year y numbre
ninio34 = xr.open_dataset(RUTA + 'ninio34_index.nc')
PV = xr.open_dataset(RUTA + 'PV_index.nc')

S = xr.open_dataset(RUTA + 'monthly_RWS.nc', chunks={'latitude':10})
S = S.__xarray_dataarray_variable__.rename('S')
uchi = xr.open_dataset(RUTA + 'monthly_uchi.nc', chunks={'latitude':10})
#uchi = uchi.__xarray_dataarray_variable__.rename('uchi')
vchi = xr.open_dataset(RUTA + 'monthly_vchi.nc', chunks={'latitude':10})
#vchi = vchi.__xarray_dataarray_variable__.rename('vchi')

month = ['Aug', 'Sep', 'Oct', 'Nov']
seas = ['ASO', 'SON']

index_PV_upper = PV.PV_mon >= PV.PV_mon.quantile(0.90, dim='dim_0', interpolation='linear')
ninio34_WPV = ninio34.sel(dim_0=index_PV_upper.values)
S_PV = S.sel(realiz=index_PV_upper.values)
uchi_PV = uchi.sel(realiz=index_PV_upper.values)
vchi_PV = vchi.sel(realiz=index_PV_upper.values)

index_monthly_upper = ninio34_WPV.ninio34_mon >= ninio34_WPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34_WPV.ninio34_mon <= ninio34_WPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(ninio34_WPV.ninio34_mon < ninio34_WPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_WPV.ninio34_mon > ninio34_WPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0,4):
	var = np.mean(S_PV.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(S_PV.values[i, index_monthly_lower.values, :, :], axis=0)
	uchi_1 = np.mean(uchi_PV.u_chi.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(uchi_PV.u_chi.values[i, index_monthly_lower.values, :, :], axis=0)
	vchi_1 = np.mean(vchi_PV.v_chi.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(vchi_PV.v_chi.values[i, index_monthly_lower.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa differences EN-LN Years - Weak PV - ' + month[i]
	filename = './figures_decile/200_RWS_chi_composites_diff_NINIO_' + month[i] +'_WPV.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

index_PV_lower = PV.PV_mon <= PV.PV_mon.quantile(0.10, dim='dim_0', interpolation='linear')
ninio34_SPV = ninio34.sel(dim_0=index_PV_lower.values)
S_PV = S.sel(realiz=index_PV_lower.values)
uchi_PV = uchi.sel(realiz=index_PV_lower.values)
vchi_PV = vchi.sel(realiz=index_PV_lower.values)

index_monthly_upper = ninio34_SPV.ninio34_mon >= ninio34_SPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34_SPV.ninio34_mon <= ninio34_SPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(ninio34_SPV.ninio34_mon < ninio34_SPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_SPV.ninio34_mon > ninio34_SPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0,4):
	var = np.mean(S_PV.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(S_PV.values[i, index_monthly_lower.values, :, :], axis=0)
	uchi_1 = np.mean(uchi_PV.u_chi.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(uchi_PV.u_chi.values[i, index_monthly_lower.values, :, :], axis=0)
	vchi_1 = np.mean(vchi_PV.v_chi.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(vchi_PV.v_chi.values[i, index_monthly_lower.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa differences EN-LN Years - Strong PV - ' + month[i]
	filename = './figures_decile/200_RWS_chi_composites_diff_NINIO_' + month[i] +'_SPV.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

S = xr.open_dataset(RUTA + 'seasonal_RWS.nc', chunks={'latitude':10})
S = S.__xarray_dataarray_variable__.rename('S')
uchi = xr.open_dataset(RUTA + 'seasonal_uchi.nc', chunks={'latitude':10})
#uchi = uchi.__xarray_dataarray_variable__.rename('uchi')
vchi = xr.open_dataset(RUTA + 'seasonal_vchi.nc', chunks={'latitude':10})
#vchi = vchi.__xarray_dataarray_variable__.rename('vchi')
S_PV = S.sel(realiz=index_PV_upper.values)
uchi_PV = uchi.sel(realiz=index_PV_upper.values)
vchi_PV = vchi.sel(realiz=index_PV_upper.values)

index_monthly_upper = ninio34_WPV.ninio34_mon >= ninio34_WPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34_WPV.ninio34_mon <= ninio34_WPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(ninio34_WPV.ninio34_mon < ninio34_WPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_WPV.ninio34_mon > ninio34_WPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0,2):
	var = np.mean(S_PV.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(S_PV.values[i, index_monthly_lower.values, :, :], axis=0)
	uchi_1 = np.mean(uchi_PV.u_chi.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(uchi_PV.u_chi.values[i, index_monthly_lower.values, :, :], axis=0)
	vchi_1 = np.mean(vchi_PV.v_chi.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(vchi_PV.v_chi.values[i, index_monthly_lower.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa differences EN-LN Years - Weak PV - ' + seas[i]
	filename = './figures_decile/200_RWS_chi_composites_diff_NINIO_' + seas[i] +'_WPV.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

S_PV = S.sel(realiz=index_PV_lower.values)
uchi_PV = uchi.sel(realiz=index_PV_lower.values)
vchi_PV = vchi.sel(realiz=index_PV_lower.values)

index_monthly_upper = ninio34_SPV.ninio34_mon >= ninio34_SPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34_SPV.ninio34_mon <= ninio34_SPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(ninio34_SPV.ninio34_mon < ninio34_SPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_SPV.ninio34_mon > ninio34_SPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0,2):
	var = np.mean(S_PV[i, index_monthly_upper.values, :, :], axis=0) - np.mean(S_PV[i, index_monthly_lower.values, :, :], axis=0)
	uchi_1 = np.mean(uchi_PV.u_chi.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(uchi_PV.u_chi.values[i, index_monthly_lower.values, :, :], axis=0)
	vchi_1 = np.mean(vchi_PV.v_chi.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(vchi_PV.v_chi.values[i, index_monthly_lower.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa differences EN-LN Years - Strong PV - ' + seas[i]
	filename = './figures_decile/200_RWS_chi_composites_diff_NINIO_' + seas[i] +'_SPV.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

