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
PV = xr.open_dataset(RUTA + 'PV_index.nc')
month = ['Aug', 'Sep', 'Oct', 'Nov']
seas = ['ASO', 'SON']
index_monthly_upper = PV.PV_mon >= PV.PV_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = PV.PV_mon <= PV.PV_mon.quantile(0.10, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(PV.PV_mon < PV.PV_mon.quantile(0.90, dim='dim_0', interpolation='linear'), PV.PV_mon > PV.PV_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0,4):
	var_pos = np.mean(S.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(S.values[i, index_monthly_normal.values, :, :], axis=0)
	var_neg = np.mean(S.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(S.values[i, index_monthly_normal.values, :, :], axis=0)
	uchi_pos = np.mean(uchi.u_chi.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(uchi.u_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	uchi_neg = np.mean(uchi.u_chi.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(uchi.u_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	vchi_pos = np.mean(vchi.v_chi.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(vchi.v_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	vchi_neg = np.mean(vchi.v_chi.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(vchi.v_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 RWS and Divergent wind 200hPa - ' + month[i]
	filename = './figures_decile/200_RWS_chi_composites_PV_' + month[i] +'.png'
	plots.PlotCompositesDivPlumbPV(var_pos, var_neg, uchi_pos, uchi_neg, vchi_pos, vchi_neg,
				     S.latitude.values, S.longitude.values, tit, filename)
	var = np.mean(S.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(S.values[i, index_monthly_upper.values, :, :], axis=0)
	uchi_1 = np.mean(uchi.u_chi.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(uchi.u_chi.values[i, index_monthly_upper.values, :, :], axis=0)
	vchi_1 = np.mean(vchi.v_chi.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(vchi.v_chi.values[i, index_monthly_upper.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa differences SPV- WPV Years - ' + month[i]
	filename = './figures_decile/200_RWS_chi_composites_diff_PV_' + month[i] +'.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

S = xr.open_dataset(RUTA + 'seasonal_RWS.nc', chunks={'latitude':10})
S = S.__xarray_dataarray_variable__.rename('S')
uchi = xr.open_dataset(RUTA + 'seasonal_uchi.nc', chunks={'latitude':10})
vchi = xr.open_dataset(RUTA + 'seasonal_vchi.nc', chunks={'latitude':10})

for i in np.arange(0, 2):
	var_pos = np.mean(S.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(S.values[i, index_monthly_normal.values, :, :], axis=0)
	var_neg = np.mean(S.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(S.values[i, index_monthly_normal.values, :, :], axis=0)
	uchi_pos = np.mean(uchi.u_chi.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(uchi.u_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	uchi_neg = np.mean(uchi.u_chi.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(uchi.u_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	vchi_pos = np.mean(vchi.v_chi.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(vchi.v_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	vchi_neg = np.mean(vchi.v_chi.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(vchi.v_chi.values[i, index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 RWS and Divergent wind 200hPa - ' + seas[i]
	filename = './figures_decile/200_RWS_chi_composites_PV_' + seas[i] +'.png'
	plots.PlotCompositesDivPlumbPV(var_pos, var_neg, uchi_pos, uchi_neg, vchi_pos, vchi_neg,
				     S.latitude.values, S.longitude.values, tit, filename)
	var = np.mean(S.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(S.values[i, index_monthly_lower.values, :, :], axis=0)
	uchi_1 = np.mean(uchi.u_chi.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(uchi.u_chi.values[i, index_monthly_upper.values, :, :], axis=0)
	vchi_1 = np.mean(vchi.v_chi.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(vchi.v_chi.values[i, index_monthly_upper.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa differences SPV- WPV Years - ' + seas[i]
	filename = './figures_decile/200_RWS_chi_composites_diff_PV_' + seas[i] +'.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

