#sort el ninio events
import numpy as np
import xarray as xr
import pandas as pd
import plots
import os
from windspharm.xarray import VectorWind
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/quartile_new/'

PV = xr.open_dataset(RUTA + 'SPV_index.nc4')
ninio34 = xr.open_dataset(RUTA + 'ninio34_monthly.nc4')
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

#search for years with EN events
index_EN_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
index_LN_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')

index_SPV_all = PV.SPV_index <= PV.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_WPV_all = PV.SPV_index >= PV.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')

index_EN_SPV = np.logical_and(index_SPV_all.values, index_EN_all.values)
index_LN_SPV = np.logical_and(index_SPV_all.values, index_LN_all.values)
index_EN_WPV = np.logical_and(index_WPV_all.values, index_EN_all.values)
index_LN_WPV = np.logical_and(index_WPV_all.values, index_LN_all.values)

S = xr.open_dataset(RUTA + 'monthly_RWS_50.nc4', chunks={'latitude':10})
S = S.transpose('month', 'latitude', 'longitude', 'realiz')
uchi = xr.open_dataset(RUTA + 'monthly_uchi.nc4', chunks={'latitude':10})
vchi = xr.open_dataset(RUTA + 'monthly_vchi.nc4', chunks={'latitude':10})
S = S.sel(latitude=slice(20, -90)).compute()
vchi = vchi.transpose('month', 'latitude', 'longitude', 'realiz')
uchi = uchi.transpose('month', 'latitude', 'longitude', 'realiz')

uchi = uchi.sel(latitude=slice(20, -90)).compute()
vchi = vchi.sel(latitude=slice(20, -90)).compute()
for i in np.arange(7):
	var_EN_all = np.mean(S.RWS[i, :, :, index_EN_all.values], axis=2)
	var_LN_all = np.mean(S.RWS[i, :, :, index_LN_all.values], axis=2)
	var_EN_SPV = np.mean(S.RWS[i, :, :, index_EN_SPV], axis=2)
	var_EN_WPV = np.mean(S.RWS[i, :, :, index_EN_WPV], axis=2)
	var_LN_SPV = np.mean(S.RWS[i, :, :, index_LN_SPV], axis=2)
	var_LN_WPV = np.mean(S.RWS[i, :, :, index_LN_WPV], axis=2)
	var_all = np.mean(S.RWS[i, :, :, :], axis=2)
	uchi_EN_all = np.mean(uchi.u_chi[i, :, :, index_EN_all.values], axis=2)
	vchi_EN_all = np.mean(vchi.v_chi[i, :, :, index_EN_all.values], axis=2)
	uchi_LN_all = np.mean(uchi.u_chi[i, :, :, index_LN_all.values], axis=2)
	vchi_LN_all = np.mean(vchi.v_chi[i, :, :, index_LN_all.values], axis=2)
	uchi_EN_SPV = np.mean(uchi.u_chi[i, :, :, index_EN_SPV], axis=2)
	vchi_EN_SPV = np.mean(vchi.v_chi[i, :, :, index_EN_SPV], axis=2)
	uchi_EN_WPV = np.mean(uchi.u_chi[i, :, :, index_EN_WPV], axis=2)
	vchi_EN_WPV = np.mean(vchi.v_chi[i, :, :, index_EN_WPV], axis=2)
	uchi_LN_SPV = np.mean(uchi.u_chi[i, :, :, index_LN_SPV], axis=2)
	vchi_LN_SPV = np.mean(vchi.v_chi[i, :, :, index_LN_SPV], axis=2)
	uchi_LN_WPV = np.mean(uchi.u_chi[i, :, :, index_LN_WPV], axis=2)
	vchi_LN_WPV = np.mean(vchi.v_chi[i, :, :, index_LN_WPV], axis=2)
	uchi_all = np.mean(uchi.u_chi[i, :, :, :], axis=2)
	vchi_all = np.mean(vchi.v_chi[i, :, :, :], axis=2)
	tit = 'Composites S4 RWS and Divergent wind 50hPa - ' + month[i] + ' ENSO'
	filename = FIG_PATH + '50_RWS_chi_composites_ENSO_' + month[i] +'_SPoV.png'
	plots.PlotCompositesRWSChiWENSOPV(var_EN_all - var_all,
					       var_EN_WPV - var_all,
					       var_EN_SPV - var_all,
					       var_LN_all - var_all,
					       var_LN_WPV - var_all,
					       var_LN_SPV - var_all,
					       uchi_EN_all - uchi_all,
					       uchi_EN_WPV - uchi_all,
					       uchi_EN_SPV - uchi_all,
					       uchi_LN_all - uchi_all,
					       uchi_LN_WPV - uchi_all,
					       uchi_LN_SPV - uchi_all,
					       vchi_EN_all - vchi_all,
					       vchi_EN_WPV - vchi_all,
					       vchi_EN_SPV - vchi_all,
					       vchi_LN_all - vchi_all,
					       vchi_LN_WPV - vchi_all,
					       vchi_LN_SPV - vchi_all,
					       S.latitude.values, S.longitude.values,
					       tit, filename)

S = xr.open_dataset(RUTA + 'seasonal_RWS_50.nc4', chunks={'latitude':10})
S = S.__xarray_dataarray_variable__.rename('RWS')
S = S.transpose('season', 'latitude', 'longitude', 'realiz')

uchi = xr.open_dataset(RUTA + 'seasonal_uchi.nc4', chunks={'latitude':10})
vchi = xr.open_dataset(RUTA + 'seasonal_vchi.nc4', chunks={'latitude':10})
vchi = vchi.transpose('season', 'latitude', 'longitude', 'realiz')
uchi = uchi.transpose('season', 'latitude', 'longitude', 'realiz')

S = S.sel(latitude=slice(20, -90)).compute()
uchi = uchi.sel(latitude=slice(20, -90)).compute()
vchi = vchi.sel(latitude=slice(20, -90)).compute()

for i in np.arange(5):
	var_EN_all = np.mean(S[i, :, :, index_EN_all.values], axis=2)
	var_LN_all = np.mean(S[i, :, :, index_LN_all.values], axis=2)
	var_EN_SPV = np.mean(S[i, :, :, index_EN_SPV], axis=2)
	var_EN_WPV = np.mean(S[i, :, :, index_EN_WPV], axis=2)
	var_LN_SPV = np.mean(S[i, :, :, index_LN_SPV], axis=2)
	var_LN_WPV = np.mean(S[i, :, :, index_LN_WPV], axis=2)
	var_all = np.mean(S[i, :, :, :], axis=2)
	uchi_EN_all = np.mean(uchi.u_chi[i, :, :, index_EN_all.values], axis=2)
	vchi_EN_all = np.mean(vchi.v_chi[i, :, :, index_EN_all.values], axis=2)
	uchi_LN_all = np.mean(uchi.u_chi[i, :, :, index_LN_all.values], axis=2)
	vchi_LN_all = np.mean(vchi.v_chi[i, :, :, index_LN_all.values], axis=2)
	uchi_EN_SPV = np.mean(uchi.u_chi[i, :, :, index_EN_SPV], axis=2)
	vchi_EN_SPV = np.mean(vchi.v_chi[i, :, :, index_EN_SPV], axis=2)
	uchi_EN_WPV = np.mean(uchi.u_chi[i, :, :, index_EN_WPV], axis=2)
	vchi_EN_WPV = np.mean(vchi.v_chi[i, :, :, index_EN_WPV], axis=2)
	uchi_LN_SPV = np.mean(uchi.u_chi[i, :, :, index_LN_SPV], axis=2)
	vchi_LN_SPV = np.mean(vchi.v_chi[i, :, :, index_LN_SPV], axis=2)
	uchi_LN_WPV = np.mean(uchi.u_chi[i, :, :, index_LN_WPV], axis=2)
	vchi_LN_WPV = np.mean(vchi.v_chi[i, :, :, index_LN_WPV], axis=2)
	uchi_all = np.mean(uchi.u_chi[i, :, :, :], axis=2)
	vchi_all = np.mean(vchi.v_chi[i, :, :, :], axis=2)
	tit = 'Composites S4 RWS and Divergent wind 50hPa - ' + seas[i] + ' ENSO'
	filename = FIG_PATH + '50_RWS_chi_composites_ENSO_' + seas[i] +'_SPoV.png'
	plots.PlotCompositesRWSChiWENSOPV(var_EN_all - var_all,
					  var_EN_WPV - var_all,
					  var_EN_SPV - var_all,
					  var_LN_all - var_all,
					  var_LN_WPV - var_all,
					  var_LN_SPV - var_all,
					  uchi_EN_all - uchi_all,
					  uchi_EN_WPV - uchi_all,
					  uchi_EN_SPV - uchi_all,
					  uchi_LN_all - uchi_all,
					  uchi_LN_WPV - uchi_all,
					  uchi_LN_SPV - uchi_all,
					  vchi_EN_all - vchi_all,
					  vchi_EN_WPV - vchi_all,
					  vchi_EN_SPV - vchi_all,
					  vchi_LN_all - vchi_all,
					  vchi_LN_WPV - vchi_all,
					  vchi_LN_SPV - vchi_all,
					  S.latitude.values, S.longitude.values,
					  tit, filename)
