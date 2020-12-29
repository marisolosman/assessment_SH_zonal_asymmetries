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
index_EN = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
index_LN = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')

index_SPV_all = PV.SPV_index <= PV.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_WPV_all = PV.SPV_index >= PV.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')

index_SPV_EN = np.logical_and(index_SPV_all.values, index_EN.values)
index_SPV_LN = np.logical_and(index_SPV_all.values, index_LN.values)
index_WPV_EN = np.logical_and(index_WPV_all.values, index_EN.values)
index_WPV_LN = np.logical_and(index_WPV_all.values, index_LN.values)

S = xr.open_dataset(RUTA + 'monthly_RWS_50.nc4', chunks={'latitude':10})
S = S.transpose('month', 'latitude', 'longitude', 'realiz')
uchi = xr.open_dataset(RUTA + 'monthly_uchi_50.nc4', chunks={'latitude':10})
vchi = xr.open_dataset(RUTA + 'monthly_vchi_50.nc4', chunks={'latitude':10})
uchi = uchi.transpose('month', 'latitude', 'longitude', 'realiz')
vchi = vchi.transpose('month', 'latitude', 'longitude', 'realiz')

S = S.sel(latitude=slice(20, -90)).compute()
uchi = uchi.sel(latitude=slice(20, -90)).compute()
vchi = vchi.sel(latitude=slice(20, -90)).compute()

for i in np.arange(7):
	var_SPV_all = np.mean(S.RWS[i, :, :, index_SPV_all.values], axis=2)
	var_WPV_all = np.mean(S.RWS[i, :, :, index_WPV_all.values], axis=2)
	var_SPV_EN = np.mean(S.RWS[i, :, :, index_SPV_EN], axis=2)
	var_WPV_EN = np.mean(S.RWS[i, :, :, index_WPV_EN], axis=2)
	var_SPV_LN = np.mean(S.RWS[i, :, :, index_SPV_LN], axis=2)
	var_WPV_LN = np.mean(S.RWS[i, :, :, index_WPV_LN], axis=2)
	var_all = np.mean(S.RWS[i, :, :, :], axis=2)
	uchi_SPV_all = np.mean(uchi.u_chi[i, :, :, index_SPV_all.values], axis=2)
	vchi_SPV_all = np.mean(vchi.v_chi[i, :, :, index_SPV_all.values], axis=2)
	uchi_WPV_all = np.mean(uchi.u_chi[i, :, :, index_WPV_all.values], axis=2)
	vchi_WPV_all = np.mean(vchi.v_chi[i, :, :, index_WPV_all.values], axis=2)
	uchi_SPV_EN = np.mean(uchi.u_chi[i, :, :, index_SPV_EN], axis=2)
	vchi_SPV_EN = np.mean(vchi.v_chi[i, :, :, index_SPV_EN], axis=2)
	uchi_WPV_EN = np.mean(uchi.u_chi[i, :, :, index_WPV_EN], axis=2)
	vchi_WPV_EN = np.mean(vchi.v_chi[i, :, :, index_WPV_EN], axis=2)
	uchi_SPV_LN = np.mean(uchi.u_chi[i, :, :, index_SPV_LN], axis=2)
	vchi_SPV_LN = np.mean(vchi.v_chi[i, :, :, index_SPV_LN], axis=2)
	uchi_WPV_LN = np.mean(uchi.u_chi[i, :, :, index_WPV_LN], axis=2)
	vchi_WPV_LN = np.mean(vchi.v_chi[i, :, :, index_WPV_LN], axis=2)
	uchi_all = np.mean(uchi.u_chi[i, :, :, :], axis=2)
	vchi_all = np.mean(vchi.v_chi[i, :, :, :], axis=2)
	tit = 'Composites S4 RWS and Divergent wind 50hPa - ' + month[i] + ' SPoV'
	filename = FIG_PATH + '50_RWS_chi_composites_PV_' + month[i] +'_ENSO.png'
	plots.PlotCompositesRWSChiWPVENSO(var_SPV_all - var_all,
					       var_SPV_EN - var_all,
					       var_SPV_LN - var_all,
					       var_WPV_all - var_all,
					       var_WPV_EN - var_all,
					       var_WPV_LN - var_all,
					       uchi_SPV_all - uchi_all,
					       uchi_SPV_EN - uchi_all,
					       uchi_SPV_LN - uchi_all,
					       uchi_WPV_all - uchi_all,
					       uchi_WPV_EN - uchi_all,
					       uchi_WPV_LN - uchi_all,
					       vchi_SPV_all - vchi_all,
					       vchi_SPV_EN - vchi_all,
					       vchi_SPV_LN - vchi_all,
					       vchi_WPV_all - vchi_all,
					       vchi_WPV_EN - vchi_all,
					       vchi_WPV_LN - vchi_all,
					       S.latitude.values, S.longitude.values,
					       tit, filename)

S = xr.open_dataset(RUTA + 'seasonal_RWS_50.nc4', chunks={'latitude':10})
S = S.__xarray_dataarray_variable__.rename('RWS')
S = S.transpose('season', 'latitude', 'longitude', 'realiz')

uchi = xr.open_dataset(RUTA + 'seasonal_uchi.nc4', chunks={'latitude':10})
vchi = xr.open_dataset(RUTA + 'seasonal_vchi.nc4', chunks={'latitude':10})
S = S.sel(latitude=slice(20, -90)).compute()
uchi = uchi.sel(latitude=slice(20, -90)).compute()
uchi = uchi.transpose('season', 'latitude', 'longitude', 'realiz')

vchi = vchi.sel(latitude=slice(20, -90)).compute()
vchi = vchi.transpose('season', 'latitude', 'longitude', 'realiz')

for i in np.arange(5):
	var_SPV_all = np.mean(S[i, :, :, index_SPV_all.values], axis=2)
	var_WPV_all = np.mean(S[i, :, :, index_WPV_all.values], axis=2)
	var_SPV_EN = np.mean(S[i, :, :, index_SPV_EN], axis=2)
	var_WPV_EN = np.mean(S[i, :, :, index_WPV_EN], axis=2)
	var_SPV_LN = np.mean(S[i, :, :, index_SPV_LN], axis=2)
	var_WPV_LN = np.mean(S[i, :, :, index_WPV_LN], axis=2)
	var_all = np.mean(S[i, :, :, :], axis=2)
	uchi_SPV_all = np.mean(uchi.u_chi[i, :, :, index_SPV_all.values], axis=2)
	vchi_SPV_all = np.mean(vchi.v_chi[i, :, :, index_SPV_all.values], axis=2)
	uchi_WPV_all = np.mean(uchi.u_chi[i, :, :, index_WPV_all.values], axis=2)
	vchi_WPV_all = np.mean(vchi.v_chi[i, :, :, index_WPV_all.values], axis=2)
	uchi_SPV_EN = np.mean(uchi.u_chi[i, :, :, index_SPV_EN], axis=2)
	vchi_SPV_EN = np.mean(vchi.v_chi[i, :, :, index_SPV_EN], axis=2)
	uchi_WPV_EN = np.mean(uchi.u_chi[i, :, :, index_WPV_EN], axis=2)
	vchi_WPV_EN = np.mean(vchi.v_chi[i, :, :, index_WPV_EN], axis=2)
	uchi_SPV_LN = np.mean(uchi.u_chi[i, :, :, index_SPV_LN], axis=2)
	vchi_SPV_LN = np.mean(vchi.v_chi[i, :, :, index_SPV_LN], axis=2)
	uchi_WPV_LN = np.mean(uchi.u_chi[i, :, :, index_WPV_LN], axis=2)
	vchi_WPV_LN = np.mean(vchi.v_chi[i, :, :, index_WPV_LN], axis=2)
	uchi_all = np.mean(uchi.u_chi[i, :, :, :], axis=2)
	vchi_all = np.mean(vchi.v_chi[i, :, :, :], axis=2)
	tit = 'Composites S4 RWS and Divergent wind 50hPa - ' + seas[i] + ' SPoV'
	filename = FIG_PATH + '50_RWS_chi_composites_PV_' + seas[i] +'_ENSO.png'
	plots.PlotCompositesRWSChiWPVENSO(var_SPV_all - var_all,
					       var_SPV_EN - var_all,
					       var_SPV_LN - var_all,
					       var_WPV_all - var_all,
					       var_WPV_EN - var_all,
					       var_WPV_LN - var_all,
					       uchi_SPV_all - uchi_all,
					       uchi_SPV_EN - uchi_all,
					       uchi_SPV_LN - uchi_all,
					       uchi_WPV_all - uchi_all,
					       uchi_WPV_EN - uchi_all,
					       uchi_WPV_LN - uchi_all,
					       vchi_SPV_all - vchi_all,
					       vchi_SPV_EN - vchi_all,
					       vchi_SPV_LN - vchi_all,
					       vchi_WPV_all - vchi_all,
					       vchi_WPV_EN - vchi_all,
					       vchi_WPV_LN - vchi_all,
					       S.latitude.values, S.longitude.values,
					       tit, filename)


