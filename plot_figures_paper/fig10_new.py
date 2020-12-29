#plot composites RWS using new layout suggested by Ted
import numpy as np
import xarray as xr
import pandas as pd
import plots_paper as plots
import os
from windspharm.xarray import VectorWind
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures_paper/'
PV = xr.open_dataset(RUTA + 'SPV_index.nc4')
ninio34 = xr.open_dataset(RUTA + 'ninio34_monthly.nc4')
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

#search for years with EN events
index_EN_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
index_LN_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')

index_SPV_all = PV.SPV_index <= PV.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_WPV_all = PV.SPV_index >= PV.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')

index_SPV_EN = np.logical_and(index_SPV_all.values, index_EN_all.values)
index_SPV_LN = np.logical_and(index_SPV_all.values, index_LN_all.values)
index_WPV_EN = np.logical_and(index_WPV_all.values, index_EN_all.values)
index_WPV_LN = np.logical_and(index_WPV_all.values, index_LN_all.values)

winds = xr.open_dataset(RUTA + 'winds200_aug_feb.nc4', chunks={'latitude':10})
#winds = winds.transpose('month', 'latitude', 'longitude', 'realiz')
uwind = winds.u.sel(latitude=(slice(0, -90))).compute()

S = xr.open_dataset(RUTA + 'monthly_RWS.nc4', chunks={'latitude':10})
uchi = xr.open_dataset(RUTA + 'monthly_uchi.nc4', chunks={'latitude':10})
vchi = xr.open_dataset(RUTA + 'monthly_vchi.nc4', chunks={'latitude':10})
S = S.transpose('month', 'latitude', 'longitude', 'realiz')
uchi = uchi.transpose('month', 'latitude', 'longitude', 'realiz')
vchi = vchi.transpose('month', 'latitude', 'longitude', 'realiz')

S = S.sel(latitude=slice(0, -90)).compute()
uchi = uchi.sel(latitude=slice(0, -90)).compute()
vchi = vchi.sel(latitude=slice(0, -90)).compute()

for i in [2]: #np.arange(7):
	var_EN_all = np.mean(S.RWS[i, :, :, index_EN_all.values], axis=2)
	var_LN_all = np.mean(S.RWS[i, :, :, index_LN_all.values], axis=2)
	var_SPV_all = np.mean(S.RWS[i, :, :, index_SPV_all.values], axis=2)
	var_WPV_all = np.mean(S.RWS[i, :, :, index_WPV_all.values], axis=2)
	var_LN_all = np.mean(S.RWS[i, :, :, index_LN_all.values], axis=2)
	var_SPV_EN = np.mean(S.RWS[i, :, :, index_SPV_EN], axis=2)
	var_WPV_EN = np.mean(S.RWS[i, :, :, index_WPV_EN], axis=2)
	var_SPV_LN = np.mean(S.RWS[i, :, :, index_SPV_LN], axis=2)
	var_WPV_LN = np.mean(S.RWS[i, :, :, index_WPV_LN], axis=2)
	var_all = np.mean(S.RWS[i, :, :, :], axis=2)

	uwind_EN_all = np.mean(uwind[i, :, :, index_EN_all.values], axis=2)
	uwind_LN_all = np.mean(uwind[i, :, :, index_LN_all.values], axis=2)
	uwind_SPV_all = np.mean(uwind[i, :, :, index_SPV_all.values], axis=2)
	uwind_WPV_all = np.mean(uwind[i, :, :, index_WPV_all.values], axis=2)
	uwind_LN_all = np.mean(uwind[i, :, :, index_LN_all.values], axis=2)
	uwind_SPV_EN = np.mean(uwind[i, :, :, index_SPV_EN], axis=2)
	uwind_WPV_EN = np.mean(uwind[i, :, :, index_WPV_EN], axis=2)
	uwind_SPV_LN = np.mean(uwind[i, :, :, index_SPV_LN], axis=2)
	uwind_WPV_LN = np.mean(uwind[i, :, :, index_WPV_LN], axis=2)
	uwind_all = np.mean(uwind[i, :, :, :], axis=2)

	uchi_EN_all = np.mean(uchi.u_chi[i, :, :, index_EN_all.values], axis=2)
	vchi_EN_all = np.mean(vchi.v_chi[i, :, :, index_EN_all.values], axis=2)
	uchi_SPV_all = np.mean(uchi.u_chi[i, :, :, index_SPV_all.values], axis=2)
	vchi_SPV_all = np.mean(vchi.v_chi[i, :, :, index_SPV_all.values], axis=2)
	uchi_WPV_all = np.mean(uchi.u_chi[i, :, :, index_WPV_all.values], axis=2)
	vchi_WPV_all = np.mean(vchi.v_chi[i, :, :, index_WPV_all.values], axis=2)
	uchi_LN_all = np.mean(uchi.u_chi[i, :, :, index_LN_all.values], axis=2)
	vchi_LN_all = np.mean(vchi.v_chi[i, :, :, index_LN_all.values], axis=2)
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
	tit = 'Composites S4 RWS (s-1) and Divergent wind (m/s) 200hPa - ' + month[i]
	filename = FIG_PATH + '200_RWS_chi_' + month[i] +'_new.eps'
	plots.PlotCompositesRWSChiWENSOPV(var_WPV_EN - var_SPV_EN,
					  var_WPV_EN - var_WPV_LN,
					  var_WPV_EN - var_all,
					  var_SPV_EN - var_all,
					  var_WPV_LN - var_SPV_LN,
					  var_SPV_EN - var_SPV_LN,
					  var_WPV_LN - var_all,
					  var_SPV_LN - var_all,
					  uchi_WPV_all - uchi_all,
					  uchi_EN_all - uchi_all,
					  uchi_WPV_EN - uchi_all,
					  uchi_SPV_EN - uchi_all,
					  uchi_SPV_all - uchi_all,
					  uchi_LN_all - uchi_all,
					  uchi_WPV_LN - uchi_all,
					  uchi_SPV_LN - uchi_all,
					  vchi_WPV_all - vchi_all,
					  vchi_EN_all - vchi_all,
					  vchi_WPV_EN - vchi_all,
					  vchi_SPV_EN - vchi_all,
					  vchi_SPV_all - vchi_all,
					  vchi_LN_all - vchi_all,
					  vchi_WPV_LN - vchi_all,
					  vchi_SPV_LN - vchi_all,

					  uwind_EN_all,
					  uwind_WPV_all,
					  uwind_WPV_EN,
					  uwind_SPV_EN,
					  uwind_LN_all,
					  uwind_SPV_all,
					  uwind_WPV_LN,
					  uwind_SPV_LN,
					  S.latitude.values, S.longitude.values,
					  tit, filename)


