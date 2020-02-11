#composites of PV events conditioned on ENSO phase
import numpy as np
import xarray as xr
import os
import plots
import plumb_flux

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
FILE_WINDS = 'fogt/winds200_aug_feb.nc4'
hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
hgt = hgt - hgt.mean(dim='longitude')
winds = xr.open_dataset(PATH_DATA + FILE_WINDS)
winds = winds.transpose('month', 'realiz', 'latitude', 'longitude')
winds_clm = winds.mean(dim='realiz')
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for EN years 
index_EN = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for LN years
index_LN = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')

# compute SPoV composites conditioned on ENSO phase
PV_index_EN = PV_index.SPV_index.sel(dim_0 = index_EN.values)
hgt_EN = hgt.sel(realiz = index_EN.values)
PV_index_LN = PV_index.SPV_index.sel(dim_0 = index_LN.values)
hgt_LN = hgt.sel(realiz = index_LN.values)

#PV during all phases
index_SPV_all = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_WPV_all = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')
#PV  during EN
index_SPV_EN = PV_index_EN <= PV_index_EN.quantile(0.25, dim='dim_0', interpolation='linear')
index_WPV_EN = PV_index_EN >= PV_index_EN.quantile(0.75, dim='dim_0', interpolation='linear')
#PV  during LN
index_SPV_LN = PV_index_LN <= PV_index_LN.quantile(0.25, dim='dim_0', interpolation='linear')
index_WPV_LN = PV_index_LN >= PV_index_LN.quantile(0.75, dim='dim_0', interpolation='linear')

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

for i in np.arange(0, 7):
	var_WPV_EN = np.mean(hgt_EN.z.values[i, index_WPV_EN.values, :, :], axis=0)
	var_WPV_LN = np.mean(hgt_LN.z.values[i, index_WPV_LN.values, :, :], axis=0)	
	var_SPV_EN = np.mean(hgt_EN.z.values[i, index_SPV_EN.values, :, :], axis=0)
	var_SPV_LN = np.mean(hgt_LN.z.values[i, index_SPV_LN.values, :, :], axis=0)	
	var_WPV_all = np.mean(hgt.z.values[i, index_WPV_all.values, :, :], axis=0)
	var_SPV_all = np.mean(hgt.z.values[i, index_SPV_all.values, :, :], axis=0)
	px_EN, py_EN, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_SPV_EN - var_WPV_EN, np.zeros_like(var_WPV_EN), hgt.latitude.values, hgt.longitude.values)
	px_LN, py_LN, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_SPV_LN - var_WPV_LN, np.zeros_like(var_WPV_LN), hgt.latitude.values, hgt.longitude.values)
	px_all, py_all, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_SPV_all - var_WPV_all, np.zeros_like(var_WPV_all), hgt.latitude.values, hgt.longitude.values)
	var ={'z1': var_SPV_all-var_WPV_all, 'px1': px_all, 'py1': py_all,
	      'z2': var_SPV_EN-var_WPV_EN, 'px2': px_EN, 'py2': py_EN,
	      'z3': var_SPV_LN-var_WPV_LN, 'px3': px_LN, 'py3': py_LN}
	tit = 'Composites S4 Z* and Plumb Fluxes 200hPa S SPoV-Weak SPoV Conditioned - ENSO - ' + month[i]
	filename = FIG_PATH + 'z200_pf_composites_SPoV_' + month[i] +'_ENSO_q.png'

	plots.PlotPoVCompositesDiffENSOZPF(var, hgt.latitude, hgt.longitude, tit, filename)

for i in np.arange(0, 5):
	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
	hgt_s_EN = hgt_s.sel(realiz=index_EN.values)
	hgt_s_LN = hgt_s.sel(realiz=index_LN.values)
	var_WPV_EN = np.mean(hgt_s_EN.z.values[index_WPV_EN.values, :, :], axis=0)
	var_WPV_LN = np.mean(hgt_s_LN.z.values[index_WPV_LN.values, :, :], axis=0)	
	var_SPV_EN = np.mean(hgt_s_EN.z.values[index_SPV_EN.values, :, :], axis=0)
	var_SPV_LN = np.mean(hgt_s_LN.z.values[index_SPV_LN.values, :, :], axis=0)	
	var_WPV_all = np.mean(hgt_s.z.values[index_WPV_all.values, :, :], axis=0)
	var_SPV_all = np.mean(hgt_s.z.values[index_SPV_all.values, :, :], axis=0)
	tit = 'Composites S4 Z* 200hPa Strong SPoV - Weak SPoV Conditioned - ENSO - ' + seas[i]
	filename = FIG_PATH + 'z200_composites_SPoV_' + seas[i] +'_ENSO_q.png'
	plots.PlotPoVCompositesDiffENSO(var_SPV_all-var_WPV_all, var_SPV_EN - var_WPV_EN,
				    var_SPV_LN - var_WPV_LN, hgt.latitude, hgt.longitude,
				    tit, filename)


