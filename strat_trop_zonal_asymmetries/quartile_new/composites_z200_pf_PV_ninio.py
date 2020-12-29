#composites of PV events conditioned on ENSO phase
import numpy as np
import xarray as xr
import os
import plots
import plumb_flux

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/quartile_new/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
FILE_WINDS = 'fogt/winds200_aug_feb.nc4'
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for EN years 
index_EN = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for LN years
index_LN = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')

index_neutral = np.logical_and(ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'), ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'))

#PV during all phases
index_SPV_all = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_WPV_all = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(PV_index.SPV_index < PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear'), PV_index.SPV_index > PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear'))
#PV  during EN
index_SPV_EN = np.logical_and(index_SPV_all.values, index_EN.values)
index_WPV_EN = np.logical_and(index_WPV_all.values, index_EN.values)

#PV  during LN
index_SPV_LN = np.logical_and(index_SPV_all.values, index_LN.values)
index_WPV_LN = np.logical_and(index_WPV_all.values, index_LN.values)

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
hgt = hgt - hgt.mean(dim='longitude')
winds = xr.open_dataset(PATH_DATA + FILE_WINDS)
winds = winds.transpose('month', 'realiz', 'latitude', 'longitude')
winds_clm = winds.mean(dim='realiz')

for i in np.arange(7):
	var_WPV_EN = np.mean(hgt.z.values[i, index_WPV_EN, :, :], axis=0)
	var_WPV_LN = np.mean(hgt.z.values[i, index_WPV_LN, :, :], axis=0)
	var_normal = np.mean(hgt.z.values[i, :, :, :], axis=0)
	var_SPV_EN = np.mean(hgt.z.values[i, index_SPV_EN, :, :], axis=0)
	var_SPV_LN = np.mean(hgt.z.values[i, index_SPV_LN, :, :], axis=0)	
	var_WPV_all = np.mean(hgt.z.values[i, index_WPV_all.values, :, :], axis=0)
	var_SPV_all = np.mean(hgt.z.values[i, index_SPV_all.values, :, :], axis=0)
	px_WPV, py_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_WPV_all - var_normal, np.zeros_like(var_WPV_EN), hgt.latitude.values, hgt.longitude.values)
	px_SPV, py_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_SPV_all - var_normal, np.zeros_like(var_WPV_EN), hgt.latitude.values, hgt.longitude.values)
	px_WPV_EN, py_WPV_EN, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_WPV_EN - var_normal, np.zeros_like(var_WPV_EN), hgt.latitude.values, hgt.longitude.values)
	px_SPV_EN, py_SPV_EN, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_SPV_EN - var_normal, np.zeros_like(var_WPV_EN), hgt.latitude.values, hgt.longitude.values)
	px_WPV_LN, py_WPV_LN, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_WPV_LN - var_normal, np.zeros_like(var_WPV_EN), hgt.latitude.values, hgt.longitude.values)
	px_SPV_LN, py_SPV_LN, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_SPV_LN - var_normal, np.zeros_like(var_WPV_EN), hgt.latitude.values, hgt.longitude.values)
	var ={'z1': var_SPV_all - var_normal, 'px1': px_SPV, 'py1': py_SPV,
	      'z2': var_SPV_EN - var_normal, 'px2': px_SPV_EN, 'py2': py_SPV_EN,
	      'z3': var_SPV_LN - var_normal, 'px3': px_SPV_LN, 'py3': py_SPV_LN,
	      'z4': var_WPV_all - var_normal, 'px4': px_WPV, 'py4': py_WPV,
	      'z5': var_WPV_EN - var_normal, 'px5': px_WPV_EN, 'py5': py_WPV_EN,
	      'z6': var_WPV_LN - var_normal, 'px6': px_WPV_LN, 'py6': py_WPV_LN}
	tit = 'Composites S4 Z* and Plumb Fluxes 200hPa SPoV Conditioned - ENSO - ' + month[i]
	filename = FIG_PATH + 'z200_pf_composites_SPoV_' + month[i] +'_ENSO_q.png'

	plots.PlotPoVCompositesENSOZPF(var, hgt.latitude, hgt.longitude, tit, filename)

for i in np.arange(0, 5):
	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
	winds_clm_s = winds_clm.isel(month=range(i, i + 3)).mean(dim='month')
	var_WPV_EN = np.mean(hgt_s.z.values[index_WPV_EN, :, :], axis=0)
	var_WPV_LN = np.mean(hgt_s.z.values[index_WPV_LN, :, :], axis=0)	
	var_SPV_EN = np.mean(hgt_s.z.values[index_SPV_EN, :, :], axis=0)
	var_SPV_LN = np.mean(hgt_s.z.values[index_SPV_LN, :, :], axis=0)	
	var_WPV_all = np.mean(hgt_s.z.values[index_WPV_all.values, :, :], axis=0)
	var_SPV_all = np.mean(hgt_s.z.values[index_SPV_all.values, :, :], axis=0)
	var_normal = np.mean(hgt_s.z.values[:, :, :], axis=0)
	px_WPV, py_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values, winds_clm_s.v.values, var_WPV_all - var_normal, np.zeros_like(var_WPV_EN), hgt.latitude.values, hgt.longitude.values)
	px_SPV, py_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values, winds_clm_s.v.values, var_SPV_all - var_normal, np.zeros_like(var_WPV_EN), hgt.latitude.values, hgt.longitude.values)
	px_WPV_EN, py_WPV_EN, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values, winds_clm_s.v.values, var_WPV_EN - var_normal, np.zeros_like(var_WPV_EN), hgt.latitude.values, hgt.longitude.values)
	px_SPV_EN, py_SPV_EN, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values, winds_clm_s.v.values, var_SPV_EN - var_normal, np.zeros_like(var_WPV_EN), hgt.latitude.values, hgt.longitude.values)
	px_WPV_LN, py_WPV_LN, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values, winds_clm_s.v.values, var_WPV_LN - var_normal, np.zeros_like(var_WPV_EN), hgt.latitude.values, hgt.longitude.values)
	px_SPV_LN, py_SPV_LN, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values, winds_clm_s.v.values, var_SPV_LN - var_normal, np.zeros_like(var_WPV_EN), hgt.latitude.values, hgt.longitude.values)
	var ={'z1': var_SPV_all - var_normal, 'px1': px_SPV, 'py1': py_SPV,
	      'z2': var_SPV_EN - var_normal, 'px2': px_SPV_EN, 'py2': py_SPV_EN,
	      'z3': var_SPV_LN - var_normal, 'px3': px_SPV_LN, 'py3': py_SPV_LN,
	      'z4': var_WPV_all - var_normal, 'px4': px_WPV, 'py4': py_WPV,
	      'z5': var_WPV_EN - var_normal, 'px5': px_WPV_EN, 'py5': py_WPV_EN,
	      'z6': var_WPV_LN - var_normal, 'px6': px_WPV_LN, 'py6': py_WPV_LN}
	filename = FIG_PATH + 'z200_pf_composites_SPoV_' + seas[i] +'_ENSO_q.png'
	tit = 'Composites S4 Z* and Plum Fluxes 200hPa SPoV Conditioned - ENSO - ' + seas[i]
	plots.PlotPoVCompositesENSOZPF(var, hgt.latitude, hgt.longitude, tit, filename)

