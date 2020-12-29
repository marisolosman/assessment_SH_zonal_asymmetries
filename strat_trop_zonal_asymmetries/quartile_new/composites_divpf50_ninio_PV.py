#composites of z200 and plumb fluxes for EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import plots
import random
import plumb_flux
from metpy import calc
#================================================
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '~/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/datos/figures/strat_trop_zonal_asymmetries/quartile_new/'
FILE_HGT_S4 = 'monthly_hgt50_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
FILE_WINDS = 'fogt/monthly_winds50_aug_feb.nc4'
hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for years with weak PV
index_SPV_upper = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')
#search for years with strong PV
index_SPV_normal = np.logical_and(PV_index.SPV_index > PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear'), PV_index.SPV_index < PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear'))

#enso during all years
index_ninio_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
index_ninia_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'))
#enso during weak PoV

index_ninio_WPV = np.logical_and(index_ninio_all.values, index_SPV_upper.values)
index_ninia_WPV = np.logical_and(index_ninia_all.values, index_SPV_upper.values)

#enso during strong PoV
index_ninio_SPV = np.logical_and(index_ninio_all.values, index_SPV_lower.values)
index_ninia_SPV = np.logical_and(index_ninia_all.values, index_SPV_lower.values)

index_normal = np.logical_and(index_SPV_normal, index_normal_all)

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

hgt = hgt - hgt.mean(dim='longitude')
winds = xr.open_dataset(PATH_DATA + FILE_WINDS, chunks={'latitude':10})
winds = winds.transpose('month', 'realiz', 'latitude', 'longitude')
winds_clm = winds.mean(dim='realiz')


for i in np.arange(0, 7):
	var_normal = np.mean(hgt.z.values[i, :, :, :], axis=0)
	print(var_normal.shape)
	var_ninio_WPV = np.mean(hgt.z.values[i, index_ninio_WPV, :, :], axis=0)
	print(var_ninio_WPV.shape)
	px_ninio_WPV, py_ninio_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_WPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
	dx, dy = calc.lat_lon_grid_deltas(hgt.longitude.values, lat)
	div_ninio_WPV = calc.divergence(px_ninio_WPV, py_ninio_WPV, dx, dy)

	var_ninia_WPV = np.mean(hgt.z.values[i, index_ninia_WPV, :, :], axis=0)	
	px_ninia_WPV, py_ninia_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninia_WPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
	div_ninia_WPV = calc.divergence(px_ninia_WPV, py_ninia_WPV, dx, dy)

	var_ninio_SPV = np.mean(hgt.z.values[i, index_ninio_SPV, :, :], axis=0)
	px_ninio_SPV, py_ninio_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_SPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
	div_ninio_SPV = calc.divergence(px_ninio_SPV, py_ninio_SPV, dx, dy)

	var_ninia_SPV = np.mean(hgt.z.values[i, index_ninia_SPV, :, :], axis=0)	
	px_ninia_SPV, py_ninia_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninia_SPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
	div_ninia_SPV = calc.divergence(px_ninia_SPV, py_ninia_SPV, dx, dy)

	var_ninio_all = np.mean(hgt.z.values[i, index_ninio_all, :, :], axis=0)
	px_ninio_all, py_ninio_all, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_all-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
	div_ninio_all = calc.divergence(px_ninio_all, py_ninio_all, dx, dy)

	var_ninia_all = np.mean(hgt.z.values[i, index_ninia_all, :, :], axis=0)	
	px_ninia_all, py_ninia_all, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninia_all-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
	div_ninia_all = calc.divergence(px_ninia_all, py_ninia_all, dx, dy)

	var ={'z1': div_ninio_all, 'px1': px_ninio_all, 'py1': py_ninio_all,
	      'z2': div_ninia_all, 'px2': px_ninia_all, 'py2': py_ninia_all,
	      'z3': div_ninio_WPV, 'px3': px_ninio_WPV, 'py3': py_ninio_WPV,
	      'z4': div_ninia_WPV, 'px4': px_ninia_WPV, 'py4': py_ninia_WPV,
	      'z5': div_ninio_SPV, 'px5': px_ninio_SPV, 'py5': py_ninio_SPV,
	      'z6': div_ninia_SPV, 'px6': px_ninia_SPV, 'py6': py_ninia_SPV}
	tit = 'Composites S4 Plumb Fluxes and divergencePF (*1e7) 50hPa Conditioned - SPoV - ' + month[i]
	filename = FIG_PATH + 'div_pf_50_composites_ENSO_' + month[i] +'_SPoV_q.png'
	plots.PlotEnsoCompositesPoVdivPF(var, lat, hgt.longitude, tit, filename)

#for i in np.arange(0, 5):
#	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
#	winds_clm_s = winds_clm.isel(month=range(i, i+3)).mean(dim='month')
#	var_ninio_WPV = np.mean(hgt_s.z.values[index_ninio_WPV, :, :], axis=0)
#	var_normal = np.mean(hgt_s.z.values[:, :, :], axis=0)
#	px_ninio_WPV, py_ninio_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninio_WPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
#	div_ninio_WPV = calc.divergence(px_ninio_WPV, py_ninio_WPV, dx, dy)
#
#	var_ninia_WPV = np.mean(hgt_s.z.values[index_ninia_WPV, :, :], axis=0)
#	px_ninia_WPV, py_ninia_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninia_WPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
#	div_ninia_WPV = calc.divergence(px_ninia_WPV, py_ninia_WPV, dx, dy)
#
#	var_ninio_SPV = np.mean(hgt_s.z.values[index_ninio_SPV, :, :], axis=0)
#	px_ninio_SPV, py_ninio_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninio_SPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
#	div_ninio_SPV = calc.divergence(px_ninio_SPV, py_ninio_SPV, dx, dy)
#
#	var_ninia_SPV = np.mean(hgt_s.z.values[index_ninia_SPV, :, :], axis=0)
#	px_ninia_SPV, py_ninia_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninia_SPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
#	div_ninia_SPV = calc.divergence(px_ninia_SPV, py_ninia_SPV, dx, dy)
#
#	var_ninio_all = np.mean(hgt_s.z.values[index_ninio_all, :, :], axis=0)
#	px_ninio_all, py_ninio_all, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninio_all-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
#	div_ninio_all = calc.divergence(px_ninio_all, py_ninio_all, dx, dy)
#
#	var_ninia_all = np.mean(hgt_s.z.values[index_ninia_all, :, :], axis=0)
#	px_ninia_all, py_ninia_all, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninia_all-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values)
#	div_ninia_all = calc.divergence(px_ninia_all, py_ninia_all, dx, dy)
#	
#	var ={'z1': div_ninio_all, 'px1': px_ninio_all, 'py1': py_ninio_all,
#	      'z2': div_ninia_all, 'px2': px_ninia_all, 'py2': py_ninia_all,
#	      'z3': div_ninio_WPV, 'px3': px_ninio_WPV, 'py3': py_ninio_WPV,
#	      'z4': div_ninia_WPV, 'px4': px_ninia_WPV, 'py4': py_ninia_WPV,
#	      'z5': div_ninio_SPV, 'px5': px_ninio_SPV, 'py5': py_ninio_SPV,
#	      'z6': div_ninia_SPV, 'px6': px_ninia_SPV, 'py6': py_ninia_SPV}
#	tit = 'Composites S4 Plumb Fluxes and divergencePF (*1e6) 50hPa Conditioned - SPoV - ' + month[i]
#	filename = FIG_PATH + 'div_pf_50_composites_ENSO_' + seas[i] +'_SPoV_q.png'
#	plots.PlotEnsoCompositesPoVdivPF(var, lat, hgt.longitude, tit, filename)
#
