#composites of z50 and plumb fluxes for EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import plots
import plumb_flux
#================================================
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '~/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/quartile_new/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
FILE_WINDS = 'fogt/monthly_winds200_aug_feb.nc4'
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
hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4, chunks={'latitude':10})
hgt = hgt - hgt.mean(dim='longitude')
hgt = hgt.sel(latitude=slice(10, -90)).compute()
winds = xr.open_dataset(PATH_DATA + FILE_WINDS, chunks={'latitude':10})
winds = winds.transpose('month', 'realiz', 'latitude', 'longitude')
winds = winds.sel(latitude=slice(10, -90)).compute()
winds_clm = winds.mean(dim='realiz')

#for i np.arange(0, 7):
#	var_normal = np.mean(hgt.z.values[i, :, :, :], axis=0)
#	var_ninio_WPV = np.mean(hgt.z.values[i, index_ninio_WPV, :, :], axis=0)
#	px_ninio_WPV, py_ninio_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_WPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
#	var_ninia_WPV = np.mean(hgt.z.values[i, index_ninia_WPV, :, :], axis=0)	
#	px_ninia_WPV, py_ninia_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninia_WPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
#	var_ninio_SPV = np.mean(hgt.z.values[i, index_ninio_SPV, :, :], axis=0)
#	px_ninio_SPV, py_ninio_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_SPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
#	var_ninia_SPV = np.mean(hgt.z.values[i, index_ninia_SPV, :, :], axis=0)	
#	px_ninia_SPV, py_ninia_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninia_SPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
#	var_ninio_all = np.mean(hgt.z.values[i, index_ninio_all, :, :], axis=0)
#	px_ninio_all, py_ninio_all, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_all-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
#	var_ninia_all = np.mean(hgt.z.values[i, index_ninia_all, :, :], axis=0)	
#	px_ninia_all, py_ninia_all, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninia_all-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
#	var ={'z1': var_ninio_all-var_normal, 'px1': px_ninio_all, 'py1': py_ninio_all,
#	      'z2': var_ninia_all-var_normal, 'px2': px_ninia_all, 'py2': py_ninia_all,
#	      'z3': var_ninio_WPV-var_normal, 'px3': px_ninio_WPV, 'py3': py_ninio_WPV,
#	      'z4': var_ninia_WPV-var_normal, 'px4': px_ninia_WPV, 'py4': py_ninia_WPV,
#	      'z5': var_ninio_SPV-var_normal, 'px5': px_ninio_SPV, 'py5': py_ninio_SPV,
#	      'z6': var_ninia_SPV-var_normal, 'px6': px_ninia_SPV, 'py6': py_ninia_SPV}
#	tit = 'Composites S4 Z* and Plumb Fluxes 50hPa Conditioned - SPoV - ' + month[i]
#	filename = FIG_PATH + 'z50_pf_composites_ENSO_' + month[i] +'_SPoV_q.png'
#	plots.PlotEnsoCompositesPoVZPF(var, hgt.latitude, hgt.longitude, tit, filename)

for i in np.arange(0, 7):
	var_normal = np.mean(hgt.z.values[i, :, :, :], axis=0)
	var_ninio_WPV = hgt.z.values[i, index_ninio_WPV, :, :]
	pf_ninio_WPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_WPV, var_normal, hgt.latitude.values, hgt.longitude.values, 20000)
	var_ninia_WPV = hgt.z.values[i, index_ninia_WPV, :, :]
	pf_ninia_WPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninia_WPV, var_normal, hgt.latitude.values, hgt.longitude.values, 20000)
	var_ninio_SPV = hgt.z.values[i, index_ninio_SPV, :, :]
	pf_ninio_SPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_SPV, var_normal, hgt.latitude.values, hgt.longitude.values, 20000)
	var_ninia_SPV = hgt.z.values[i, index_ninia_SPV, :, :]
	pf_ninia_SPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninia_SPV, var_normal, hgt.latitude.values, hgt.longitude.values, 20000)
	var_ninio_all = hgt.z.values[i, index_ninio_all, :, :]
	pf_ninio_all = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninio_all, var_normal, hgt.latitude.values, hgt.longitude.values, 20000)
	var_ninia_all = hgt.z.values[i, index_ninia_all, :, :]
	pf_ninia_all = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_ninia_all, var_normal, hgt.latitude.values, hgt.longitude.values, 20000)
	var ={'z1': np.mean(var_ninio_all, axis=0)-var_normal,
	      'px1': pf_ninio_all['pxl'], 'py1': pf_ninio_all['pyl'],
	      'z2': np.mean(var_ninia_all, axis=0) - var_normal,
	      'px2': pf_ninia_all['pxl'], 'py2': pf_ninia_all['pyl'],
	      'z3': np.mean(var_ninio_WPV, axis=0) - var_normal,
	      'px3': pf_ninio_WPV['pxl'], 'py3': pf_ninio_WPV['pyl'],
	      'z4': np.mean(var_ninia_WPV, axis=0) - var_normal,
	      'px4': pf_ninia_WPV['pxl'], 'py4': pf_ninia_WPV['pyl'],
	      'z5': np.mean(var_ninio_SPV, axis=0) - var_normal, 
	      'px5': pf_ninio_SPV['pxl'], 'py5': pf_ninio_SPV['pyl'],
	      'z6': np.mean(var_ninia_SPV, axis=0) - var_normal,
	      'px6': pf_ninia_SPV['pxl'], 'py6': pf_ninia_SPV['pyl']}
	tit = 'Composites S4 Z* and Plumb Fluxes 200hPa Conditioned - SPoV - ' + month[i]
	filename = FIG_PATH + 'z200_pf_composites_ENSO_' + month[i] +'_SPoV_lt.png'
	plots.PlotEnsoCompositesPoVZPF(var, hgt.latitude, hgt.longitude, tit, filename)
	var ={'z1': np.mean(var_ninio_all, axis=0)-var_normal,
	      'px1': pf_ninio_all['pxnl'], 'py1': pf_ninio_all['pynl'],
	      'z2': np.mean(var_ninia_all, axis=0) - var_normal,
	      'px2': pf_ninia_all['pxnl'], 'py2': pf_ninia_all['pynl'],
	      'z3': np.mean(var_ninio_WPV, axis=0) - var_normal,
	      'px3': pf_ninio_WPV['pxnl'], 'py3': pf_ninio_WPV['pynl'],
	      'z4': np.mean(var_ninia_WPV, axis=0) - var_normal,
	      'px4': pf_ninia_WPV['pxnl'], 'py4': pf_ninia_WPV['pynl'],
	      'z5': np.mean(var_ninio_SPV, axis=0) - var_normal, 
	      'px5': pf_ninio_SPV['pxnl'], 'py5': pf_ninio_SPV['pynl'],
	      'z6': np.mean(var_ninia_SPV, axis=0) - var_normal,
	      'px6': pf_ninia_SPV['pxnl'], 'py6': pf_ninia_SPV['pynl']}
	tit = 'Composites S4 Z* and Plumb Fluxes 200hPa Conditioned - SPoV - ' + month[i]
	filename = FIG_PATH + 'z200_pf_composites_ENSO_' + month[i] +'_SPoV_nlt.png'
	plots.PlotEnsoCompositesPoVZPF(var, hgt.latitude, hgt.longitude, tit, filename)

#for i in np.arange(0, 5):
#	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
#	winds_clm_s = winds_clm.isel(month=range(i, i+3)).mean(dim='month')
#	var_ninio_WPV = np.mean(hgt_s.z.values[index_ninio_WPV, :, :], axis=0)
#	var_normal = np.mean(hgt_s.z.values[:, :, :], axis=0)
#	px_ninio_WPV, py_ninio_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninio_WPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
#	var_ninia_WPV = np.mean(hgt_s.z.values[index_ninia_WPV, :, :], axis=0)
#	px_ninia_WPV, py_ninia_WPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninia_WPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
#	
#	var_ninio_SPV = np.mean(hgt_s.z.values[index_ninio_SPV, :, :], axis=0)
#	px_ninio_SPV, py_ninio_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninio_SPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
#	var_ninia_SPV = np.mean(hgt_s.z.values[index_ninia_SPV, :, :], axis=0)
#	px_ninia_SPV, py_ninia_SPV, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninia_SPV-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
#	var_ninio_all = np.mean(hgt_s.z.values[index_ninio_all, :, :], axis=0)
#	px_ninio_all, py_ninio_all, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninio_all-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
#	var_ninia_all = np.mean(hgt_s.z.values[index_ninia_all, :, :], axis=0)
#	px_ninia_all, py_ninia_all, lat = plumb_flux.ComputePlumbFluxes(winds_clm_s.u.values[:, :], winds_clm_s.v.values[:, :], var_ninia_all-var_normal, np.zeros_like(var_ninio_WPV), hgt.latitude.values, hgt.longitude.values, 5000)
#	tit = 'Composites S4 Z* and Plumb Fluxes 50hPa Conditioned - SPoV - ' + seas[i]
#	filename = FIG_PATH + 'z50_pf_composites_ENSO_' + seas[i] +'_SPoV_q.png'
#	var ={'z1': var_ninio_all-var_normal, 'px1': px_ninio_all, 'py1': py_ninio_all,
#	      'z2': var_ninia_all-var_normal, 'px2': px_ninia_all, 'py2': py_ninia_all,
#	      'z3': var_ninio_WPV-var_normal, 'px3': px_ninio_WPV, 'py3': py_ninio_WPV,
#	      'z4': var_ninia_WPV-var_normal, 'px4': px_ninia_WPV, 'py4': py_ninia_WPV,
#	      'z5': var_ninio_SPV-var_normal, 'px5': px_ninio_SPV, 'py5': py_ninio_SPV,
#	      'z6': var_ninia_SPV-var_normal, 'px6': px_ninia_SPV, 'py6': py_ninia_SPV}
#	plots.PlotEnsoCompositesPoVZPF(var, hgt.latitude, hgt.longitude, tit, filename)
#
