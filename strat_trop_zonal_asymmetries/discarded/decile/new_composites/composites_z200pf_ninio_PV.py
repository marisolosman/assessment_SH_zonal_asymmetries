#sort el ninio events
import numpy as np
import xarray as xr
import plots
import os
import plumb_flux
RUTA='~/datos/data/'
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
winds = xr.open_dataset(RUTA + 'monthly_winds200_aug_feb.nc')
winds_clm = winds.mean(dim='realiz')
#abro el archivo de geopotencial y junto la coordenada year y numbre
ds = xr.open_dataset(RUTA + 'monthly_hgt200_aug_feb.nc')
ds = ds - ds.mean(dim='longitude')
ninio34 = xr.open_dataset(RUTA + 'ninio34_index.nc')
PV = xr.open_dataset(RUTA + 'PV_index.nc')
index_PV_upper = PV.PV_mon >= PV.PV_mon.quantile(0.90, dim='dim_0', interpolation='linear')
ninio34_WPV = ninio34.sel(dim_0=index_PV_upper.values)

ds_PV = ds.sel(realiz=index_PV_upper.values)

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

index_monthly_upper = ninio34_WPV.ninio34_mon >= ninio34_WPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34_WPV.ninio34_mon <= ninio34_WPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(ninio34_WPV.ninio34_mon < ninio34_WPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_WPV.ninio34_mon > ninio34_WPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0,7):
	var_pos = np.mean(ds_PV.z.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(ds_PV.z.values[i, index_monthly_normal.values, :, :], axis=0)
	px_pos, py_pos, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_pos, np.zeros_like(var_pos), ds.latitude.values, ds.longitude.values)
	var_neg = np.mean(ds_PV.z.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(ds_PV.z.values[i, index_monthly_normal.values, :, :], axis=0)
	px_neg, py_neg, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_neg, np.zeros_like(var_neg), ds.latitude.values, ds.longitude.values)
	tit = 'Composites S4 Z* Plumb 200hPa - Weak SPV - ' + month[i]
	filename = './new_figures_decile/z200_plumb_composites_NINIO_' + month[i] +'_WPV.png'
	plots.PlotCompositesPlumb(var_pos, var_neg, px_pos, px_neg, py_pos, py_neg, ds.latitude, ds.longitude, tit, filename)
#

for i in np.arange(0,5):
	var = ds_PV.isel(month=range(i, i +3 )).mean(dim='month')
	var_wnd = winds_clm.isel(month=range(i, i +3 )).mean(dim='month')
	var_pos = np.mean(var.z.values[index_monthly_upper.values, :, :], axis=0) - np.mean(var.z.values[index_monthly_normal.values, :, :], axis=0)
	px_pos, py_pos, lat = plumb_flux.ComputePlumbFluxes(var_wnd.u.values[:, :], var_wnd.v.values[:, :], var_pos, np.zeros_like(var_pos), ds.latitude.values, ds.longitude.values)
	var_neg = np.mean(var.z.values[index_monthly_lower.values, :, :], axis=0) - np.mean(var.z.values[index_monthly_normal.values, :, :], axis=0)

	px_neg, py_neg, lat = plumb_flux.ComputePlumbFluxes(var_wnd.u.values[:, :], var_wnd.v.values[:, :], var_neg, np.zeros_like(var_neg), ds.latitude.values, ds.longitude.values)
	tit = 'Composites S4 Z* Plumb 200hPa - Weak SPV - ' + seas[i]
	filename = './new_figures_decile/z200_plumb_composites_NINIO_' + seas[i] +'_WPV.png'
	plots.PlotCompositesPlumb(var_pos, var_neg, px_pos, px_neg, py_pos, py_neg, ds.latitude, ds.longitude, tit, filename)


index_PV_lower = PV.PV_mon <= PV.PV_mon.quantile(0.10, dim='dim_0', interpolation='linear')
ninio34_SPV = ninio34.sel(dim_0=index_PV_lower.values)
ds_PV = ds.sel(realiz=index_PV_lower.values)

index_monthly_upper = ninio34_SPV.ninio34_mon >= ninio34_SPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34_SPV.ninio34_mon <= ninio34_SPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(ninio34_SPV.ninio34_mon < ninio34_SPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_SPV.ninio34_mon > ninio34_SPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0,7):
	var_pos = np.mean(ds_PV.z.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(ds_PV.z.values[i, index_monthly_normal.values, :, :], axis=0)
	px_pos, py_pos, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_pos, np.zeros_like(var_pos), ds.latitude.values, ds.longitude.values)
	var_neg = np.mean(ds_PV.z.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(ds_PV.z.values[i, index_monthly_normal.values, :, :], axis=0)
	px_neg, py_neg, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var_neg, np.zeros_like(var_neg), ds.latitude.values, ds.longitude.values)
	tit = 'Composites S4 Z* Plumb 200hPa - Strong SPV - ' + month[i]
	filename = './new_figures_decile/z200_plumb_composites_NINIO_' + month[i] +'_SPV.png'
	plots.PlotCompositesPlumb(var_pos, var_neg, px_pos, px_neg, py_pos, py_neg, ds.latitude, ds.longitude, tit, filename)
#

for i in np.arange(0,5):
	var = ds_PV.isel(month=range(i, i +3 )).mean(dim='month')
	var_wnd = winds_clm.isel(month=range(i, i +3 )).mean(dim='month')
	var_pos = np.mean(var.z.values[index_monthly_upper.values, :, :], axis=0) - np.mean(var.z.values[index_monthly_normal.values, :, :], axis=0)
	px_pos, py_pos, lat = plumb_flux.ComputePlumbFluxes(var_wnd.u.values[:, :], var_wnd.v.values[:, :], var_pos, np.zeros_like(var_pos), ds.latitude.values, ds.longitude.values)
	var_neg = np.mean(var.z.values[index_monthly_lower.values, :, :], axis=0) - np.mean(var.z.values[index_monthly_normal.values, :, :], axis=0)
	px_neg, py_neg, lat = plumb_flux.ComputePlumbFluxes(var_wnd.u.values[:, :], var_wnd.v.values[:, :], var_neg, np.zeros_like(var_neg), ds.latitude.values, ds.longitude.values)
	tit = 'Composites S4 Z* Plumb 200hPa - Strong SPV - ' + seas[i]
	filename = './new_figures_decile/z200_plumb_composites_NINIO_' + seas[i] +'_SPV.png'
	plots.PlotCompositesPlumb(var_pos, var_neg, px_pos, px_neg, py_pos, py_neg, ds.latitude, ds.longitude, tit, filename)


