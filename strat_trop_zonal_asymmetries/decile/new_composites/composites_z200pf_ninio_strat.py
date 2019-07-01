import numpy as np
import xarray as xr
import pandas as pd
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
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

index_monthly_normal = np.logical_and(ninio34.ninio34_mon < ninio34.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34.ninio34_mon > ninio34.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear'))

for i in np.arange(0,7):
	var = np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0)
	px_c, py_c, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var, np.zeros_like(var), ds.latitude.values, ds.longitude.values)
	tit = 'S4 Z* Plumb Fluxes Neutral ENSO Years - ' + month[i]
	filename = './new_figures_decile/z200_plumb_neutral_NINIO_' + month[i] +'.png'
	plots.PlotCompPlumbDiff(var, px_c, py_c, ds.latitude, ds.longitude, tit, filename)

for i in np.arange(0,5):
	var = ds.isel(month=range(i, i + 3)).mean(dim='month')
	var_wnd = winds_clm.isel(month=range(i, i + 3)).mean(dim='month')
	var = np.mean(var.z.values[index_monthly_normal.values, :, :], axis=0)
	px_c, py_c, lat = plumb_flux.ComputePlumbFluxes(var_wnd.u.values[:, :], var_wnd.v.values[:, :], var, np.zeros_like(var), ds.latitude.values, ds.longitude.values)
	tit = 'S4 Z* Plumb Fluxes Neutral ENSO Years - ' + seas[i]
	filename = './new_figures_decile/z200_plumb_neutral_NINIO_' + seas[i] +'.png'
	plots.PlotCompPlumbDiff(var, px_c, py_c, ds.latitude, ds.longitude, tit, filename)

