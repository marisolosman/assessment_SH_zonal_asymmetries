#sort el ninio events
import numpy as np
import xarray as xr
import pandas as pd
import plots
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA='~/datos/data/'
#abro el archivo de geopotencial y junto la coordenada year y numbre
ds = xr.open_dataset(RUTA + 'monthly_hgt200_aug_feb.nc')
ds = ds - ds.mean(dim='longitude')
ninio34 = xr.open_dataset(RUTA + 'ninio34_index.nc')
month = ['Aug', 'Sep', 'Oct', 'Nov','Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

index_monthly_upper = ninio34.ninio34_mon >= ninio34.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34.ninio34_mon <= ninio34.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear')

index_monthly_normal = np.logical_and(ninio34.ninio34_mon < ninio34.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34.ninio34_mon > ninio34.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0,7):
	var = np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0)
	tit = 'S4 Z* 200hPa Neutral ENSO years - ' + month[i]
	filename = './new_figures_decile/z200_neutral_NINIO_' + month[i] +'.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)

for i in np.arange(0,5):
	var = ds.isel(month=range(i, i + 3)).mean(dim='month')
	var = np.mean(var.z[index_monthly_normal.values, :, :], axis=0)
	tit = 'S4 Z* 200hPa Neutral ENSO years - ' + seas[i]
	filename = './new_figures_decile/z200_neutral_NINIO_' + seas[i] +'.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)

