#sort el ninio events
import numpy as np
import xarray as xr
import plots
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA='~/datos/data/'
ds = xr.open_dataset(RUTA + 'monthly_hgt50_aug_feb.nc')
ds = ds - ds.mean(dim='longitude')
ninio34 = xr.open_dataset(RUTA + 'ninio34_index.nc')
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

index_monthly_upper = ninio34.ninio34_mon >= ninio34.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34.ninio34_mon <= ninio34.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(ninio34.ninio34_mon < ninio34.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34.ninio34_mon > ninio34.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0,7):
	var_neg = np.mean(ds.z.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0)
	var_pos = np.mean(ds.z.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 Z* 50hPa - ' + month[i]
	filename = './new_figures_decile/z50_composites_NINIO_' + month[i] +'.png'
	plots.PlotEnsoComposites(var_pos, var_neg, ds.latitude, ds.longitude, tit, filename)
	var = np.mean(ds.z.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(ds.z.values[i, index_monthly_lower.values, :, :], axis=0)
	tit = 'Composites Z* 50hPa diff EN-LN Years - ' + month[i]
	filename = './new_figures_decile/z50_composites_diff_NINIO_' + month[i] +'.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)

for i in np.arange(0, 5):
	var = ds.isel(month=range(i, i+3)).mean(dim='month')
	var_neg = np.mean(var.z[index_monthly_lower.values, :, :], axis=0) - np.mean(var.z[index_monthly_normal.values, :, :], axis=0)
	var_pos = np.mean(var.z[index_monthly_upper.values, :, :], axis=0) - np.mean(var.z[index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 Z* 50hPa - ' + seas[i]
	filename = './new_figures_decile/z50_composites_NINIO_' + seas[i] +'.png'
	plots.PlotEnsoComposites(var_pos, var_neg, ds.latitude, ds.longitude, tit, filename)
	var = np.mean(var.z[index_monthly_upper.values, :, :], axis=0) - np.mean(var.z[index_monthly_lower.values, :, :], axis=0)
	tit = 'Composites Z* 50hPa diff EN-LN Years - ' + seas[i]
	filename = './new_figures_decile/z50_composites_diff_NINIO_' + seas[i] +'.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)

