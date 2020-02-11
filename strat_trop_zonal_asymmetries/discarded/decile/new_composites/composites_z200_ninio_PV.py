#composites of EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import plots

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA='~/datos/data/'
#abro el archivo de geopotencial y junto la coordenada year y numbre
ds = xr.open_dataset(RUTA + 'monthly_hgt200_aug_feb.nc')
ds = ds - ds.mean(dim='longitude')
ninio34 = xr.open_dataset(RUTA + 'ninio34_index.nc')
PV_index = xr.open_dataset(RUTA + 'PV_index.nc')

#search for years with weak PV
index_monthly_upper = PV_index.PV_mon >= PV_index.PV_mon.quantile(0.90, dim='dim_0', interpolation='linear')

# compute EN-LA composites conditioned on PV anomalies
ninio34_WPV = ninio34.sel(dim_0 = index_monthly_upper.values)
ds_PV = ds.sel(realiz = index_monthly_upper.values)

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
index_monthly_upper = ninio34_WPV.ninio34_mon >= ninio34_WPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34_WPV.ninio34_mon <= ninio34_WPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(ninio34_WPV.ninio34_mon < ninio34_WPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_WPV.ninio34_mon > ninio34_WPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0, 7):
	var_pos = np.mean(ds_PV.z.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(ds_PV.z.values[i, index_monthly_normal.values, :, :], axis=0)
	var_neg = np.mean(ds_PV.z.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(ds_PV.z.values[i, index_monthly_normal.values, :, :], axis=0)

	tit = 'Composites S4 Z* 200hPa Conditioned - WPV - ' + month[i]
	filename = './new_figures_decile/z200_composites_NINIO_' + month[i] +'_WPV.png'
	plots.PlotEnsoComposites(var_pos, var_neg, ds.latitude, ds.longitude, tit, filename)

for i in np.arange(0, 5):
	var = ds_PV.isel(month=range(i, i+3)).mean(dim='month')
	var_pos = np.mean(var.z[index_monthly_upper.values, :, :], axis=0) - np.mean(var.z[index_monthly_normal.values, :, :], axis=0)
	var_neg = np.mean(var.z[index_monthly_lower.values, :, :], axis=0) - np.mean(var.z[index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 Z* 200hPa Conditioned - WPV - ' + seas[i]
	filename = './new_figures_decile/z200_composites_NINIO_' + seas[i] +'_WPV.png'
	plots.PlotEnsoComposites(var_pos, var_neg, ds.latitude, ds.longitude, tit, filename)
#search for years with strong PV

index_monthly_lower = PV_index.PV_mon <= PV_index.PV_mon.quantile(0.10, dim='dim_0', interpolation='linear')

# compute EN-LA composites conditioned on PV anomalies
ninio34_SPV = ninio34.sel(dim_0 = index_monthly_lower.values)
ds_PV = ds.sel(realiz = index_monthly_lower.values)

index_monthly_upper = ninio34_SPV.ninio34_mon >= ninio34_SPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34_SPV.ninio34_mon <= ninio34_SPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(ninio34_SPV.ninio34_mon < ninio34_SPV.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_SPV.ninio34_mon > ninio34_SPV.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0, 7):
	var_pos = np.mean(ds_PV.z.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(ds_PV.z.values[i, index_monthly_normal.values, :, :], axis=0)
	var_neg = np.mean(ds_PV.z.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(ds_PV.z.values[i, index_monthly_normal.values, :, :], axis=0)

	tit = 'Composites S4 Z* 200hPa Conditioned - SPV - ' + month[i]
	filename = './new_figures_decile/z200_composites_NINIO_' + month[i] +'_SPV.png'
	plots.PlotEnsoComposites(var_pos, var_neg, ds.latitude, ds.longitude, tit, filename)

for i in np.arange(0, 5):
	var = ds_PV.isel(month=range(i, i+3)).mean(dim='month')
	var_pos = np.mean(var.z[index_monthly_upper.values, :, :], axis=0) - np.mean(var.z[index_monthly_normal.values, :, :], axis=0)
	var_neg = np.mean(var.z[index_monthly_lower.values, :, :], axis=0) - np.mean(var.z[index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 Z* 200hPa Conditioned - SPV - ' + seas[i]
	filename = './new_figures_decile/z200_composites_NINIO_' + seas[i] +'_SPV.png'
	plots.PlotEnsoComposites(var_pos, var_neg, ds.latitude, ds.longitude, tit, filename)

