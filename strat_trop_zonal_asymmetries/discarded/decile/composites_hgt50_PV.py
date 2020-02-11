#composites HGT50 based on intensity of PV
import numpy as np
import xarray as xr
import os
import plots
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/data/'
ds = xr.open_dataset(RUTA + 'monthly_hgt50_aug_feb.nc')
PV = xr.open_dataset(RUTA + 'PV_index.nc')
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
index_monthly_upper = PV.PV_mon >= PV.PV_mon.quantile(0.9, dim='dim_0', interpolation='linear')
index_monthly_lower = PV.PV_mon <= PV.PV_mon.quantile(0.10, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(PV.PV_mon < PV.PV_mon.quantile(0.9, dim='dim_0', interpolation='linear'), PV.PV_mon > PV.PV_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0,7):
	var_neg = np.mean(ds.z[i, index_monthly_upper.values, :, :], axis=0) -\
			np.mean(ds.z[i, index_monthly_normal.values, :, :], axis=0)
	var_pos = np.mean(ds.z[i, index_monthly_lower.values, :, :], axis=0) -\
			np.mean(ds.z[i, index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 HGT 50hPa - ' + month[i]
	filename = './figures_decile/hgt_50_composites_PV_' + month[i] +'.png'
	plots.PlotPVComposites(var_pos, var_neg, ds.latitude, ds.longitude, tit, filename)
	var = np.mean(ds.z[i, index_monthly_lower.values, :, :], axis=0) -\
			np.mean(ds.z[i, index_monthly_upper.values, :, :], axis=0)
	tit = 'Composites differences SPV - WPV Years - ' + month[i]
	filename = './figures_decile/hgt_50_composites_diff_PV_' + month[i] +'.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)

for i in np.arange(0,5):
	var = ds.isel(month=range(i, i + 3)).mean(dim='month')
	var_neg = np.mean(var.z[index_monthly_upper.values, :, :], axis=0) -\
			np.mean(var.z[index_monthly_normal.values, :, :], axis=0)
	var_pos = np.mean(var.z[index_monthly_lower.values, :, :], axis=0) -\
			np.mean(var.z[index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 HGT 50hPa - ' + seas[i]
	filename = './figures_decile/hgt_50_composites_PV_' + seas[i] +'.png'
	plots.PlotPVComposites(var_pos, var_neg, ds.latitude, ds.longitude, tit, filename)
	var = np.mean(var.z[index_monthly_lower.values, :, :], axis=0) -\
			np.mean(var.z[index_monthly_upper.values, :, :], axis=0)
	tit = 'Composites differences SPV - WPV Years - ' + seas[i]
	filename = './figures_decile/hgt_50_composites_diff_PV_' + seas[i] +'.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)

