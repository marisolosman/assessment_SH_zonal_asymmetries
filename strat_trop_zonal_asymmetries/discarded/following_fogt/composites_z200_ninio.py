#sort el ninio events
import numpy as np
import xarray as xr
import numpy.ma as ma
from scipy import stats
import plots
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA='~/datos/data/'
ds = xr.open_dataset(RUTA + 'monthly_hgt200_aug_feb.nc')
ds = ds - ds.mean(dim='longitude')
ninio34 = xr.open_dataset(RUTA + 'fogt/' + 'ninio34_index.nc')
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

index_monthly_upper = ninio34.ninio34_mon >= ninio34.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34.ninio34_mon <= ninio34.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(ninio34.ninio34_mon < ninio34.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34.ninio34_mon > ninio34.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0,7):
	var_neg = np.mean(ds.z.values[i, index_monthly_lower.values, :, :] - np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0), axis=0)
	ttest_neg, pval_neg = stats.ttest_ind(ds.z.values[i, index_monthly_lower.values, :, :],
					      ds.z.values[i, index_monthly_normal.values, :, :])
	var_pos = np.mean(ds.z.values[i, index_monthly_upper.values, :, :] - np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0), axis=0)
	ttest_pos, pval_pos = stats.ttest_ind(ds.z.values[i, index_monthly_upper.values, :, :],
					      ds.z.values[i, index_monthly_normal.values, :, :])
	tit = 'Composites S4 Z* 200hPa - ' + month[i]
	filename = './figures/z200_composites_NINIO_' + month[i] +'.png'
	var_pos = ma.masked_array(var_pos, mask=np.abs(pval_pos)>0.025)
	var_neg = ma.masked_array(var_neg, mask=np.abs(pval_neg)>0.025)
	plots.PlotEnsoComposites(var_pos, var_neg, ds.latitude, ds.longitude, tit, filename)
	var = np.mean(ds.z.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(ds.z.values[i, index_monthly_lower.values, :, :], axis=0)
	tit = 'Composites S4 Z* 200hPa differences EN-LN Years - ' + month[i]
	filename = './figures/z200_composites_diff_NINIO_' + month[i] +'.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)

for i in np.arange(0, 5):
	var = ds.isel(month=range(i, i+3)).mean(dim='month')
	var_neg = np.mean(var.z[index_monthly_lower.values, :, :] - np.mean(var.z[index_monthly_normal.values, :, :], axis=0), axis=0)
	var_pos = np.mean(var.z[index_monthly_upper.values, :, :] - np.mean(var.z[index_monthly_normal.values, :, :], axis=0), axis=0)
	ttest_pos, pval_pos = stats.ttest_ind(var.z[index_monthly_upper.values, :, :],
					      var.z.values[index_monthly_normal.values, :, :])
	ttest_neg, pval_neg = stats.ttest_ind(var.z[index_monthly_lower.values, :, :],
					      var.z.values[index_monthly_normal.values, :, :])
	var_pos = ma.masked_array(var_pos, mask=np.abs(pval_pos)>0.025)
	var_neg = ma.masked_array(var_neg, mask=np.abs(pval_neg)>0.025)
	tit = 'Composites S4 Z* 200hPa - ' + seas[i]
	filename = './figures/z200_composites_NINIO_' + seas[i] +'.png'
	plots.PlotEnsoComposites(var_pos, var_neg, ds.latitude, ds.longitude, tit, filename)
	var = np.mean(var.z[index_monthly_upper.values, :, :], axis=0) - np.mean(var.z[index_monthly_lower.values, :, :], axis=0)
	tit = 'Composites S4 Z* 200hPa differences EN-LN Years - ' + seas[i]
	filename = './figures/z200_composites_diff_NINIO_' + seas[i] +'.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)

