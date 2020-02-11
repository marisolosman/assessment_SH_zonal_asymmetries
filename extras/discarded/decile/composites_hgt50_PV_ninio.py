#composites of PV events conditioned on ninio events
import numpy as np
import xarray as xr
import os
import plots

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '/storage/shared/glusterfs/acrcc/users/vg140344/data/'
#abro el archivo de geopotencial y junto la coordenada year y numbre
ds = xr.open_dataset(RUTA + 'monthly_hgt50_aug_feb.nc')
ninio34 = xr.open_dataset(RUTA + 'ninio34_index.nc')
PV_index = xr.open_dataset(RUTA + 'PV_index.nc')
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

#search for years with EN events
index_monthly_upper = ninio34.ninio34_mon >= ninio34.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear')

# compute PV composites conditioned on EN anomalies
PV_index_EN = PV_index.sel(dim_0 = index_monthly_upper.values)
ds_EN = ds.sel(realiz = index_monthly_upper.values)

index_monthly_upper = PV_index_EN.PV_mon >= PV_index_EN.PV_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = PV_index_EN.PV_mon <= PV_index_EN.PV_mon.quantile(0.10, dim='dim_0', interpolation='linear')
for i in np.arange(0,7):
	var = np.mean(ds_EN.z.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(ds_EN.z.values[i, index_monthly_upper.values, :, :], axis=0)
	tit = 'Composites HGT 50hPa differences SPW-WPV Years - El Niño - ' + month[i]
	filename = './figures_decile/hgt_50_composites_diff_PV_' + month[i] +'_EN.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)

for i in np.arange(0,5):
	var = ds_EN.isel(month=range(i, i + 3)).mean(dim='month')
	var = np.mean(var.z[index_monthly_lower.values, :, :], axis=0) - np.mean(var.z[index_monthly_upper.values, :, :], axis=0)
	tit = 'Composites HGT 50hPa differences SPV-WPV Years - El Niño - ' + seas[i]
	filename = './figures_decile/hgt_50_composites_diff_PV_' + seas[i] +'_EN.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)

#search for years with LN events
index_monthly_lower = ninio34.ninio34_mon <= ninio34.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear')

# compute PV composites conditioned on LN anomalies
PV_index_EN = PV_index.sel(dim_0 = index_monthly_lower.values)
ds_EN = ds.sel(realiz = index_monthly_lower.values)

index_monthly_upper = PV_index_EN.PV_mon >= PV_index_EN.PV_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = PV_index_EN.PV_mon <= PV_index_EN.PV_mon.quantile(0.10, dim='dim_0', interpolation='linear')

for i in np.arange(0,7):
	var = np.mean(ds_EN.z.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(ds_EN.z.values[i, index_monthly_upper.values, :, :], axis=0)
	tit = 'Composites differences SPV-WPV Years - La Niña - ' + month[i]
	filename = './figures_decile/hgt_50_composites_diff_PV_' + month[i] +'_LN.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)

for i in np.arange(0,5):
	var = ds_EN.isel(month=range(i, i + 3)).mean(dim='month')
	var = np.mean(var.z[index_monthly_lower.values, :, :], axis=0) - np.mean(var.z[index_monthly_upper.values, :, :], axis=0)
	tit = 'Composites differences SPV-WPV Years - La Niña - ' + seas[i]
	filename = './figures_decile/hgt_50_composites_diff_PV_' + seas[i] +'_LN.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)

