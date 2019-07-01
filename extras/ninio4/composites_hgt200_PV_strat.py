#composites of PV events conditioned on ninio events
import numpy as np
import xarray as xr
import pandas as pd
import plots
RUTA='~/data/'
#abro el archivo de geopotencial y junto la coordenada year y numbre
ds = xr.open_dataset(RUTA + 'monthly_hgt200_aug_feb.nc')
ninio34 = xr.open_dataset(RUTA + 'ninio34_index.nc')
PV_index = xr.open_dataset(RUTA + 'PV_index.nc')
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

index_monthly_normal = np.logical_and(ninio34.ninio34_mon < ninio34.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34.ninio34_mon > ninio34.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear'))
# compute PV composites conditioned on non-EN anomalies
PV_index_EN = PV_index.sel(dim_0 = ~index_monthly_normal.values)
ds_EN = ds.sel(realiz = ~index_monthly_normal.values)

index_monthly_upper = PV_index_EN.PV_mon >= PV_index_EN.PV_mon.quantile(0.75, dim='dim_0', interpolation='linear')
index_monthly_lower = PV_index_EN.PV_mon <= PV_index_EN.PV_mon.quantile(0.25, dim='dim_0', interpolation='linear')
for i in np.arange(0,7):
	var = np.mean(ds_EN.z.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(ds_EN.z.values[i, index_monthly_upper.values, :, :], axis=0)
	tit = 'Composites differences SPV-WPV Years - ' + month[i] + '- No ENSO'
	filename = './figures/hgt_200_composites_diff_PV_' + month[i] +'_NoENSO.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)

for i in np.arange(0,5):
	var = ds_EN.isel(month=range(i, i + 3)).mean(dim='month')
	var = np.mean(var.z[index_monthly_lower.values, :, :], axis=0) - np.mean(var.z[index_monthly_upper.values, :, :], axis=0)
	tit = 'Composites differences SPV-WPV Years - ' + seas[i] + '- No ENSO'
	filename = './figures/hgt_200_composites_diff_PV_' + seas[i] +'_noENSO.png'
	plots.PlotCompDiff(var, ds.latitude, ds.longitude, tit, filename)


