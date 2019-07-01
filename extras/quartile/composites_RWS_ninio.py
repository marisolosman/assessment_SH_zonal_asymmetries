#sort el ninio events
import numpy as np
import xarray as xr
import pandas as pd
import plots
import os
from windspharm.xarray import VectorWind
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '/storage/shared/glusterfs/acrcc/users/vg140344/data/'
#abro el archivo de geopotencial y junto la coordenada year y numbre
ninio34 = xr.open_dataset(RUTA + 'ninio34_index.nc')
ds = xr.open_dataset(RUTA + 'winds200.nc', chunks={'latitude':10})
ds = ds.stack(realiz=['year', 'number'])
ds = ds.transpose('month', 'realiz', 'latitude', 'longitude')
ds = ds.drop(['IC', 'isobaricInhPa'])
ds = ds.reset_index('realiz')
month = ['Aug', 'Sep', 'Oct', 'Nov']
seas = ['ASO', 'SON']
index_monthly_upper = ninio34.ninio34_mon >= ninio34.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34.ninio34_mon <= ninio34.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear')
index_monthly_normal = np.logical_and(ninio34.ninio34_mon < ninio34.ninio34_mon.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34.ninio34_mon > ninio34.ninio34_mon.quantile(0.25, dim='dim_0', interpolation='linear'))

for i in np.arange(0,4):
	#create vectorwind instance
	w = VectorWind(ds['u'][i, :, :, :], ds['v'][i, :, :, :])
	eta = w.absolutevorticity()
	div = w.divergence()
	uchi, vchi = w.irrotationalcomponent()
	etax, etay = w.gradient(eta)
	etax.attrs['units'] = 'm**-1 s**-1'
	etay.attrs['units'] = 'm**-1 s**-1'
	# Combine the components to form the Rossby wave source term.
	S = eta * -1. * div - (uchi * etax + vchi * etay)
	S *= 1e11 
	var_neg = np.mean(S[index_monthly_lower.values, :, :], axis=0) - np.mean(S[index_monthly_normal.values, :, :], axis=0)
	var_pos = np.mean(S[index_monthly_upper.values, :, :], axis=0) - np.mean(S[index_monthly_normal.values, :, :], axis=0)
	uchi_neg = np.mean(uchi[index_monthly_lower.values, :, :], axis=0) - np.mean(uchi[index_monthly_normal.values, :, :], axis=0)
	uchi_pos = np.mean(uchi[index_monthly_upper.values, :, :], axis=0) - np.mean(uchi[index_monthly_normal.values, :, :], axis=0)
	vchi_neg = np.mean(vchi[index_monthly_lower.values, :, :], axis=0) - np.mean(vchi[index_monthly_normal.values, :, :], axis=0)
	vchi_pos = np.mean(vchi[index_monthly_upper.values, :, :], axis=0) - np.mean(vchi[index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 RWS and Divergent wind 200hPa - ' + month[i]
	filename = './figures/200_RWS_chi_composites_NINIO_' + month[i] +'.png'
	plots.PlotCompositesDivPlumb(var_pos, var_neg, uchi_pos, uchi_neg, vchi_pos, vchi_neg,
				     S.latitude.values, S.longitude.values, tit, filename)
	var = np.mean(S[index_monthly_upper.values, :, :], axis=0) - np.mean(S[index_monthly_lower.values, :, :], axis=0)
	uchi_1 = np.mean(uchi[index_monthly_upper.values, :, :], axis=0) - np.mean(uchi[index_monthly_lower.values, :, :], axis=0)
	vchi_1 = np.mean(vchi[index_monthly_upper.values, :, :], axis=0) - np.mean(vchi[index_monthly_lower.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa differences EN-LN Years - ' + month[i]
	filename = './figures/200_RWS_chi_composites_diff_NINIO_' + month[i] +'.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)
for i in np.arange(0,2):
	ds1 = ds.isel(month=range(i, i +3)).mean(dim='month')
	#create vectorwind instance
	w = VectorWind(ds1['u'][:, :, :], ds1['v'][:, :, :])
	eta = w.absolutevorticity()
	div = w.divergence()
	uchi, vchi = w.irrotationalcomponent()
	etax, etay = w.gradient(eta)
	etax.attrs['units'] = 'm**-1 s**-1'
	etay.attrs['units'] = 'm**-1 s**-1'
#	# Combine the components to form the Rossby wave source term.
	S = eta * -1. * div - (uchi * etax + vchi * etay)
	S *= 1e11 
	var_neg = np.mean(S[index_monthly_lower.values, :, :], axis=0) - np.mean(S[index_monthly_normal.values, :, :], axis=0)
	var_pos = np.mean(S[index_monthly_upper.values, :, :], axis=0) - np.mean(S[index_monthly_normal.values, :, :], axis=0)
	uchi_neg = np.mean(uchi[index_monthly_lower.values, :, :], axis=0) - np.mean(uchi[index_monthly_normal.values, :, :], axis=0)
	uchi_pos = np.mean(uchi[index_monthly_upper.values, :, :], axis=0) - np.mean(uchi[index_monthly_normal.values, :, :], axis=0)
	vchi_neg = np.mean(vchi[index_monthly_lower.values, :, :], axis=0) - np.mean(vchi[index_monthly_normal.values, :, :], axis=0)
	vchi_pos = np.mean(vchi[index_monthly_upper.values, :, :], axis=0) - np.mean(vchi[index_monthly_normal.values, :, :], axis=0)
	tit = 'Composites S4 RWS and Divergent wind 200hPa - ' + seas[i]
	filename = './figures/200_RWS_chi_composites_NINIO_' + seas[i] +'.png'
	plots.PlotCompositesDivPlumb(var_pos, var_neg, uchi_pos, uchi_neg, vchi_pos, vchi_neg,
				     S.latitude.values, S.longitude.values, tit, filename)
	var = np.mean(S[index_monthly_upper.values, :, :], axis=0) - np.mean(S[index_monthly_lower.values, :, :], axis=0)
	uchi_1 = np.mean(uchi[index_monthly_upper.values, :, :], axis=0) - np.mean(uchi[index_monthly_lower.values, :, :], axis=0)
	vchi_1 = np.mean(vchi[index_monthly_upper.values, :, :], axis=0) - np.mean(vchi[index_monthly_lower.values, :, :], axis=0)
	tit = 'Composites S4 RWS and divergent wind 200 hPa differences EN-LN Years - ' + seas[i]
	filename = './figures/200_RWS_chi_composites_diff_NINIO_' + seas[i] +'.png'
	plots.PlotCompDivPlumbDiff(var, uchi_1, vchi_1, S.latitude, S.longitude, tit, filename)

