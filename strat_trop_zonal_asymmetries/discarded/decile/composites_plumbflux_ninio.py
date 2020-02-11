#sort el ninio events
import numpy as np
import xarray as xr
import pandas as pd
import plots
import os
import plumb_flux
RUTA='~/datos/data/'
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
winds = xr.open_dataset(RUTA + 'winds200.nc')
winds.coords['year'] = np.arange(1981, 2017)
winds = winds.stack(realiz = ['year', 'number'])
winds = winds.transpose('month', 'realiz', 'latitude', 'longitude')
winds_clm = winds.mean(dim='realiz')
#abro el archivo de geopotencial y junto la coordenada year y numbre
ds = xr.open_dataset(RUTA + 'hgt200.nc')
ds.coords['year'] = np.arange(1981, 2017)
ds = ds.stack(realiz = ['year', 'number'])
ds = ds.transpose('month', 'realiz', 'latitude', 'longitude')
ds.z.values = ds.z.values / 10
ninio34 = xr.open_dataset(RUTA + 'ninio34_index.nc')
px = xr.open_dataset(RUTA + 'monthly_plumb_xflux.nc')
px = px.__xarray_dataarray_variable__.rename('pfx')
py = xr.open_dataset(RUTA + 'monthly_plumb_yflux.nc')
py = py.__xarray_dataarray_variable__.rename('pfy')

month = ['Aug', 'Sep', 'Oct', 'Nov']
seas = ['ASO', 'SON']
index_monthly_upper = ninio34.ninio34_mon >= ninio34.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear')
index_monthly_lower = ninio34.ninio34_mon <= ninio34.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear')

index_monthly_normal = np.logical_and(ninio34.ninio34_mon < ninio34.ninio34_mon.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34.ninio34_mon > ninio34.ninio34_mon.quantile(0.10, dim='dim_0', interpolation='linear'))

for i in np.arange(0,4):
	var_neg = np.mean(ds.z.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0)
	var_pos = np.mean(ds.z.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0)
	px_neg, py_neg, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], np.mean(ds.z.values[i, index_monthly_lower.values, :, :], axis=0), np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0),ds.latitude.values, ds.longitude.values )
	px_pos, py_pos, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], np.mean(ds.z.values[i, index_monthly_upper.values, :, :], axis=0), np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0),ds.latitude.values, ds.longitude.values )
	tit = 'Composites S4 Plumb Fluxes 200hPa - ' + month[i]
	filename = './figures_decile/hgt_200_plumb_composites_NINIO_' + month[i] +'.png'
	plots.PlotCompositesPlumb(var_pos, var_neg,px_pos, px_neg, py_pos, py_neg, ds.latitude,
			    ds.longitude, tit, filename)
	var = np.mean(ds.z.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(ds.z.values[i, index_monthly_lower.values, :, :], axis=0)
	px_c, py_c, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var, np.zeros_like(var), ds.latitude.values, ds.longitude.values)
	tit = 'Composites differences EN-LN Years - ' + month[i]
	filename = './figures_decile/hgt_200_plumb_composites_diff_NINIO_' + month[i] +'.png'
	plots.PlotCompPlumbDiff(var, px_c, py_c, ds.latitude, ds.longitude, tit, filename)
#
px = xr.open_dataset(RUTA + 'seasonal_plumb_xflux.nc')
px = px.__xarray_dataarray_variable__.rename('pfx')
py = xr.open_dataset(RUTA + 'seasonal_plumb_yflux.nc')
py = py.__xarray_dataarray_variable__.rename('pfy')


for i in np.arange(0,2):
	var = ds.isel(month=range(i, i +3 )).mean(dim='month')
	var_wnd = winds_clm.isel(month=range(i, i +3 )).mean(dim='month')
	
	var_neg = np.mean(var.z.values[index_monthly_lower.values, :, :], axis=0) - np.mean(var.z.values[index_monthly_normal.values, :, :], axis=0)
	var_pos = np.mean(var.z.values[index_monthly_upper.values, :, :], axis=0) - np.mean(var.z.values[index_monthly_normal.values, :, :], axis=0)
	px_neg, py_neg, lat = plumb_flux.ComputePlumbFluxes(var_wnd.u.values[:, :], var_wnd.v.values[:, :], np.mean(var.z.values[index_monthly_lower.values, :, :], axis=0), np.mean(var.z.values[index_monthly_normal.values, :, :], axis=0),ds.latitude.values, ds.longitude.values)
	px_pos, py_pos, lat = plumb_flux.ComputePlumbFluxes(var_wnd.u.values[:, :], var_wnd.v.values[:, :], np.mean(var.z.values[index_monthly_upper.values, :, :], axis=0), np.mean(var.z.values[index_monthly_normal.values, :, :], axis=0),ds.latitude.values, ds.longitude.values )
	tit = 'Composites S4 Plum Fluxes 200hPa - ' + seas[i]
	filename = './figures_decile/hgt_200_plumb_composites_NINIO_' + seas[i] +'.png'
	plots.PlotCompositesPlumb(var_pos, var_neg,px_pos, px_neg, py_pos, py_neg, ds.latitude,
				 ds.longitude, tit, filename)
	var = np.mean(var.z.values[index_monthly_upper.values, :, :], axis=0) - np.mean(var.z.values[index_monthly_lower.values, :, :], axis=0)
	px_c, py_c, lat = plumb_flux.ComputePlumbFluxes(var_wnd.u.values[:, :], var_wnd.v.values[:, :], var, np.zeros_like(var), ds.latitude.values, ds.longitude.values)
	tit = 'Composites differences EN-LN Years - ' + seas[i]
	filename = './figures_decile/hgt_200_plumb_composites_diff_NINIO_' + seas[i] +'.png'
	plots.PlotCompPlumbDiff(var, px_c, py_c, ds.latitude, ds.longitude, tit, filename)

