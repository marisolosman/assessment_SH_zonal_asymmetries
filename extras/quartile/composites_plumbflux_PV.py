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
PV = xr.open_dataset(RUTA + 'PV_index.nc')

month = ['Aug', 'Sep', 'Oct', 'Nov']
seas = ['ASO', 'SON']
index_monthly_upper = PV.PV_mon >= PV.PV_mon.quantile(0.75, dim='dim_0', interpolation='linear')
index_monthly_lower = PV.PV_mon <= PV.PV_mon.quantile(0.25, dim='dim_0', interpolation='linear')

index_monthly_normal = np.logical_and(PV.PV_mon < PV.PV_mon.quantile(0.75, dim='dim_0', interpolation='linear'), PV.PV_mon > PV.PV_mon.quantile(0.25, dim='dim_0', interpolation='linear'))

for i in np.arange(0,4):
	var_neg = np.mean(ds.z.values[i, index_monthly_upper.values, :, :], axis=0) - np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0)
	var_pos = np.mean(ds.z.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0)
	px_neg, py_neg, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], np.mean(ds.z.values[i, index_monthly_upper.values, :, :], axis=0), np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0),ds.latitude.values, ds.longitude.values )
	px_pos, py_pos, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], np.mean(ds.z.values[i, index_monthly_lower.values, :, :], axis=0), np.mean(ds.z.values[i, index_monthly_normal.values, :, :], axis=0),ds.latitude.values, ds.longitude.values )
	tit = 'Composites S4 Plumb Fluxes 200hPa - ' + month[i]
	filename = './figures/hgt_200_plumb_composites_PV_' + month[i] +'.png'
	plots.PlotCompositesPlumbPV(var_pos, var_neg,px_pos, px_neg, py_pos, py_neg, ds.latitude,
			    ds.longitude, tit, filename)
	var = np.mean(ds.z.values[i, index_monthly_lower.values, :, :], axis=0) - np.mean(ds.z.values[i, index_monthly_upper.values, :, :], axis=0)
	px_c, py_c, lat = plumb_flux.ComputePlumbFluxes(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var, np.zeros_like(var), ds.latitude.values, ds.longitude.values)
	tit = 'Composites differences SPV-WPV Years - ' + month[i]
	filename = './figures/hgt_200_plumb_composites_diff_PV_' + month[i] +'.png'
	plots.PlotCompPlumbDiff(var, px_c, py_c, ds.latitude, ds.longitude, tit, filename)

for i in np.arange(0,2):
	var = ds.isel(month=range(i, i +3 )).mean(dim='month')
	var_wnd = winds_clm.isel(month=range(i, i +3 )).mean(dim='month')
	
	var_neg = np.mean(var.z.values[index_monthly_upper.values, :, :], axis=0) - np.mean(var.z.values[index_monthly_normal.values, :, :], axis=0)
	var_pos = np.mean(var.z.values[index_monthly_lower.values, :, :], axis=0) - np.mean(var.z.values[index_monthly_normal.values, :, :], axis=0)
	px_neg, py_neg, lat = plumb_flux.ComputePlumbFluxes(var_wnd.u.values[:, :], var_wnd.v.values[:, :], np.mean(var.z.values[index_monthly_upper.values, :, :], axis=0), np.mean(var.z.values[index_monthly_normal.values, :, :], axis=0),ds.latitude.values, ds.longitude.values)
	px_pos, py_pos, lat = plumb_flux.ComputePlumbFluxes(var_wnd.u.values[:, :], var_wnd.v.values[:, :], np.mean(var.z.values[index_monthly_lower.values, :, :], axis=0), np.mean(var.z.values[index_monthly_normal.values, :, :], axis=0),ds.latitude.values, ds.longitude.values )
	tit = 'Composites S4 Plum Fluxes 200hPa - ' + seas[i]
	filename = './figures/hgt_200_plumb_composites_PV_' + seas[i] +'.png'
	plots.PlotCompositesPlumbPV(var_pos, var_neg,px_pos, px_neg, py_pos, py_neg, ds.latitude,
				 ds.longitude, tit, filename)
	var = np.mean(var.z.values[index_monthly_lower.values, :, :], axis=0) - np.mean(var.z.values[index_monthly_upper.values, :, :], axis=0)
	px_c, py_c, lat = plumb_flux.ComputePlumbFluxes(var_wnd.u.values[:, :], var_wnd.v.values[:, :], var, np.zeros_like(var), ds.latitude.values, ds.longitude.values)
	tit = 'Composites differences SPV-WPV Years - ' + seas[i]
	filename = './figures/hgt_200_plumb_composites_diff_PV_' + seas[i] +'.png'
	plots.PlotCompPlumbDiff(var, px_c, py_c, ds.latitude, ds.longitude, tit, filename)

