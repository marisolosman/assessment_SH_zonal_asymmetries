#sort el ninio events
import numpy as np
import xarray as xr
import pandas as pd
import plumb_flux
import os

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA='/storage/shared/glusterfs/acrcc/users/vg140344/data/'
#abro el archivo de geopotencial y junto la coordenada year y numbre
winds = xr.open_dataset(RUTA + 'winds200.nc')
winds.coords['year'] = np.arange(1981, 2017)
winds = winds.stack(realiz = ['year', 'number'])
winds = winds.transpose('month', 'realiz', 'latitude', 'longitude')
hgt = xr.open_dataset(RUTA + 'hgt200.nc')
hgt.coords['year'] = np.arange(1981, 2017)
hgt = hgt.stack(realiz = ['year', 'number'])
hgt = hgt.transpose('month', 'realiz', 'latitude', 'longitude')
hgt.z.values = hgt.z.values / 9.8
month = ['Aug', 'Sep', 'Oct', 'Nov']
seas = ['ASO', 'SON']
px = []
py = []
for i in np.arange(0,4):
	pfx, pfy, lat = plumb_flux.ComputePlumbFluxes(winds.u.values[i, :, :, :],
						 winds.v.values[i,:,:,:], hgt.z.values[i, :, :, :],
						 hgt.latitude.values, hgt.longitude.values)
	pfx = xr.DataArray(pfx, coords=[np.arange(51 * 36), lat,
					hgt.longitude.values], dims=['realiz', 'latitude', 'logitude'])
	pfy = xr.DataArray(pfy, coords=[np.arange(51 * 36), lat,
					hgt.longitude.values], dims=['realiz', 'latitude', 'logitude'])
	px.append(pfx)
	py.append(pfy)

px = xr.concat(px, dim='month')
px ['month'] = np.array([8, 9, 10, 11])
py = xr.concat(py, dim='month')
py ['month'] = np.array([8, 9, 10, 11])

px.to_netcdf(RUTA + 'monthly_plumb_xflux.nc')
py.to_netcdf(RUTA + 'monthly_plumb_yflux.nc')


px = []
py = []		
for i in np.arange(0,2):
	if i == 0:
		var1 = hgt.sel(**{'month':slice(8, 10)}).mean(dim='month')
		var2 = winds.sel(**{'month':slice(8, 10)}).mean(dim='month')

	else:
		var1 = hgt.sel(**{'month':slice(9, 11)}).mean(dim='month')
		var2 = winds.sel(**{'month':slice(9, 11)}).mean(dim='month')
	pfx, pfy, lat = plumb_flux.ComputePlumbFluxes(var2.u.values, var2.v.values,
						      var1.z.values,
						      var1.latitude.values, var1.longitude.values)
	pfx = xr.DataArray(pfx, coords=[np.arange(51 * 36), lat,
					hgt.longitude.values], dims=['realiz', 'latitude', 'logitude'])
	pfy = xr.DataArray(pfy, coords=[np.arange(51 * 36), lat,
					hgt.longitude.values], dims=['realiz', 'latitude', 'logitude'])
	px.append(pfx)
	py.append(pfy)

px = xr.concat(px, dim='season')
px ['season'] = np.array(['ASO', 'SON'])
py = xr.concat(py, dim='season')
py ['season'] = np.array(['ASO', 'SON'])

px.to_netcdf(RUTA + 'seasonal_plumb_xflux.nc')
py.to_netcdf(RUTA + 'seasonal_plumb_yflux.nc')


