#sort el ninio events
import numpy as np
import xarray as xr
import pandas as pd
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA='~/datos/data/'
##abro el archivo de geopotencial en ASON y junto la coordenada year y number
ds = xr.open_dataset(RUTA + 'hgt200.nc')
ds.coords['year'] = np.arange(1981, 2017)
ds = ds.stack(realiz = ['year', 'number'])
ds = ds.transpose('month', 'realiz', 'latitude', 'longitude')
ds.z.values = ds.z.values / 10
ds = ds.sel(**{'month':slice(8,11)})
ds = ds.drop(['isobaricInhPa', 'IC'])
ds = ds.reset_index('realiz')
ds = ds.drop(['year'])
#abro el archivo de geopotencial en DEF y junto la coordenada year and number
ds1 = xr.open_dataset(RUTA + 'hgt_djf_200.nc')
ds1.coords['time'] = np.arange(1981, 2017)
ds1 = ds1.stack(realiz = ['time', 'number'])
ds1 = ds1.transpose('month', 'realiz', 'latitude', 'longitude')
ds1.z.values = ds1.z.values / 10
ds1 = ds1.sel(**{'month': [1, 2, 12]})
ds1 = ds1.reset_index('realiz')
ds1.coords['latitude'] = ds.latitude.values
ds1.coords['longitude'] = ds.longitude.values
ds1 = ds1.drop(['time'])

ds2 = xr.concat([ds, ds1], dim='month')
ds2 = ds2.loc[{'month': [8, 9, 10, 11, 12, 1, 2]}]
ds2.to_netcdf('~/datos/data/monthly_hgt200_aug_feb.nc')
#
#os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
#RUTA='/storage/shared/glusterfs/acrcc/users/vg140344/data/'
##abro el archivo de geopotencial en ASON y junto la coordenada year y number
#ds = xr.open_dataset(RUTA + 'hgt50.nc')
#ds.coords['year'] = np.arange(1981, 2017)
#ds = ds.stack(realiz = ['year', 'number'])
#ds = ds.transpose('month', 'realiz', 'latitude', 'longitude')
#ds.z.values = ds.z.values / 10
#ds = ds.sel(**{'month':slice(8,11)})
#ds = ds.drop(['isobaricInhPa', 'IC'])
#ds = ds.reset_index('realiz')
##abro el archivo de geopotencial en DEF y junto la coordenada year and number
#ds1 = xr.open_dataset(RUTA + 'hgt_levels_0801.nc')
##selecciono el nivel de 50
#ds1 = ds1.isel(isobaricInhPa=1)
#ds = ds.drop(['isobaricInhPa', 'IC', 'time'])
#ds1.coords['year'] = np.arange(1981, 2017)
#ds1 = ds1.stack(realiz = ['year', 'number'])
#ds1 = ds1.transpose('month', 'realiz', 'latitude', 'longitude')
#ds1.z.values = ds1.z.values / 10
#ds1 = ds1.sel(**{'month': [1, 2, 12]})
#ds1 = ds1.reset_index('realiz')
#ds1.coords['latitude'] = ds.latitude.values
#ds1.coords['longitude'] = ds.longitude.values
#
#ds2 = xr.concat([ds, ds1], dim='month')
#ds2 = ds2.loc[{'month': [8, 9, 10, 11, 12, 1, 2]}]
#ds2.to_netcdf(RUTA + 'monthly_hgt50_aug_feb.nc')
#

