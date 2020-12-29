# concate monthly wind data
import numpy as np
import xarray as xr
import pandas as pd
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/data/'
#abro el archivo de geopotencial en ASON y junto la coordenada year y number
ds = xr.open_dataset(RUTA + 'winds200.nc')
ds.coords['year'] = np.arange(1981, 2017)
ds = ds.stack(realiz = ['year', 'number'])
ds = ds.transpose('month', 'realiz', 'latitude', 'longitude')
ds = ds.sel(**{'month': slice(8, 11)})
ds = ds.drop(['isobaricInhPa', 'IC'])
ds = ds.reset_index('realiz')
ds = ds.drop(['year', 'number'])
# abro el archivo de geopotencial en DEF y junto la coordenada year and number
ds1 = xr.open_dataset(RUTA + 'winds200_djf.nc')
ds1.coords['year'] = np.arange(1981, 2017)
ds1 = ds1.stack(realiz = ['year', 'number'])
ds1 = ds1.transpose('month', 'realiz', 'latitude', 'longitude')
ds1 = ds1.sel(**{'month': [1, 2, 12]})
ds1 = ds1.drop(['isobaricInhPa', 'IC', 'year', 'number'])
ds1 = ds1.reset_index('realiz')
ds1.coords['latitude'] = ds.latitude.values
ds1.coords['longitude'] = ds.longitude.values

ds2 = xr.concat([ds, ds1], dim='month')
ds2 = ds2.loc[{'month': [8, 9, 10, 11, 12, 1, 2]}]
ds2.to_netcdf('~/datos/data/monthly_winds200_aug_feb.nc')

# abro datos de 81-92
ds = xr.open_dataset(RUTA + 'winds50_asondjf_81-92.nc4', chunks={'latitude':10} )
ds.coords['year'] = np.arange(1981, 1993)
ds = ds.sel(**{'month': [8, 9, 10, 11, 12, 1, 2]})
ds = ds.transpose('month', 'year', 'number', 'latitude', 'longitude')
ds = ds.drop(['isobaricInhPa', 'IC'])

#abro el archivo de 93-17 de ASON
ds1 = xr.open_dataset(RUTA + 'winds50_ason_93-17.nc4', chunks={'latitude':10})
ds1.coords['year'] = np.arange(1993, 2017)
ds1 = ds1.transpose('month', 'year', 'number', 'latitude', 'longitude')
ds1 = ds1.sel(**{'month': [8, 9, 10, 11]})
ds1.coords['latitude'] = ds.latitude.values
ds1.coords['longitude'] = ds.longitude.values
ds1 = ds1.drop(['isobaricInhPa', 'IC'])

#abro el archivo de 93-17 DEF
ds2 = xr.open_dataset(RUTA + 'winds50_djf_93-17.nc4', chunks={'latitude':10})
ds2.coords['year'] = np.arange(1993, 2017)
ds2 = ds2.transpose('month', 'year', 'number', 'latitude', 'longitude')
ds2 = ds2.sel(**{'month': [12, 1, 2]})
ds2.coords['latitude'] = ds.latitude.values
ds2.coords['longitude'] = ds.longitude.values
ds2 = ds2.drop(['isobaricInhPa', 'IC'])

#junto todo
ds3 = xr.concat([ds1, ds2], dim='month')
ds4 = xr.concat([ds, ds3], dim='year')
ds4 = ds4.stack(realiz=['year', 'number'])
ds4 = ds4.reset_index('realiz')
ds4 = ds4.loc[{'month': [8, 9, 10, 11, 12, 1, 2]}]
ds4.to_netcdf('~/datos/data/monthly_winds50_aug_feb.nc4')

