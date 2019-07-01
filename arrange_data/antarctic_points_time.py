#select AS_point (67.5 and 120W) and eastern antarctic point (67.5 and 110E) and save time series
import numpy as np
import xarray as xr
import pandas as pd
import os

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

RUTA = '~/datos/data/'

#open files

hgt200 = xr.open_dataset(RUTA + 'monthly_hgt200_aug_feb.nc')
hgt50 = xr.open_dataset(RUTA + 'monthly_hgt50_aug_feb.nc')
hgt50 = hgt50.drop(['year'])
ds = xr.concat([hgt200, hgt50], dim='level')
ds.coords['level'] = np.array([200, 50])

ds_as = ds.sel(latitude=-67.5, longitude=240, method='nearest')
ds_ea = ds.sel(latitude=-67.5, longitude=110, method='nearest')
ds1 = xr.concat([ds_ea, ds_as], dim='point')
ds1.coords['point'] = ['EA', 'AS']
ds1.to_netcdf('~/datos/data/hgt_points.nc')

#hgt50_as = hgt50.sel(latitude=-67.5, longitude=240, method='nearest')
#hgt50_ea = hgt50.sel(latitude=-67.5, longitude=110, method='nearest')
