#concatenate forecastr of ninio indexes made with aug ic fo aug to feb
import numpy as np
import xarray as xr
import pandas as pd
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA='/home/users/vg140344/datos/data/fogt/'
##abro el archivo de sst en ASON y junto la coordenada year y number
ds = xr.open_dataset(RUTA + 'sst_ninio34_s4.nc4')
ds.coords['year'] = np.arange(1981, 2017)
ds = ds.stack(realiz = ['year', 'number'])
ds = ds.sel(**{'month':slice(8,11)})
ds = ds.reset_index('realiz')
ds = ds.drop(['year', 'number', 'IC'])
#abro el archivo de sst en DEF y junto la coordenada year and number
ds1 = xr.open_dataset(RUTA + 'sst_ninio34_s4_dec_feb.nc4')
ds1.coords['year'] = np.arange(1981, 2017)
ds1 = ds1.stack(realiz = ['year', 'number'])
ds1 = ds1.sel(**{'month': [1, 2, 12]})
ds1 = ds1.reset_index('realiz')
ds1 = ds1.drop([ 'IC', 'year', 'number'])
#concatenate and sort months
ds2 = xr.concat([ds, ds1], dim='month')
ds2 = ds2.loc[{'month': [8, 9, 10, 11, 12, 1, 2]}]
ds2.to_netcdf(RUTA + 'sst_ninio34_aug_feb.nc4')

##abro el archivo de sst en ASON y junto la coordenada year y number
ds = xr.open_dataset(RUTA + 'sst_ninio3_s4.nc4')
ds.coords['year'] = np.arange(1981, 2017)
ds = ds.stack(realiz = ['year', 'number'])
ds = ds.sel(**{'month':slice(8,11)})
ds = ds.reset_index('realiz')
ds = ds.drop(['year', 'number', 'IC'])
#abro el archivo de sst en DEF y junto la coordenada year and number
ds1 = xr.open_dataset(RUTA + 'sst_ninio3_s4_dec_feb.nc4')
ds1.coords['year'] = np.arange(1981, 2017)
ds1 = ds1.stack(realiz = ['year', 'number'])
ds1 = ds1.sel(**{'month': [1, 2, 12]})
ds1 = ds1.reset_index('realiz')
ds1 = ds1.drop([ 'IC', 'year', 'number'])
#concatenate and sort months
ds2 = xr.concat([ds, ds1], dim='month')
ds2 = ds2.loc[{'month': [8, 9, 10, 11, 12, 1, 2]}]
ds2.to_netcdf(RUTA + 'sst_ninio3_aug_feb.nc4')

##abro el archivo de sst en ASON y junto la coordenada year y number
ds = xr.open_dataset(RUTA + 'sst_ninio4_s4.nc4')
ds.coords['year'] = np.arange(1981, 2017)
ds = ds.stack(realiz = ['year', 'number'])
ds = ds.sel(**{'month':slice(8,11)})
ds = ds.reset_index('realiz')
ds = ds.drop(['year', 'number', 'IC'])
#abro el archivo de sst en DEF y junto la coordenada year and number
ds1 = xr.open_dataset(RUTA + 'sst_ninio4_s4_dec_feb.nc4')
ds1.coords['year'] = np.arange(1981, 2017)
ds1 = ds1.stack(realiz = ['year', 'number'])
ds1 = ds1.sel(**{'month': [1, 2, 12]})
ds1 = ds1.reset_index('realiz')
ds1 = ds1.drop([ 'IC', 'year', 'number'])
#concatenate and sort months
ds2 = xr.concat([ds, ds1], dim='month')
ds2 = ds2.loc[{'month': [8, 9, 10, 11, 12, 1, 2]}]
ds2.to_netcdf(RUTA + 'sst_ninio4_aug_feb.nc4')

