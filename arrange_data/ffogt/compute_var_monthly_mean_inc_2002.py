#compute monthly geopotential height for different levels and IC including 2002
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

RUTA = '/storage/silver/acrcc/ys916780/Data_for_Sol/output'
hgt = []
for Y in np.arange(1981,2018):
	ds = xr.open_dataset(RUTA + '_' + str(Y)+ '_geopotential200', engine='cfgrib', chunks={'step': 10},
			     backend_kwargs={'indexpath': ''})
	ds = ds.rename({'time': 'IC'}).set_coords(['IC'])
	ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
	ds = ds.rename({'step': 'time'}).set_coords(['time'])
	ds = ds.groupby('time.month').mean(dim='time').compute()
	hgt.append(ds)
hgt =xr.concat(hgt, dim='year')
hgt.to_netcdf('/storage/silver/acrcc/vg140344/data/hgt200_inc_2002.nc')

hgt = []
for Y in np.arange(1981,2018):
	ds = xr.open_dataset(RUTA + '_' + str(Y)+ '_geopotential50', engine='cfgrib', chunks={'step': 10},
			     backend_kwargs={'indexpath': ''})
	ds = ds.rename({'time': 'IC'}).set_coords(['IC'])
	ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
	ds = ds.rename({'step': 'time'}).set_coords(['time'])
	ds = ds.groupby('time.month').mean(dim='time').compute()
	hgt.append(ds)

hgt =xr.concat(hgt, dim='year')

hgt.compute().to_netcdf('/storage/silver/acrcc/vg140344/data/hgt50_inc_2002.nc')

RUTA = '/storage/silver/acrcc/vg140344/HGT_S4Hindcasts_0108_2levels_24Hourly'
hgt = []
for Y in np.arange(1981,2018):
	ds = xr.open_dataset(RUTA + str(Y)+ '0801.grib', engine='cfgrib',chunks={'step': 10},
			     backend_kwargs={'indexpath': ''})
	ds = ds.rename({'time': 'IC'}).set_coords(['IC'])
	ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
	ds = ds.rename({'step': 'time'}).set_coords(['time'])
	ds = ds.groupby('time.month').mean(dim='time').compute()
	hgt.append(ds)
hgt =xr.concat(hgt, dim='year')
hgt.to_netcdf('/storage/silver/acrcc/vg140344/data/hgt_levels_0801_inc_2002.nc')

