#compute monthly el ninio 3.4 index
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

RUTA = '/storage/shared/glusterfs/acrcc/users/ys916780/Data_for_Sol/output'
#hgt = []
#for Y in np.arange(1981,2018):
#	if Y != 2002:
#		ds = xr.open_dataset(RUTA + '_' + str(Y)+ '_geopotential200', engine='cfgrib', chunks={'step':10},  backend_kwargs={'indexpath':''})
#		ds = ds.rename({'time':'IC'}).set_coords(['IC'])
#		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
#		ds = ds.rename({'step':'time'}).set_coords(['time'])
#		ds = ds.groupby('time.month').mean(dim='time')
#		hgt.append(ds)
#
#hgt =xr.concat(hgt, dim='year')
#
#hgt.to_netcdf('/storage/shared/glusterfs/acrcc/users/vg140344/data/hgt200.nc')
hgt = []
for Y in np.arange(1981,2018):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + '_' + str(Y)+ '_geopotential50', engine='cfgrib', chunks={'step':10},  backend_kwargs={'indexpath':''})
		ds = ds.rename({'time':'IC'}).set_coords(['IC'])
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step':'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		hgt.append(ds)

hgt =xr.concat(hgt, dim='year')

hgt.compute().to_netcdf('/storage/shared/glusterfs/acrcc/users/vg140344/data/hgt50.nc')

RUTA = '/storage/shared/glusterfs/acrcc/users/vg140344/HGT_S4Hindcasts_0111_200hPa_24Hourly'
hgt = []
for Y in np.arange(1981,2018):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + str(Y)+ '1101.nc', chunks={'time':10})
		ds = ds.groupby('time.month').mean(dim='time')
		hgt.append(ds)
#
hgt =xr.concat(hgt, dim='time')

hgt.compute().to_netcdf('/storage/shared/glusterfs/acrcc/users/vg140344/data/hgt_djf_200.nc')
#RUTA = '/storage/shared/glusterfs/acrcc/users/vg140344/HGT_S4Hindcasts_0108_2levels_24Hourly'
#hgt = []
#for Y in np.arange(1981,2018):
#	if Y != 2002:
#		ds = xr.open_dataset(RUTA + str(Y)+ '0801.grib', engine='cfgrib',chunks={'step':10}, backend_kwargs={'indexpath':''})
#		ds = ds.rename({'time':'IC'}).set_coords(['IC'])
#		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
#		ds = ds.rename({'step':'time'}).set_coords(['time'])
#		ds = ds.groupby('time.month').mean(dim='time')
#		hgt.append(ds)
#
#hgt =xr.concat(hgt, dim='year')
#
#hgt.to_netcdf('/storage/shared/glusterfs/acrcc/users/vg140344/data/hgt_levels_0801.nc')

