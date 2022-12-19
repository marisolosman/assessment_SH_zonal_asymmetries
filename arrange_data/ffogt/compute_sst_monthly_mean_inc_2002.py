#compute monthly sst for S4 forecasts
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import os

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

RUTA = '/storage/silver/acrcc/ys916780/Data_for_Sol/output'

sst = []
for Y in np.arange(1981,2018):
	ds = xr.open_dataset(RUTA + '_' + str(Y) + '_sst', engine='cfgrib', chunks={'step':10},
			     backend_kwargs={'indexpath':''})
	ds = ds.rename({'time':'IC'}).set_coords(['IC'])
	ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
	ds = ds.rename({'step':'time'}).set_coords(['time'])
	ds = ds.groupby('time.month').mean(dim='time')
	sst.append(ds)
#
sst =xr.concat(sst, dim='year')
sst.to_netcdf('~/datos/data/fogt/sst_s4_aug_nov_inc_2002.nc4')

RUTA = '/storage/silver/acrcc/vg140344/S4Hindcast/'
sst = []
for Y in np.arange(1981,2018):
	ds = xr.open_dataset(RUTA + 'SST_S4Hindcast_0108_24Hourly_' + str(Y) + '0801.grib',
				 engine='cfgrib', chunks={'step':10},
			     backend_kwargs={'indexpath':''})
	ds = ds.rename({'time':'IC'}).set_coords(['IC'])
	ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
	ds = ds.rename({'step':'time'}).set_coords(['time'])
	ds = ds.groupby('time.month').mean(dim='time')
	sst.append(ds)

sst =xr.concat(sst, dim='year')
sst.to_netcdf('/home/users/vg140344/datos/data/fogt/sst_s4_dec_feb_inc_2002.nc4')

