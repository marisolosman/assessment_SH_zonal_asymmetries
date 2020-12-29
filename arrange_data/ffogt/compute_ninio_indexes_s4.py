#compute monthly el ninio indexes for S4 forecasts
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import os

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '/storage/silver/acrcc/vg140344/'

sst34 = []
sst3 = [] #ninio 3 5N-5S 150W-90W 
sst4 = [] #ninio 4 5N-5S 160E-150W
for Y in np.arange(1981,2018):
	if Y != 2002:
		ds = xr.open_dataset(RUTA + 'SST_S4Hindcast_0108_24Hourly_' + str(Y) + '0801.grib',
					 engine='cfgrib', chunks={'step': 10},
				     backend_kwargs={'indexpath': ''})
		ds = ds.rename({'time': 'IC'}).set_coords(['IC'])
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step': 'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		ds_34 = ds.sel(**{'latitude': slice(5, -5), 'longitude': slice(190, 240)}).mean(
			dim=['longitude', 'latitude']).compute()
		ds_3 = ds.sel(**{'latitude': slice(5, -5), 'longitude': slice(210, 270)}).mean(
			dim=['longitude', 'latitude']).compute()
		ds_4 = ds.sel(**{'latitude': slice(5, -5), 'longitude': slice(160, 210)}).mean(
			dim=['longitude', 'latitude']).compute()
		sst4.append(ds_4)
		sst3.append(ds_3)
		sst34.append(ds_34)

sst34 =xr.concat(sst34, dim='year')
sst3 =xr.concat(sst3, dim='year')
sst4 =xr.concat(sst4, dim='year')

sst34.to_netcdf('/home/users/vg140344/datos/data/fogt/sst_ninio34_s4_dec_feb.nc4')
sst3.to_netcdf('/home/users/vg140344/datos/data/fogt/sst_ninio3_s4_dec_feb.nc4')
sst4.to_netcdf('/home/users/vg140344/datos/data/fogt/sst_ninio4_s4_dec_feb.nc4')

