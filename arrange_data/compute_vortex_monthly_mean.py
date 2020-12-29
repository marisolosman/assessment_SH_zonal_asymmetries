#compute monthly PoV
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

RUTA = '/storage/shared/glusterfs/acrcc/users/ys916780/Data_for_Sol/output'
hgt = []
for Y in np.arange(1981,2018):
	if Y !=2002:
		ds = xr.open_dataset(RUTA + '_'+ str(Y) + '_geopotential50', engine='cfgrib', chunks={'step': 10},
				     backend_kwargs={'indexpath': ''})
		ds = ds.rename({'time': 'IC'}).set_coords(['IC'])
		ds['step'] = pd.DatetimeIndex(ds.valid_time.values)
		ds = ds.rename({'step': 'time'}).set_coords(['time'])
		ds = ds.groupby('time.month').mean(dim='time')
		ds = ds.sel(**{'latitude': slice(-60, -90)}).mean(dim=['longitude', 'latitude']).compute()
		hgt.append(ds)

hgt =xr.concat(hgt, dim='year')
hgt.to_netcdf('./data/PV_monthly.nc')


