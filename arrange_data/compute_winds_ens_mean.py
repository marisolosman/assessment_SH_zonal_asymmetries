# compute climatology monthly winds
import numpy as np
import xarray as xr
import pandas as pd
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

ds = xr.open_mfdataset('/home/users/vg140344/datos/UV_S4Hindcasts_monthly_Mlevels_*.grib', parallel=False,
			 combine='nested', concat_dim=['time'], engine='cfgrib')
ds = ds.mean('number')

ds.to_netcdf('~/datos/data/UV_S4_Mlevels_mean.nc4')
ds = xr.open_dataset('~/datos/data/UV_S4_Mlevels_mean.nc4')
ds = ds.where(ds.time != np.datetime64("2002-08-01"), drop=True)
ds1 = ds.stack(realiz=['time', 'step'])

times = np.ndarray.flatten(ds.valid_time.values)

ds1.coords['realiz'] = times

ds1 = ds1.groupby('realiz.month').mean(skipna=True)

ds1 = ds1.sel(**{'latitude': slice(10, -90)}).loc[{'month': [9, 10, 11, 12, 1, 2, 3]}]

ds1.coords['month'] = [8, 9, 10, 11, 12, 1, 2]

ds1.to_netcdf('~/datos/data/UV_S4_Mlevels_clim.nc4')
