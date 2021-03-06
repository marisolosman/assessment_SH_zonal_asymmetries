#sort el ninio events
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import eofdata
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

RUTA = '~/datos/ERAI/'
hgt = []

for Y in np.arange(1981,2019):
	ds = xr.open_dataset(RUTA + 'HGT50_erai_' + str(Y) + '.grib',
			     engine='cfgrib', backend_kwargs={'indexpath':''})
	ds = ds.drop(['number','isobaricInhPa', 'step'])
	hgt.append(ds)

hgt = xr.concat(hgt, dim='time')
hgt.time.values = hgt.valid_time.values
hgt = hgt.sel(time= np.logical_or(hgt.time.values <=np.datetime64('2002-07-31'), hgt.time.values>=np.datetime64('2003-08-01')))
hgt.time.values = np.arange(0, len(hgt.time.values))
hgt.to_netcdf('~/datos/data/fogt/hgt_erai_50.nc4')
hgt.time.values = hgt.valid_time.values

hgt = hgt.sel(**{'latitude':slice(-60, -90)}).mean(dim=['longitude','latitude'])
hgt.time.values = np.arange(0, len(hgt.time.values))
hgt.to_netcdf('~/datos/data/PV_monthly_erai.nc4')
hgt.time.values = hgt.valid_time.values
hgt = hgt.sel(**{'time':slice('1981-08-01', '2018-02-28')})

SPV = hgt.sel(time=np.logical_or(hgt['time.month']>=8, hgt['time.month']<=2))

SPV = np.reshape(SPV.z.values,[36, 7])
SPV = (SPV - np.mean(SPV, axis=0))
[lamb, v, PC] = eofdata.eofdata(SPV.T, 3)
print(v[:, 0])
SPV_monthly_index = PC[0, :]
ds_new = xr.Dataset({'SPV_mon': xr.DataArray(SPV_monthly_index)})
ds_new.to_netcdf('~/datos/data/fogt/SPV_index_erai.nc4')

