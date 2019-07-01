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
	if Y != 2002:
		ds = xr.open_dataset(RUTA + 'HGT50_erai_' + str(Y) + '.grib',
				     engine='cfgrib', backend_kwargs={'indexpath':''})
		ds = ds.drop(['number','isobaricInhPa', 'step'])
		hgt.append(ds)

hgt = xr.concat(hgt, dim='time')
hgt.time.values = np.arange(0, len(hgt.time.values))
hgt.to_netcdf('~/datos/data/hgt_erai_50.nc4')
hgt = hgt.sel(**{'latitude':slice(-60, -90)}).mean(dim=['longitude','latitude'])
hgt.to_netcdf('~/datos/data/PV_monthly_erai.nc4')

hgt.time.values = hgt.valid_time.values
hgt = hgt.sel(**{'time':slice('1981-01-01', '2017-12-31')})


PV_ason = hgt.sel(time=np.logical_and(hgt['time.month']>=8, hgt['time.month']<=11))
PV = np.reshape(PV_ason.z.values,[36, 4])
PV = (PV - np.mean(PV, axis=0))
[lamb, v, PC] = eofdata.eofdata(PV.T, 3)
print(v[:, 0])
PV_monthly_index = -1 * PC[0, :]
ds_new = xr.Dataset({'PV_mon': xr.DataArray(PV_monthly_index)})
ds_new.to_netcdf('~/datos/data/PV_index_erai.nc4')

