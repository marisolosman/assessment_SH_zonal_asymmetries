#define ninio and ninia years for S4
#compute ninio indexes for erai and define ninio and ninia years
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import eofdata
import matplotlib.pyplot as plt
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
RUTA = '~/datos/data/fogt/'

ds = xr.open_dataset(RUTA + 'sst_ninio34_s4_aug_feb_inc_2002.nc4')

ninio34_monthly_lower = ds.quantile(0.25, dim='realiz', interpolation='linear')
ninio34_monthly_upper = ds.quantile(0.75, dim='realiz', interpolation='linear')
[lamb, v, PC] = eofdata.eofdata(ds.sst.values, 4)
print(v[:, 0])
ninio34_monthly_index = -PC[0, :]

ds_new = xr.Dataset({'ninio34_index': xr.DataArray(ninio34_monthly_index)})
ds_new.to_netcdf('~/datos/data/fogt/ninio34_index_s4_inc_2002.nc4')

#==================================================================================================
##open era intreim data a compute el ninio index

sst_erai = xr.open_dataset('~/datos/data/fogt/sst_erai.nc4', cache=False)
sst_erai.time.values = sst_erai.valid_time.values

period = np.logical_and(sst_erai.time.values>=np.datetime64('1981-08-01'), sst_erai.time.values<=np.datetime64('2018-02-01'))

#compute ninio 3.4 index
ninio34_erai = sst_erai.sel(time=period, latitude=slice(5, -5),
			    longitude=slice(190, 240)).mean(dim=['longitude', 'latitude'])

ninio34_erai = ninio34_erai.sel(time=np.logical_or(ninio34_erai['time.month'] >= 8, ninio34_erai['time.month'] <= 2))
ninio34_erai = np.reshape(ninio34_erai.sst.values, [37, 7])
[lamb, v, PC] = eofdata.eofdata(np.transpose(ninio34_erai), 4)
print(v[:, 0])
sst_monthly_index = -PC[0, :]
ds = xr.Dataset({'ninio34_index': xr.DataArray(sst_monthly_index)})
ds.to_netcdf('~/datos/data/fogt/ninio34_erai_index_inc_2002.nc4')

