#compute PoV index for ERAI data
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import eofdata
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

RUTA = '~/datos/data/'
FILE = 'hgt_erai_50.nc4'
hgt = xr.open_dataset(RUTA + FILE)
hgt.time.values = hgt.valid_time.values
hgt = hgt.sel(**{'latitude': slice(-60, -90)}).mean(dim=['longitude', 'latitude'])
period = np.logical_and(hgt.time >= np.datetime64('1981-08-01'), hgt.time < np.datetime64('2018-03-01'))
hgt = hgt.sel(time=period)
SPV = hgt.sel(time=np.logical_or(hgt['time.month'] >= 8, hgt['time.month'] <= 2))

SPV = np.reshape(SPV.z.values,[37, 7])
SPV = (SPV - np.mean(SPV, axis=0))
[lamb, v, PC] = eofdata.eofdata(SPV.T, 3)
print(v[:, 0])
SPV_monthly_index = PC[0, :]

ds_new = xr.Dataset({'SPV_mon': xr.DataArray(SPV_monthly_index)})
ds_new.to_netcdf(RUTA + 'fogt/SPV_index_erai_inc_2002.nc4')

