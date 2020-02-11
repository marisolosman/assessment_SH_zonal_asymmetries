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

ds = xr.open_dataset(RUTA + 'sst_ninio34_aug_feb.nc4')

ninio34_monthly_lower = ds.quantile(0.25, dim='realiz', interpolation='linear')
ninio34_monthly_upper = ds.quantile(0.75, dim='realiz', interpolation='linear')
[lamb, v, PC] = eofdata.eofdata(ds.sst.values, 4)
print(v[:, 0])
ninio34_monthly_index = -PC[0, :]

ds_new = xr.Dataset({'ninio34_index': xr.DataArray(ninio34_monthly_index)})
ds_new.to_netcdf('~/datos/data/fogt/ninio34_monthly.nc4')

#plot scatter PC against SST
#
#plt.figure(1)
#plt.scatter(sst_monthly_index, sst_seasonal_index, s=80)
#plt.xlabel('PC monthly')
#plt.ylabel('PC seasonal')
#plt.savefig('./figures/scatter_PCs.jpg')
#plt.figure(2, (12, 16))
#plt.subplot(421)
#plt.scatter(sst_monthly_index, ds.sst.values[0, :], s=80)
#plt.xlabel('PC monthly')
#plt.ylabel('SST Aug')
#plt.subplot(422)
#plt.scatter(sst_monthly_index, ds.sst.values[1, :], s=80)
#plt.xlabel('PC monthly')
#plt.ylabel('SST Sep')
#plt.subplot(423)
#plt.scatter(sst_monthly_index, ds.sst.values[2, :], s=80)
#plt.xlabel('PC monthly')
#plt.ylabel('SST Oct')
#plt.subplot(424)
#plt.scatter(sst_monthly_index, ds.sst.values[3, :], s=80)
#plt.xlabel('PC monthly')
#plt.ylabel('SST Nov')
#plt.subplot(425)
#plt.scatter(sst_monthly_index, ninio34_seasonal.sst.values[0, :], s=80)
#plt.xlabel('PC monthly')
#plt.ylabel('SST ASO')
#plt.subplot(426)
#plt.scatter(sst_monthly_index, ninio34_seasonal.sst.values[1, :], s=80)
#plt.xlabel('PC monthly')
#plt.ylabel('SST SON')
#plt.subplot(427)
#plt.scatter(sst_seasonal_index, ninio34_seasonal.sst.values[0, :], s=80, marker=(5, 0))
#plt.xlabel('PC seasonal')
#plt.ylabel('SST ASO')
#plt.subplot(428)
#plt.scatter(sst_seasonal_index, ninio34_seasonal.sst.values[1, :], s=80, marker=(5, 0))
#plt.xlabel('PC seasonal')
#plt.ylabel('SST SON')
#plt.savefig('./figures/scatter_PCs_SSTs.jpg')
#==================================================================================================
##open era intreim data a compute el ninio index
#
#sst_erai = xr.open_dataset('~/datos/data/fogt/sst_erai.nc4')
#sst_erai.time.values = sst_erai.valid_time.values
#
##remove 2002-2003 years 
#sst_erai = sst_erai.sel(time= np.logical_or(sst_erai.time.values <=np.datetime64('2002-07-31'), sst_erai.time.values>=np.datetime64('2003-08-01')))
#
##compute ninio 3.4 index
#ninio34_erai = sst_erai.sel(**{'time':slice('1981-08-01', '2018-02-01'), 'latitude':slice(5, -5), 'longitude':slice(190, 240)}).mean(
#			dim=['longitude', 'latitude'])
#
#ninio34_erai = ninio34_erai.sel(time = np.logical_or(ninio34_erai['time.month']>=8, ninio34_erai['time.month']<=2))
#
#ninio34_erai = np.reshape(ninio34_erai.sst.values,[36, 7])
#
#[lamb, v, PC] = eofdata.eofdata(np.transpose(ninio34_erai), 4)
#print(v[:, 0])
#sst_monthly_index = -PC[0, :]
#
#ds = xr.Dataset({'ninio34_index': xr.DataArray(sst_monthly_index)})
#ds.to_netcdf('~/datos/data/fogt/ninio34_erai_index.nc4')
#
