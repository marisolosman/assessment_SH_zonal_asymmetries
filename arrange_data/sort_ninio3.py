# compute ninio 3 index ad compare with monthly and seasonal values
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import eofdata
import matplotlib.pyplot as plt
import os

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

ds = xr.open_dataset('~/datos/data/sst_ninio3.nc')

ds.coords['year'] = np.arange(1981, 2017)

# junto la cordenada year y number para categorizar enventos

ds = ds.stack(realiz = ['year', 'number'])
#compute seasonal means: ASO and SON
ninio3_aso = ds.sel(**{'month': slice(8, 10)}).mean(dim='month')
ninio3_son = ds.sel(**{'month': slice(9, 11)}).mean(dim='month')
ninio3_seasonal = xr.concat([ninio3_aso, ninio3_son], dim='season')
#select upper and lower quartile
ninio3_season_lower = ninio3_aso.quantile(0.25, dim='realiz', interpolation='linear')
ninio3_season_upper = ninio3_aso.quantile(0.75, dim='realiz', interpolation='linear')

ninio3_monthly_lower = ds.quantile(0.25, dim='realiz', interpolation='linear')
ninio3_monthly_upper = ds.quantile(0.75, dim='realiz', interpolation='linear')
[lamb, v, PC] = eofdata.eofdata(ds.sst.values, 4)
print(v[:,0])
sst_monthly_index = -PC[0, :]


[lamb, v, PC] = eofdata.eofdata(ninio3_seasonal.sst.values, 8)
sst_seasonal_index = PC[0, :]


ds_new = xr.Dataset({'ninio3_mon': xr.DataArray(sst_monthly_index), 'ninio3_seas': xr.DataArray(sst_seasonal_index)})
ds_new.to_netcdf('~/datos/data/ninio3_index.nc')
ds.reset_index('realiz').to_netcdf('~/datos/data/ninio3_monthly.nc')
ninio3_seasonal.reset_index('realiz').to_netcdf('~/datos/data/ninio3_seasonal.nc')

#plot scatter PC against SST

plt.figure(1)
plt.scatter(sst_monthly_index, sst_seasonal_index, s=80)
plt.xlabel('PC monthly')
plt.ylabel('PC seasonal')
plt.savefig('/home/users/vg140344/figures/scatter_PCs_ninio3.jpg')
plt.figure(2, (12, 16))
plt.subplot(421)
plt.scatter(sst_monthly_index, ds.sst.values[0, :], s=80)
plt.xlabel('PC monthly')
plt.ylabel('SST Aug')
plt.subplot(422)
plt.scatter(sst_monthly_index, ds.sst.values[1, :], s=80)
plt.xlabel('PC monthly')
plt.ylabel('SST Sep')
plt.subplot(423)
plt.scatter(sst_monthly_index, ds.sst.values[2, :], s=80)
plt.xlabel('PC monthly')
plt.ylabel('SST Oct')
plt.subplot(424)
plt.scatter(sst_monthly_index, ds.sst.values[3, :], s=80)
plt.xlabel('PC monthly')
plt.ylabel('SST Nov')
plt.subplot(425)
plt.scatter(sst_monthly_index, ninio3_seasonal.sst.values[0, :], s=80)
plt.xlabel('PC monthly')
plt.ylabel('SST ASO')
plt.subplot(426)
plt.scatter(sst_monthly_index, ninio3_seasonal.sst.values[1, :], s=80)
plt.xlabel('PC monthly')
plt.ylabel('SST SON')
plt.subplot(427)
plt.scatter(sst_seasonal_index, ninio3_seasonal.sst.values[0, :], s=80, marker=(5, 0))
plt.xlabel('PC seasonal')
plt.ylabel('SST ASO')
plt.subplot(428)
plt.scatter(sst_seasonal_index, ninio3_seasonal.sst.values[1, :], s=80, marker=(5, 0))
plt.xlabel('PC seasonal')
plt.ylabel('SST SON')
plt.savefig('/home/users/vg140344/figures/scatter_PCs_SSTs_ninio3.jpg')

#open era intreim data a compute el ninio index

sst_erai = xr.open_dataset('~/datos/data/sst_erai.nc')
sst_erai.time.values = sst_erai.valid_time.values

#compute ninio 4 index
ninio3_erai = sst_erai.sel(**{'time': slice('1981-08-01', '2018-02-01'), 'latitude':slice(5, -5),
			      'longitude': slice(210, 270)}).mean(dim=['longitude', 'latitude'])
#compute seasonal means: ASO and SON
ninio3_aso_erai = ninio3_erai.resample(time='QS-Aug').mean(dim='time', skipna='True')
ninio3_son_erai = ninio3_erai.resample(time='QS-Sep').mean(dim='time', skipna='True')

ninio3_erai_seasonal = xr.concat([ninio3_aso_erai.sel(time=(ninio3_aso_erai['time.month'] == 8)),
				 ninio3_son_erai.sel(time=(ninio3_son_erai['time.month'] == 9))],
				 dim='time').dropna(dim='time')

ninio3_erai_monthly = ninio3_erai.sel(time = np.logical_and(ninio3_erai['time.month'] >= 8,
				      ninio3_erai['time.month'] <= 11))

[lamb, v, PC] = eofdata.eofdata(np.transpose(np.reshape(ninio3_erai_monthly.sst.values, [36, 4])), 4)
print(v[:, 0])
sst_monthly_index = -PC[0, :]

[lamb, v, PC] = eofdata.eofdata(np.reshape(ninio3_erai_seasonal.sst.values[~np.isnan(ninio3_erai_seasonal.sst.values)],
                                                                           [2, 36]), 8)
sst_seasonal_index = PC[0, :]
print(v[:, 0])
ds_erai = xr.Dataset({'ninio3_mon': xr.DataArray(sst_monthly_index), 'ninio3_seas': xr.DataArray(sst_seasonal_index)})
ds_erai.to_netcdf('~/datos/data/ninio3_erai_index.nc')

