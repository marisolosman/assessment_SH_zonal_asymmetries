#sort el ninio events
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import eofdata
import matplotlib.pyplot as plt
ds = xr.open_dataset('./data/sst.nc')

ds.coords['year'] = np.arange(1981,2017)

#junto la cordenada year y number para categorizar enventos

ds = ds.stack(realiz = ['year', 'number'])
print(ds)
print(ds.realiz.values)
#compute seasonal means: ASO and SON
ninio34_aso = ds.sel(**{'month': slice(8,10)}).mean(dim='month')
ninio34_son = ds.sel(**{'month': slice(9,11)}).mean(dim='month')
ninio34_seasonal = xr.concat([ninio34_aso, ninio34_son], dim='season')
#select upper and lower quartile
ninio34_season_lower = ninio34_aso.quantile(0.25, dim='realiz', interpolation='linear')
ninio34_season_upper = ninio34_aso.quantile(0.75, dim='realiz', interpolation='linear')

ninio34_monthly_lower = ds.quantile(0.25, dim='realiz', interpolation='linear')
ninio34_monthly_upper = ds.quantile(0.75, dim='realiz', interpolation='linear')
[lamb, v, PC] = eofdata.eofdata(ds.sst.values, 4)
sst_monthly_index = -PC[0, :]

#ninio34_monthly_quartil = np.quantile(sst_index, [0.25, 0.75])

[lamb, v, PC] = eofdata.eofdata(ninio34_seasonal.sst.values, 8)
sst_seasonal_index = PC[0, :]

#ninio34_seasonal_quartil = np.quantile(sst_seasonal_index, [0.25, 0.75])

ds_new = xr.Dataset({'ninio34_mon': xr.DataArray(sst_monthly_index),
		     'ninio34_seas': xr.DataArray(sst_seasonal_index)})
ds_new.to_netcdf('./data/ninio34_index.nc')
ds.reset_index('realiz').to_netcdf('./data/ninio34_monthly.nc')
ninio34_seasonal.reset_index('realiz').to_netcdf('./data/ninio34_seasonal.nc')

ninio34_monthly_lower.to_netcdf('./data/ninio34_monthly_lower.nc')
ninio34_monthly_upper.to_netcdf('./data/ninio34_monthly_upper.nc')
ninio34_season_lower.to_netcdf('./data/ninio34_seasonal_lower.nc')
ninio34_season_upper.to_netcdf('./data/ninio34_seasonal_upper.nc')

#plot scatter PC against SST

plt.figure(1)
plt.scatter(sst_monthly_index, sst_seasonal_index, s=80)
plt.xlabel('PC monthly')
plt.ylabel('PC seasonal')
plt.savefig('./figures/scatter_PCs.jpg')
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
plt.scatter(sst_monthly_index, ninio34_seasonal.sst.values[0, :], s=80)
plt.xlabel('PC monthly')
plt.ylabel('SST ASO')
plt.subplot(426)
plt.scatter(sst_monthly_index, ninio34_seasonal.sst.values[1, :], s=80)
plt.xlabel('PC monthly')
plt.ylabel('SST SON')
plt.subplot(427)
plt.scatter(sst_seasonal_index, ninio34_seasonal.sst.values[0, :], s=80, marker=(5, 0))
plt.xlabel('PC seasonal')
plt.ylabel('SST ASO')
plt.subplot(428)
plt.scatter(sst_seasonal_index, ninio34_seasonal.sst.values[1, :], s=80, marker=(5, 0))
plt.xlabel('PC seasonal')
plt.ylabel('SST SON')
plt.savefig('./figures/scatter_PCs_SSTs.jpg')

#open era intreim data a compute el ninio index

sst_erai = xr.open_dataset('./data/sst_erai.nc')
sst_erai.time.values = sst_erai.valid_time.values

#compute ninio 3.4 index
ninio34_erai = sst_erai.sel(**{'time':slice('1981-08-01', '2018-02-01'), 'latitude':slice(5, -5),
			       'longitude':slice(190, 240)}).mean(dim=['longitude', 'latitude'])
#compute seasonal means: ASO and SON
ninio34_aso_erai = ninio34_erai.resample(time='QS-Aug').mean(dim='time', skipna='True')
ninio34_son_erai = ninio34_erai.resample(time='QS-Sep').mean(dim='time', skipna='True')

ninio34_erai_seasonal = xr.concat([ninio34_aso_erai.sel(time=(ninio34_aso_erai['time.month'] == 8)),
				   ninio34_son_erai.sel(time=(ninio34_son_erai['time.month'] == 9))],
				   dim='time').dropna(dim='time')

ninio34_erai_monthly = ninio34_erai.sel(time=np.logical_and(ninio34_erai['time.month'] >= 8,
					 ninio34_erai['time.month'] <= 11))
print(ninio34_erai_seasonal)
print(ninio34_erai_monthly)

[lamb, v, PC] = eofdata.eofdata(np.transpose(np.reshape(ninio34_erai_monthly.sst.values, [36, 4])), 4)
print(v[:, 0])
sst_monthly_index = -PC[0, :]

[lamb, v, PC] = eofdata.eofdata(np.reshape(ninio34_erai_seasonal.sst.values[~np.isnan(ninio34_erai_seasonal.sst.values)],
									   [2, 36]), 8)
sst_seasonal_index = PC[0, :]
print(v[:, 0])
ds_erai = xr.Dataset({'ninio34_mon': xr.DataArray(sst_monthly_index),
		      'ninio34_seas': xr.DataArray(sst_seasonal_index)})
ds_erai.to_netcdf('./data/ninio34_erai_index.nc')

