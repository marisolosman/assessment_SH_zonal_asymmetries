#sort el ninio events
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import eofdata

ds = xr.open_dataset('./data/PV_monthly.nc')

ds.coords['year'] = np.arange(1981,2017)

#junto la cordenada year y number para categorizar enventos

ds = ds.stack(realiz = ['year', 'number']).compute()
print(ds)
#compute seasona means: ASO and SON
PV_aso = ds.sel(**{'month':slice(8,10)}).mean(dim='month')
PV_son = ds.sel(**{'month':slice(9,11)}).mean(dim='month')
PV_seasonal = xr.concat([PV_aso, PV_son], dim='season')
#select upper and lower quartile
PV_aso_lower = PV_aso.quantile(0.25, dim='realiz', interpolation='linear')
PV_aso_upper = PV_aso.quantile(0.75, dim='realiz', interpolation='linear')
PV_son_lower = PV_aso.quantile(0.25, dim='realiz', interpolation='linear')
PV_son_upper = PV_aso.quantile(0.75, dim='realiz', interpolation='linear')

PV_monthly_lower = ds.quantile(0.25, dim='realiz', interpolation='linear')
PV_monthly_upper = ds.quantile(0.75, dim='realiz', interpolation='linear')

ds.reset_index('realiz').to_netcdf('./data/PV_monthly2.nc')

[lamb, v, PC] = eofdata.eofdata(ds.z.values[0:4, :], 3)

print(v[:, 0])
PV_monthly_index = PC[0, :]

[lamb, v, PC] = eofdata.eofdata(PV_seasonal.z.values, 2)
print(v[0,:])
PV_seasonal_index = -1 * PC[0, :]

ds_new = xr.Dataset({'PV_mon': xr.DataArray(PV_monthly_index), 'PV_seas': xr.DataArray(PV_seasonal_index)})
ds_new.to_netcdf('./data/PV_index.nc')

PV_aso.reset_index('realiz').to_netcdf('./data/PV_aso.nc')
PV_son.reset_index('realiz').to_netcdf('./data/PV_son.nc')

PV_aso_lower.to_netcdf('./data/PV_aso_lower.nc')
PV_aso_upper.to_netcdf('./data/PV_aso_upper.nc') 
PV_son_lower.to_netcdf('./data/PV_son_lower.nc') 
PV_son_upper.to_netcdf('./data/PV_son_upper.nc')
PV_monthly_lower.to_netcdf('./data/PV_monthly_lower.nc')
PV_monthly_upper.to_netcdf('./data/PV_monthly_upper.nc')






