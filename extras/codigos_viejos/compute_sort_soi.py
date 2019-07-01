#compute monthly and seasonal soi and sort events
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd

ds = xr.open_dataset('./data/msl.nc')

#SOI MSLP differences between tahiti 17.40S 149.28W and darwin 12.26S 130.5 E
ds = ds.sel(longitude=(360 - 149.28), latitude=-17.4, method='nearest') -\
	ds.sel(longitude=130.5, latitude=-12.26, method='nearest')

ds.coords['year'] = np.arange(1981,2017)

#junto la cordenada year y number para categorizar enventos

ds = ds.stack(realiz = ['year', 'number'])
print(ds)
#compute seasona means: ASO and SON
soi_aso = ds.sel(**{'month':slice(8,10)}).mean(dim='month')
soi_son = ds.sel(**{'month':slice(9,11)}).mean(dim='month')

#select upper and lower quartile
soi_aso_lower = soi_aso.quantile(0.25, dim='realiz', interpolation='linear')
soi_aso_upper = soi_aso.quantile(0.75, dim='realiz', interpolation='linear')
soi_son_lower = soi_aso.quantile(0.25, dim='realiz', interpolation='linear')
soi_son_upper = soi_aso.quantile(0.75, dim='realiz', interpolation='linear')

soi_monthly_lower = ds.quantile(0.25, dim='realiz', interpolation='linear')
soi_monthly_upper = ds.quantile(0.75, dim='realiz', interpolation='linear')

ds.reset_index('realiz').to_netcdf('./data/soi_monthly.nc')
soi_aso.reset_index('realiz').to_netcdf('./data/soi_aso.nc')
soi_son.reset_index('realiz').to_netcdf('./data/soi_son.nc')

soi_aso_lower.to_netcdf('./data/soi_aso_lower.nc')
soi_aso_upper.to_netcdf('./data/soi_aso_upper.nc') 
soi_son_lower.to_netcdf('./data/soi_son_lower.nc') 
soi_son_upper.to_netcdf('./data/soi_son_upper.nc')
soi_monthly_lower.to_netcdf('./data/soi_monthly_lower.nc')
soi_monthly_upper.to_netcdf('./data/soi_monthly_upper.nc')






