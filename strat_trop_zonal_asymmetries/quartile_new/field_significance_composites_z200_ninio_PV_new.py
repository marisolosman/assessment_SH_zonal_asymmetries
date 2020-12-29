import glob
import xarray as xr
import numpy as np
import pandas as pd
from scipy import stats
import datetime
start = datetime.datetime.now()
RUTA = '/home/users/vg140344/datos/data/fogt/'
lista = xr.open_mfdataset(RUTA + "correlations/seasonal_correlations_z200_enso_SPoV_*_new.nc4",
			  parallel=True, combine='nested', chunks=10000, concat_dim='iter')

lista.to_netcdf(RUTA + 'seasonal_correlations_z200_enso_SPoV_10k_new.nc4')
lista = xr.open_dataset(RUTA + 'seasonal_correlations_z200_enso_SPoV_10k_new.nc4')
file = RUTA + "seasonal_correlations_enso_SPoV_polar_new.nc4"
correlations = xr.open_dataset(file)
seasons = lista.seas.values
df = pd.DataFrame()
for ii in lista.data_vars:
	Percentiles = lista[ii].load().quantile([0.05, 0.95], dim='iter', interpolation='linear')
	corr_coef = [correlations[ii].sel(**{'seas':i}).values for i in seasons]

	pvalue = [stats.percentileofscore(lista[ii][:, lista.seas==i].values,
				       correlations[ii].sel(**{'seas':i}).values) for i in seasons]
	df = df.append(pd.DataFrame({'Field': ii,'Correlation':corr_coef, 'Percentile':pvalue,
					     '5th Percentile':Percentiles[0, :],
					     '95th Percentile': Percentiles[1, :]},
					     index=seasons, columns=['Field', 'Correlation',
								     'Percentile', 
					     '5th Percentile', '95th Percentile']))
df.to_csv(RUTA + "percentile_seasonal_correlations_composites_z200_enso_SPoV_new.csv")


lista = xr.open_mfdataset(RUTA + "correlations/monthly_correlations_z200_enso_SPoV_*_new.nc4",
			  parallel=True, combine='nested', chunks=10000, concat_dim='iter')

lista.to_netcdf(RUTA + 'monthly_correlations_z200_enso_SPoV_10k_new.nc4')

file = RUTA + "monthly_correlations_enso_SPoV_polar_new.nc4"
correlations = xr.open_dataset(file)
months = lista.month.values
df = pd.DataFrame()
for ii in lista.data_vars:
	Percentiles = lista[ii].load().quantile([0.05, 0.95], dim='iter', interpolation='linear')
	corr_coef = [correlations[ii].sel(**{'month':i}).values for i in months]

	pvalue = [stats.percentileofscore(lista[ii][:, lista.month==i].values,
				       correlations[ii].sel(**{'month':i}).values) for i in months]
	df = df.append(pd.DataFrame({'Field': ii, 'Correlation':corr_coef, 'Percentile':pvalue,
					     '5th Percentile':Percentiles[0, :],
					     '95th Percentile': Percentiles[1, :]},
					     index=months, columns=['Field', 'Correlation',
							            'Percentile', 
					     '5th Percentile', '95th Percentile']))
df.to_csv(RUTA + "percentile_monthly_correlations_composites_z200_enso_SPoV_new.csv")

print(datetime-datetime.now()-start)
