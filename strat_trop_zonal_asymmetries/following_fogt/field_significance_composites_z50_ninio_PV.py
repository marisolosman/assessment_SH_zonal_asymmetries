import glob
import xarray as xr
import numpy as np
import pandas as pd
from scipy import stats
RUTA = '/home/users/vg140344/datos/data/fogt/'
#lista = xr.open_mfdataset(RUTA + "correlations/seasonal_correlations_z50_enso_SPoV_*.nc4",
#			  parallel=True, combine='nested', chunks=10000, concat_dim='iter')
#file = RUTA + "seasonal_correlations_z50_enso_SPoV_polar.nc4"
#correlations = xr.open_dataset(file)
#seasons = lista.seas.values
#df = pd.DataFrame()
#for ii in lista.data_vars:
#	Percentiles = lista[ii].load().quantile([0.05, 0.95], dim='iter', interpolation='linear')
#	pvalue = [stats.percentileofscore(lista[ii][:, lista.seas==i].values,
#				       correlations[ii].sel(**{'seas':i}).values) for i in seasons]
#	df = df.append(pd.DataFrame({'Field': ii, 'Percentile':pvalue,
#					     '5th Percentile':Percentiles[0, :],
#					     '95th Percentile': Percentiles[1, :]},
#					     index=seasons, columns=['Field', 'Percentile', 
#					     '5th Percentile', '95th Percentile']))
#
#df.to_csv(RUTA + "percentile_seasonal_correlations_composites_z50_enso_SPoV.csv")
#
lista = xr.open_mfdataset(RUTA + "correlations/monthly_correlations_z50_enso_SPoV_*.nc4",
			  parallel=True, combine='nested', chunks=10000, concat_dim='iter')
file = RUTA + "monthly_correlations_z50_enso_SPoV_polar.nc4"
correlations = xr.open_dataset(file)

months = lista.month.values
df = pd.DataFrame()
for ii in lista.data_vars:
	Percentiles = lista[ii].load().quantile([0.05, 0.95], dim='iter', interpolation='linear')
	pvalue = [stats.percentileofscore(lista[ii][:, lista.month==i].values,
				       correlations[ii].sel(**{'month':i}).values) for i in months]
	df = df.append(pd.DataFrame({'Field': ii, 'Percentile':pvalue,
					     '5th Percentile':Percentiles[0, :],
					     '95th Percentile': Percentiles[1, :]},
					     index=months, columns=['Field', 'Percentile', 
					     '5th Percentile', '95th Percentile']))
df.to_csv(RUTA + "percentile_monthly_correlations_composites_z50_enso_SPoV.csv")


