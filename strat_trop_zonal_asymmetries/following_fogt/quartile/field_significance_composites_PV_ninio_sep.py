import glob
import xarray as xr
import numpy as np
import pandas as pd
from scipy import stats
RUTA = '/home/users/vg140344/datos/data/fogt/'
lista = xr.open_mfdataset(RUTA + "correlations/seasonal*_SPoV_enso*_sep_q.nc4", chunks=None)
lista.to_netcdf(RUTA + 'seasonal_correlations_z200_SPoV_enso_polar_10k_sep_q.nc4')

file = RUTA + "seasonal_correlations_SPoV_enso_polar_sep_q.nc4"
correlations = xr.open_dataset(file)
seasons = np.unique(lista.seas.values)
df = pd.DataFrame() 
for i in seasons:
	for ii in lista.data_vars:
		aux = np.percentile(lista[ii].values[lista[ii].seas == i], [5, 95])
		pvalue = stats.percentileofscore(lista[ii].values[lista.seas == i],
					      correlations[ii].sel(**{'seas':i}).values)
		df = df.append(pd.DataFrame({'Field': ii,
					     'Correlation':[correlations[ii].sel(**{'seas':i}).values],
					     'Percentile':[pvalue], '5th Percentile': [aux[0]],
					     '95th Percentile': [aux[1]]},
					     index=[i], columns=['Field', 'Correlation',
								 'Percentile', '5th Percentile',
								 '95th Percentile']))
df = df.loc[['ASO', 'SON', 'OND', 'NDJ', 'DJF']]

df.to_csv(RUTA + "percentile_seasonal_correlations_composites_SPoV_enso_sep_q.csv")

lista = xr.open_mfdataset(RUTA + "correlations/monthly*_SPoV_enso*_sep_q.nc4", chunks=70000)
lista.to_netcdf(RUTA + 'monthly_correlations_z200_SPoV_enso_polar_10k_sep_q.nc4')
file = RUTA + "monthly_correlations_SPoV_enso_polar_sep_q.nc4"
correlations = xr.open_dataset(file)

months = np.unique(lista.month.values)
df = pd.DataFrame()
for i in months:
	for ii in lista.data_vars:
		aux = np.percentile(lista[ii].values[lista[ii].month == i], [5, 95])
		pvalue = stats.percentileofscore(lista[ii].values[lista.month == i],
					      correlations[ii].sel(**{'month':i}).values)
		df = df.append(pd.DataFrame({'Field': ii, 'Correlation':[correlations[ii].sel(**{'month':i}).values], 'Percentile':[pvalue],
					     'Percentile':[pvalue], '5th Percentile': [aux[0]],
					     '95th Percentile': [aux[1]]},
					     index=[i], columns=['Field', 'Correlation',
								 'Percentile', '5th Percentile',
								 '95th Percentile']))

df = df.loc[['Aug', 'Sep', 'Oct', 'Nov', 'Ded', 'Jan', 'Feb']]

df.to_csv(RUTA + "percentile_monthly_correlations_composites_SPoV_enso_sep_q.csv")


