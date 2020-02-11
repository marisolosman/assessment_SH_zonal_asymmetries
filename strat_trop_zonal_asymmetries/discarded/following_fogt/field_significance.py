import glob
import xarray as xr
import numpy as np
from scipy import stats
import pandas as pd
RUTA = '/home/users/vg140344/datos/data/fogt/'

#lista = xr.open_mfdataset(RUTA + "correlations/seasonal*_enso_SPoV_*.nc4")

#seasons = np.unique(lista.seas.values)
#df = pd.DataFrame()

file = RUTA +  "seasonal_correlations_enso_SPoV_polar.nc4"
correlations = xr.open_dataset(file)

#for i in seasons:
#	for ii in lista.data_vars:
#		aux = stats.percentileofscore(lista[ii].values[lista.seas == i],
#					      correlations[ii].sel(**{'seas':i}).values)
#		df = df.append(pd.DataFrame({'Season': [i], 'Field': [ii], 'Percentile': [aux]},
#					    columns=['Season', 'Field', 'Percentile']))
#
#df.to_csv(RUTA + "percentile_seasonal_correlations_composites_enso_SPoV.csv")

correlations.to_dataframe().to_csv(RUTA + "seasonal_corelations_composites_enso_SPoV.csv")

lista = xr.open_mfdataset(RUTA + "correlations/monthly*_enso_SPoV_*.nc4")
months = np.unique(lista.month.values)

file = RUTA + "monthly_correlations_enso_SPoV_polar.nc4"
correlations = xr.open_dataset(file)

df = pd.DataFrame()
for i in months:
	for ii in lista.data_vars:
		aux = stats.percentileofscore(lista[ii].values[lista.month == i],
					      correlations[ii].sel(**{'month':i}).values)
		df = df.append(pd.DataFrame({'Month': [i], 'Field': [ii], 'Percentile': [aux]},
					     columns=['Month', 'Field', 'Percentile']))

df.to_csv(RUTA + "percentile_monthly_correlations_composites_enso_SPoV.csv")

correlations.to_dataframe().to_csv(RUTA + "monthly_corelations_composites_enso_SPoV.csv")

