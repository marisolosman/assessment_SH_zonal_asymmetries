import glob
import xarray as xr
import numpy as np
import pandas as pd
from scipy import stats

RUTA = '/home/users/vg140344/datos/data/fogt/'
lista = xr.open_mfdataset(RUTA + "correlations/monthly_correlations_z50_new_*.nc4",
                          parallel=True, combine='nested', chunks=10000, concat_dim="iter", coords="different", compat="broadcast_equals")
lista['iter'] = np.arange(lista.iter.values.shape[0])
lista.to_netcdf(RUTA + 'monthly_correlations_z50_new.nc4')

file = RUTA + "monthly_correlations_z50_new.nc4"
samples = xr.open_dataset(file)
file2 = RUTA + 'monthly_z50_new_correlations.nc4'
correlations = xr.open_dataset(file2)
months = samples.month.values
df = pd.DataFrame()
for ii in samples.data_vars:
    Percentiles = samples[ii].load().quantile([0.05, 0.95], dim='iter', interpolation='linear')
    corr_coef = [correlations[ii].sel(**{'month': i}).values for i in months]
    pvalue = [stats.percentileofscore(samples[ii][:, samples.month == i].values,
                                      -1 * correlations[ii].sel(**{'month': i}).values) for i in months]
    df = df.append(pd.DataFrame({'Field': ii, 'Correlation': corr_coef, 'Percentile': pvalue,
                                 '5th Percentile': Percentiles[0, :],
                                 '95th Percentile': Percentiles[1, :]},
                                index=months, columns=['Field', 'Correlation',
                                                       'Percentile',
                                                       '5th Percentile', '95th Percentile']))
df.to_csv(RUTA + "percentile_monthly_correlations_z50_new.csv")


RUTA = '/home/users/vg140344/datos/data/fogt/'
lista = xr.open_mfdataset(RUTA + "correlations/monthly_correlations_z200_new_*.nc4",
                          parallel=True, combine='nested', chunks=10000, concat_dim="iter", coords="different", compat="broadcast_equals")
lista['iter'] = np.arange(lista.iter.values.shape[0])
lista.to_netcdf(RUTA + 'monthly_correlations_z200_new.nc4')

file = RUTA + "monthly_correlations_z200_new.nc4"
samples = xr.open_dataset(file)
file2 = RUTA + 'monthly_z200_new_correlations.nc4'
correlations = xr.open_dataset(file2)
months = samples.month.values
df = pd.DataFrame()
for ii in samples.data_vars:
    Percentiles = samples[ii].load().quantile([0.05, 0.95], dim='iter', interpolation='linear')
    corr_coef = [correlations[ii].sel(**{'month': i}).values for i in months]
    pvalue = [stats.percentileofscore(samples[ii][:, samples.month == i].values,
                                      -1 * correlations[ii].sel(**{'month': i}).values) for i in months]
    df = df.append(pd.DataFrame({'Field': ii, 'Correlation': corr_coef, 'Percentile': pvalue,
                                 '5th Percentile': Percentiles[0, :],
                                 '95th Percentile': Percentiles[1, :]},
                                index=months, columns=['Field', 'Correlation',
                                                       'Percentile',
                                                       '5th Percentile', '95th Percentile']))
df.to_csv(RUTA + "percentile_monthly_correlations_z200_new.csv")
