import glob
import xarray as xr
import numpy as np
import pandas as pd
from scipy import stats
RUTA = '/home/users/vg140344/datos/data/fogt/'
lista = xr.open_mfdataset(RUTA + "correlations/seasonal_correlations_z50_SPoV_enso_*.nc4",
			  parallel=True, combine='nested', chunks=10000, concat_dim='iter')

lista.to_netcdf(RUTA + 'seasonal_correlations_z50_SPoV_enso_polar_10k.nc4')

lista = xr.open_mfdataset(RUTA + "correlations/monthly_correlations_z50_SPoV_enso_*.nc4",
			  parallel=True, combine='nested', chunks=10000, concat_dim='iter')
lista.to_netcdf(RUTA + 'monthly_correlations_z50_SPoV_enso_polar_10k.nc4')

