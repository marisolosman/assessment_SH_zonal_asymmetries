import glob
import numpy as np
import os
import xarray as xr
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH = "/home/users/vg140344/datos/data/fogt/correlations/"

correlations = xr.open_mfdataset(PATH + "monthly_correlations_enso_SPoV_*.nc4")

intervals_ninio_WPV = np.empty([7, 2])
intervals_ninia_WPV = np.empty([7, 2])
intervals_ninio_SPV = np.empty([7, 2])
intervals_ninia_SPV = np.empty([7, 2])
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
for i in range(7):
	intervals_ninio_WPV[i, :] = np.percentile(correlations.correl_ninio_WPV[correlations.month==month[i]], [5, 95], axis=0)
	intervals_ninio_SPV[i, :] = np.percentile(correlations.correl_ninio_SPV[correlations.month==month[i]], [5, 95], axis=0)
	intervals_ninia_WPV[i, :] = np.percentile(correlations.correl_ninia_WPV[correlations.month==month[i]], [5, 95], axis=0)
	intervals_ninia_SPV[i, :] = np.percentile(correlations.correl_ninia_SPV[correlations.month==month[i]], [5, 95], axis=0)

ds = xr.Dataset({'intervals_ninio_WPV': (['month', 'quant'], intervals_ninio_WPV),
		 'intervals_ninia_WPV': (['month', 'quant'], intervals_ninia_WPV),
		 'intervals_ninio_SPV': (['month', 'quant'], intervals_ninio_SPV),
		 'intervals_ninia_SPV': (['month', 'quant'], intervals_ninia_SPV)},
		 coords={'seas': (['month'], month),
		 	 'quant': (['quant'], [5, 95])})
ds.to_netcdf('~/datos/data/fogt/monthly_statistics_correlations.nc4')
correlations = xr.open_mfdataset(PATH + "seasonal_correlations_enso_SPoV_*.nc4")

intervals_ninio_WPV = np.empty([5, 2])
intervals_ninia_WPV = np.empty([5, 2])
intervals_ninio_SPV = np.empty([5, 2])
intervals_ninia_SPV = np.empty([5, 2])

seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
for i in range(5):
	intervals_ninio_WPV[i, :] = np.percentile(correlations.correl_ninio_WPV[correlations.seas==seas[i]], [5, 95], axis=0)
	intervals_ninio_SPV[i, :] = np.percentile(correlations.correl_ninio_SPV[correlations.seas==seas[i]], [5, 95], axis=0)
	intervals_ninia_WPV[i, :] = np.percentile(correlations.correl_ninia_WPV[correlations.seas==seas[i]], [5, 95], axis=0)
	intervals_ninia_SPV[i, :] = np.percentile(correlations.correl_ninia_SPV[correlations.seas==seas[i]], [5, 95], axis=0)

ds = xr.Dataset({'intervals_ninio_WPV': (['seas', 'quant'], intervals_ninio_WPV),
		 'intervals_ninia_WPV': (['seas', 'quant'], intervals_ninia_WPV),
		 'intervals_ninio_SPV': (['seas', 'quant'], intervals_ninio_SPV),
		 'intervals_ninia_SPV': (['seas', 'quant'], intervals_ninia_SPV)},
		 coords={'seas': (['seas'], seas),
		 	 'quant': (['quant'], [5, 95])})
ds.to_netcdf('~/datos/data/fogt/seasonal_statistics_correlations.nc4')




	
