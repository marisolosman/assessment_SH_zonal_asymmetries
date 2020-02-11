#composites of EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import plots
import sys
from datetime import datetime

def TestCorrelation(var1, var2, var3, var4):
	aux1 = np.concatenate((var1, var3), axis=0)
	aux2 = np.concatenate((var2, var4), axis=0)
	j = np.random.permutation(aux1.shape[0])
	k = np.random.permutation(aux2.shape[0])
	aux11 = aux1[j, :, :]
	aux22 = aux2[k, :, :]
	var11 = aux11[0: var1.shape[0], :, :]
	var33 = aux11[var1.shape[0]: , :, :]
	var22 = aux22[0: var2.shape[0], :, :]
	var44 = aux22[var2.shape[0]: , :, :]
	varx = np.mean(var11, axis=0) - np.mean(var22, axis=0)
	vary = np.mean(var33, axis=0) - np.mean(var44, axis=0)
	varx = np.ravel(varx)
	vary = np.ravel(vary)
	correlacion = np.corrcoef(varx, vary)[0, 1]
	return correlacion

#================================================
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
startTime=datetime.now()

iter = sys.argv[1]

PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'

if ~os.path.isfile(PATH_DATA_2 + 'correlations/monthly_correlations_SPoV_enso_' + str(iter) +'_q.nc4'):
	correl_EN = np.empty([7])
	correl_LN = np.empty([7])

	month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
	seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

	for i in np.arange(0, 7):
		FILE = np.load(PATH_DATA_2 + 'z200_PV_conditioned_ENSO_' + month[i] + '_q.npz')
		#test correlation
		correl_EN[i] = TestCorrelation(FILE['var1'],  FILE['var2'], FILE['var3'], FILE['var4'])
		correl_LN[i] = TestCorrelation(FILE['var1'],  FILE['var2'], FILE['var5'], FILE['var6'])

	ds = xr.Dataset({'correl_EN': (['month'], correl_EN),
			 'correl_LN': (['month'], correl_LN)},
			 coords={'month': (['month'], month)})
	ds.to_netcdf(PATH_DATA_2 + 'correlations/monthly_correlations_SPoV_enso_' + str(iter) + '_q.nc4')
if ~os.path.isfile(PATH_DATA_2 + 'correlations/seasonal_correlations_SPoV_enso_' + str(iter) +'_q.nc4'):

	correl_EN = np.empty([5])
	correl_LN = np.empty([5])

	for i in np.arange(0, 5):
		FILE = np.load(PATH_DATA_2 + 'z200_PV_conditioned_ENSO_' + seas[i] + '_q.npz')
		#test correlation
		correl_EN[i] = TestCorrelation(FILE['var1'], FILE['var2'], FILE['var3'], FILE['var4'])
		correl_LN[i] = TestCorrelation(FILE['var1'], FILE['var2'], FILE['var5'], FILE['var6'])
	ds = xr.Dataset({'correl_EN': (['seas'], correl_EN),
			 'correl_LN': (['seas'], correl_LN)},
			 coords={'seas': (['seas'], seas)})
	ds.to_netcdf(PATH_DATA_2 + 'correlations/seasonal_correlations_SPoV_enso_' + str(iter) +'_q.nc4')

	print(datetime.now()- startTime)


