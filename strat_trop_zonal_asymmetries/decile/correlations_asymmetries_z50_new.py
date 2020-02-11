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

if ~os.path.isfile(PATH_DATA_2 + 'correlations/monthly_correlations_z50_enso_SPoV_' + str(iter) +'_d_new.nc4'):
	correl_ninio_WPV = np.empty([7])
	correl_ninia_WPV = np.empty([7])
	correl_ninio_SPV = np.empty([7])
	correl_ninia_SPV = np.empty([7])

	month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
	seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

	for i in np.arange(0, 7):
		FILE = np.load(PATH_DATA_2 + 'z50_conditioned_' + month[i] + '_d_new.npz')
		#test correlation
		correl_ninio_WPV[i] = TestCorrelation(FILE['var1'],  FILE['var2'], FILE['var4'], FILE['var5'])
		correl_ninia_WPV[i] = TestCorrelation(FILE['var3'],  FILE['var2'], FILE['var6'], FILE['var5'])
		correl_ninio_SPV[i] = TestCorrelation(FILE['var1'], FILE['var2'], FILE['var7'], FILE['var8'])
		correl_ninia_SPV[i] = TestCorrelation(FILE['var3'], FILE['var2'], FILE['var9'], FILE['var8'])

	ds = xr.Dataset({'correl_ninio_WPV': (['month'], correl_ninio_WPV),
			 'correl_ninia_WPV': (['month'], correl_ninia_WPV),
			 'correl_ninio_SPV': (['month'], correl_ninio_SPV),
			 'correl_ninia_SPV': (['month'], correl_ninia_SPV),},
			 coords={'month': (['month'], month), 'iter': iter })
	ds.to_netcdf(PATH_DATA_2 + 'correlations/monthly_correlations_z50_enso_SPoV_' + str(iter) + '_d_new.nc4')
if ~os.path.isfile(PATH_DATA_2 + 'correlations/seasonal_correlations_z50_enso_SPoV_' + str(iter) +'_d_new.nc4'):

	correl_ninio_WPV = np.empty([5])
	correl_ninia_WPV = np.empty([5])
	correl_ninio_SPV = np.empty([5])
	correl_ninia_SPV = np.empty([5])


	for i in np.arange(0, 5):
		FILE = np.load(PATH_DATA_2 + 'z50_conditioned_' + seas[i] + '_d_new.npz')
		#test correlation
		correl_ninio_WPV[i] = TestCorrelation(FILE['var1'], FILE['var2'], FILE['var4'], FILE['var5'])
		correl_ninia_WPV[i] = TestCorrelation(FILE['var3'], FILE['var2'], FILE['var6'], FILE['var5'])
		correl_ninio_SPV[i] = TestCorrelation(FILE['var1'], FILE['var2'], FILE['var7'], FILE['var8'])
		correl_ninia_SPV[i] = TestCorrelation(FILE['var3'], FILE['var2'], FILE['var9'], FILE['var8'])

	ds = xr.Dataset({'correl_ninio_WPV': (['seas'], correl_ninio_WPV),
			 'correl_ninia_WPV': (['seas'], correl_ninia_WPV),
			 'correl_ninio_SPV': (['seas'], correl_ninio_SPV),
			 'correl_ninia_SPV': (['seas'], correl_ninia_SPV),},
			 coords={'seas': (['seas'], seas), 'iter': iter})
	ds.to_netcdf(PATH_DATA_2 + 'correlations/seasonal_correlations_z50_enso_SPoV_' + str(iter) +'_d_new.nc4')

	print(datetime.now()- startTime)


