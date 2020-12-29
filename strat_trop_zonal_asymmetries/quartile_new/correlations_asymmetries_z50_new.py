# composites of EN events conditioned on PV strength
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
    var33 = aux11[var1.shape[0]:, :, :]
    var22 = aux22[0: var2.shape[0], :, :]
    var44 = aux22[var2.shape[0]:, :, :]
    varx = np.mean(var11, axis=0) - np.mean(var22, axis=0)
    vary = np.mean(var33, axis=0) - np.mean(var44, axis=0)
    varx = np.ravel(varx)
    vary = np.ravel(vary)
    correlacion = np.corrcoef(varx, vary)[0, 1]
    return correlacion

def ComputeCorrelation(var1, var2, var3, var4):
    varx = np.mean(var1, axis=0) - np.mean(var2, axis=0)
    vary = np.mean(var3, axis=0) - np.mean(var4, axis=0)
    varx = np.ravel(varx)
    vary = np.ravel(vary)
    correlacion = np.corrcoef(varx, vary)[0, 1]
    return correlacion

# ================================================
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
startTime = datetime.now()
iter = sys.argv[1]

PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'

if ~os.path.isfile(PATH_DATA_2 + 'correlations/monthly_correlations_z50_new_' + str(iter) + '.nc4'):
    correl_PV_cond_ENSO = np.empty([7])
    correl_ENSO_cond_PV = np.empty([7])
    correl_in_phase = np.empty([7])
    correl_out_phase = np.empty([7])
    month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']

    for i in np.arange(0, 7):
        FILE = np.load(PATH_DATA_2 + 'z50_conditioned_' + month[i] + '_new.npz')
        # var4 ninio WPV, var5 ninia WPV, var6 ninio SPV, var7 ninia SPV var2 all
        # test correlation
        correl_PV_cond_ENSO[i] = TestCorrelation(FILE['var4'], FILE['var6'], FILE['var5'], FILE['var7'])
        correl_ENSO_cond_PV[i] = TestCorrelation(FILE['var4'], FILE['var5'], FILE['var6'], FILE['var7'])
        correl_in_phase[i] = TestCorrelation(FILE['var4'], FILE['var2'], -1 * FILE['var7'], -1 * FILE['var2'])
        correl_out_phase[i] = TestCorrelation(FILE['var6'], FILE['var2'], -1 * FILE['var5'], -1 * FILE['var2'])

    ds = xr.Dataset(dict(correl_PV_cond_ENSO=(['month'], correl_PV_cond_ENSO),
                         correl_ENSO_cond_PV=(['month'], correl_ENSO_cond_PV),
                         correl_in_phase=(['month'], correl_in_phase),
                         correl_out_pase=(['month'], correl_out_phase)),
                    coords={'month': (['month'], month), 'iter': iter})
    ds.to_netcdf(PATH_DATA_2 + 'correlations/monthly_correlations_z50_new_' + str(iter) + '.nc4')

if ~os.path.isfile(PATH_DATA_2 + 'monthly_z50_new_correlations.nc4'):
    correl_PV_cond_ENSO = np.empty([7])
    correl_ENSO_cond_PV = np.empty([7])
    correl_in_phase = np.empty([7])
    correl_out_phase = np.empty([7])
    month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']

    for i in np.arange(0, 7):
        FILE = np.load(PATH_DATA_2 + 'z50_conditioned_' + month[i] + '_new.npz')
        # var4 ninio WPV, var5 ninia WPV, var6 ninio SPV, var7 ninia SPV var2 all
        # test correlation
        correl_PV_cond_ENSO[i] = ComputeCorrelation(FILE['var4'], FILE['var6'], FILE['var5'], FILE['var7'])
        correl_ENSO_cond_PV[i] = ComputeCorrelation(FILE['var4'], FILE['var5'], FILE['var6'], FILE['var7'])
        correl_in_phase[i] = ComputeCorrelation(FILE['var4'], FILE['var2'], FILE['var7'], FILE['var2'])
        correl_out_phase[i] = ComputeCorrelation(FILE['var6'], FILE['var2'], FILE['var5'], FILE['var2'])

    ds = xr.Dataset(dict(correl_PV_cond_ENSO=(['month'], correl_PV_cond_ENSO),
                         correl_ENSO_cond_PV=(['month'], correl_ENSO_cond_PV),
                         correl_in_phase=(['month'], correl_in_phase),
                         correl_out_pase=(['month'], correl_out_phase)),
                    coords={'month': (['month'], month)})
    ds.to_netcdf(PATH_DATA_2 + 'monthly_z50_new_correlations.nc4')

if ~os.path.isfile(PATH_DATA_2 + 'correlations/monthly_correlations_z200_new_' + str(iter) + '.nc4'):
    correl_PV_cond_ENSO = np.empty([7])
    correl_ENSO_cond_PV = np.empty([7])
    correl_in_phase = np.empty([7])
    correl_out_phase = np.empty([7])
    month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']

    for i in np.arange(0, 7):
        FILE = np.load(PATH_DATA_2 + 'z200_conditioned_' + month[i] + '_new.npz')
        # var4 ninio WPV, var5 ninia WPV, var6 ninio SPV, var7 ninia SPV var2 all
        # test correlation
        correl_PV_cond_ENSO[i] = TestCorrelation(FILE['var4'], FILE['var6'], FILE['var5'], FILE['var7'])
        correl_ENSO_cond_PV[i] = TestCorrelation(FILE['var4'], FILE['var5'], FILE['var6'], FILE['var7'])
        correl_in_phase[i] = TestCorrelation(FILE['var4'], FILE['var2'], -1 * FILE['var7'], -1 * FILE['var2'])
        correl_out_phase[i] = TestCorrelation(FILE['var6'], FILE['var2'], -1 * FILE['var5'], -1 * FILE['var2'])

    ds = xr.Dataset(dict(correl_PV_cond_ENSO=(['month'], correl_PV_cond_ENSO),
                         correl_ENSO_cond_PV=(['month'], correl_ENSO_cond_PV),
                         correl_in_phase=(['month'], correl_in_phase),
                         correl_out_pase=(['month'], correl_out_phase)),
                    coords={'month': (['month'], month), 'iter': iter})
    ds.to_netcdf(PATH_DATA_2 + 'correlations/monthly_correlations_z200_new_' + str(iter) + '.nc4')

if ~os.path.isfile(PATH_DATA_2 + 'monthly_z200_new_correlations.nc4'):
    correl_PV_cond_ENSO = np.empty([7])
    correl_ENSO_cond_PV = np.empty([7])
    correl_in_phase = np.empty([7])
    correl_out_phase = np.empty([7])
    month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']

    for i in np.arange(0, 7):
        FILE = np.load(PATH_DATA_2 + 'z200_conditioned_' + month[i] + '_new.npz')
        # var4 ninio WPV, var5 ninia WPV, var6 ninio SPV, var7 ninia SPV var2 all
        # test correlation
        correl_PV_cond_ENSO[i] = ComputeCorrelation(FILE['var4'], FILE['var6'], FILE['var5'], FILE['var7'])
        correl_ENSO_cond_PV[i] = ComputeCorrelation(FILE['var4'], FILE['var5'], FILE['var6'], FILE['var7'])
        correl_in_phase[i] = ComputeCorrelation(FILE['var4'], FILE['var2'], FILE['var7'], FILE['var2'])
        correl_out_phase[i] = ComputeCorrelation(FILE['var6'], FILE['var2'], FILE['var5'], FILE['var2'])

    ds = xr.Dataset(dict(correl_PV_cond_ENSO=(['month'], correl_PV_cond_ENSO),
                         correl_ENSO_cond_PV=(['month'], correl_ENSO_cond_PV),
                         correl_in_phase=(['month'], correl_in_phase),
                         correl_out_pase=(['month'], correl_out_phase)),
                    coords={'month': (['month'], month)})
    ds.to_netcdf(PATH_DATA_2 + 'monthly_z200_new_correlations.nc4')
print(datetime.now() - startTime)
