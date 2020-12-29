#composites of EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import plots_stereo as plots
import matplotlib.pyplot as plt
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures_paper/'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for years with weak PV
index_SPV_upper = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')
#search for years with normal PV
index_SPV_normal = np.logical_and(PV_index.SPV_index > PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear'), PV_index.SPV_index < PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear'))

#enso during all years
index_ninio_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
index_ninia_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'))

#enso during weak PoV
index_ninio_WPV = np.logical_and(index_ninio_all.values, index_SPV_upper.values)
index_ninia_WPV = np.logical_and(index_ninia_all.values, index_SPV_upper.values)

#enso during strong PoV
index_ninio_SPV = np.logical_and(index_ninio_all.values, index_SPV_lower.values)
index_ninia_SPV = np.logical_and(index_ninia_all.values, index_SPV_lower.values)

index_normal = np.logical_and(index_normal_all.values, index_SPV_normal.values)
index_ninio_normal = np.logical_and(index_ninio_all.values, index_SPV_normal.values)
index_ninia_normal = np.logical_and(index_ninia_all.values, index_SPV_normal.values)
index_SPV_lower_normal = np.logical_and(index_normal_all.values, index_SPV_lower.values)
index_SPV_upper_normal = np.logical_and(index_normal_all.values, index_SPV_upper.values)

print(sum(index_ninio_all), sum(index_ninia_all))
print(sum(index_ninio_WPV), sum(index_ninia_WPV))
print(sum(index_ninio_SPV), sum(index_ninia_SPV))
print(sum(index_normal))
print(sum(index_ninio_normal),sum(index_ninia_normal) )
print(sum(index_SPV_upper), sum(index_SPV_lower))
print(sum(index_SPV_upper_normal), sum(index_SPV_lower_normal))

fig = plt.figure(figsize=(10,10), dpi=300)

plt.scatter(ninio34.ninio34_index[index_normal], PV_index.SPV_index[index_normal], color='black',
	   label='Neutral')
plt.scatter(ninio34.ninio34_index[index_SPV_upper_normal], PV_index.SPV_index[index_SPV_upper_normal],
	    c='tab:red', marker='v', label='Weak SPV only')
plt.scatter(ninio34.ninio34_index[index_SPV_lower_normal], PV_index.SPV_index[index_SPV_lower_normal],
	    c='tab:blue', marker='v', label='Strong SPV only')
plt.scatter(ninio34.ninio34_index[index_ninio_normal], PV_index.SPV_index[index_ninio_normal],
	    c='tab:green', marker='o', label='Ninio only')
plt.scatter(ninio34.ninio34_index[index_ninio_WPV], PV_index.SPV_index[index_ninio_WPV],
	    c='tab:red', marker='o', label='Ninio and Weak SPV')
plt.scatter(ninio34.ninio34_index[index_ninio_SPV], PV_index.SPV_index[index_ninio_SPV],
	    c='tab:blue', marker='o', label='Ninio and Strong SPV')
plt.scatter(ninio34.ninio34_index[index_ninia_normal], PV_index.SPV_index[index_ninia_normal],
	    c='tab:green', marker='d', label='Ninia only')
plt.scatter(ninio34.ninio34_index[index_ninia_WPV], PV_index.SPV_index[index_ninia_WPV],
	    c='tab:red', marker='d', label='Ninia and Weak SPV')
plt.scatter(ninio34.ninio34_index[index_ninia_SPV], PV_index.SPV_index[index_ninia_SPV],
	    c='tab:blue', marker='d', label='Ninia and Strong SPV')
plt.axhline(0, color='black')
plt.axvline(0, color='black')
plt.axhline(PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear'), linestyle='--')
plt.axhline(PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear'), linestyle='--')
plt.axvline(ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'), linestyle='--')
plt.axvline(ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'), linestyle='--')
plt.xlabel('Ninio 3.4', fontsize=12)
plt.ylabel('SPV index', fontsize=12)
plt.xlim([-10.2, 10.2])
plt.tick_params(axis="x", labelsize=11)
plt.tick_params(axis="y", labelsize=11)
plt.legend(prop={'size': 12})
plt.ylim([-1200, 1200])
plt.title('Scatter SPV vs Ninio 3.4 - S4', fontsize=13)
plt.tight_layout()
plt.savefig(FIG_PATH + 'scatter_ninio_spov_s4.eps', dpi=300)
