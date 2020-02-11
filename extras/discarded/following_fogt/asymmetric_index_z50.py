#correlations of extratropical asymmetries fo z* at 50hPa
import numpy as np
import xarray as xr
import os

def TestCorrelation(var1, var2, var3, var4):
	varx = np.mean(var1, axis=0) - np.mean(var2, axis=0)
	vary = np.mean(var3, axis=0) - np.mean(var4, axis=0)
	varx = np.ravel(varx)
	vary = np.ravel(vary)
	correlacion = np.corrcoef(varx, vary)[0, 1]
	return correlacion
#================================================
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/'
FILE_HGT_S4 = 'monthly_hgt50_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
hgt = hgt - hgt.mean(dim='longitude')
hgt = hgt.sel(**{'latitude':slice(-45, -90)})
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for years with weak PV
index_SPV_upper = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.90, dim='dim_0', interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.10, dim='dim_0', interpolation='linear')

# compute EN-LA composites conditioned on PV anomalies
ninio34_WPV = ninio34.sel(dim_0 = index_SPV_upper.values)
hgt_WPV = hgt.sel(realiz = index_SPV_upper.values)
ninio34_SPV = ninio34.sel(dim_0 = index_SPV_lower.values)
hgt_SPV = hgt.sel(realiz = index_SPV_lower.values)

#enso during all years
index_ninio_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear')
index_ninia_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear'))
#enso during weak PoV
index_ninio_WPV = ninio34_WPV.ninio34_index >= ninio34_WPV.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear')
index_ninia_WPV = ninio34_WPV.ninio34_index <= ninio34_WPV.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear')
index_normal_WPV = np.logical_and(ninio34_WPV.ninio34_index < ninio34_WPV.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_WPV.ninio34_index > ninio34_WPV.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear'))
#enso during strong PoV
index_ninio_SPV = ninio34_SPV.ninio34_index >= ninio34_SPV.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear')
index_ninia_SPV = ninio34_SPV.ninio34_index <= ninio34_SPV.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear')
index_normal_SPV = np.logical_and(ninio34_SPV.ninio34_index < ninio34_SPV.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34_SPV.ninio34_index > ninio34_SPV.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear'))

correl_ninio_WPV = np.empty([7])
correl_ninia_WPV = np.empty([7])
correl_ninio_SPV = np.empty([7])
correl_ninia_SPV = np.empty([7])
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

for i in np.arange(0, 7):
	var_ninio_WPV = np.mean(hgt_WPV.z.values[i, index_ninio_WPV.values, :, :], axis=0)
	var_normal_WPV = np.mean(hgt_WPV.z.values[i, index_normal_WPV.values, :, :], axis=0)
	var_ninia_WPV = np.mean(hgt_WPV.z.values[i, index_ninia_WPV.values, :, :], axis=0)	
	var_ninio_SPV = np.mean(hgt_SPV.z.values[i, index_ninio_SPV.values, :, :], axis=0)
	var_normal_SPV = np.mean(hgt_SPV.z.values[i, index_normal_SPV.values, :, :], axis=0)
	var_ninia_SPV = np.mean(hgt_SPV.z.values[i, index_ninia_SPV.values, :, :], axis=0)	
	var_ninio_all = np.mean(hgt.z.values[i, index_ninio_all.values, :, :], axis=0)
	var_normal_all = np.mean(hgt.z.values[i, index_normal_all.values, :, :], axis=0)
	var_ninia_all = np.mean(hgt.z.values[i, index_ninia_all.values, :, :], axis=0)	
	np.savez(PATH_DATA_2 + 'z50_conditioned_' + month[i] + '.npz', var1=hgt.z.values[i, index_ninio_all.values, :, :], var2=hgt.z.values[i, index_normal_all.values, :, :], var3=hgt.z.values[i, index_ninia_all.values, :, :], var4=hgt_WPV.z.values[i, index_ninio_WPV.values, :, :], var5=hgt_WPV.z.values[i, index_normal_WPV.values, :, :], var6=hgt_WPV.z.values[i, index_ninia_WPV.values, :, :], var7=hgt_SPV.z.values[i, index_ninio_SPV.values, :, :], var8=hgt_SPV.z.values[i, index_normal_SPV.values, :, :], var9=hgt_SPV.z.values[i, index_ninia_SPV.values, :, :])

#testear si los campos son distintos:
	#test correlation
	correl_ninio_WPV[i] = TestCorrelation(var_ninio_all, var_normal_all, var_ninio_WPV, var_normal_WPV)
	correl_ninia_WPV[i] = TestCorrelation(var_ninia_all, var_normal_all, var_ninia_WPV, var_normal_WPV)
	correl_ninio_SPV[i] = TestCorrelation(var_ninio_all, var_normal_all, var_ninio_SPV, var_normal_SPV)
	correl_ninia_SPV[i] = TestCorrelation(var_ninia_all, var_normal_all, var_ninia_SPV, var_normal_SPV)


ds = xr.Dataset({'correl_ninio_WPV': (['month'], correl_ninio_WPV),
		 'correl_ninia_WPV': (['month'], correl_ninia_WPV),
		 'correl_ninio_SPV': (['month'], correl_ninio_SPV),
		 'correl_ninia_SPV': (['month'], correl_ninia_SPV),},
		 coords={'month': (['month'], month)})
ds.to_netcdf(PATH_DATA_2 + 'monthly_correlations_z50_enso_SPoV_polar.nc4')

correl_ninio_WPV = np.empty([5])
correl_ninia_WPV = np.empty([5])
correl_ninio_SPV = np.empty([5])
correl_ninia_SPV = np.empty([5])

for i in np.arange(0, 5):
	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
	hgt_s_WPV = hgt_s.sel(realiz=index_SPV_upper.values)
	hgt_s_SPV = hgt_s.sel(realiz=index_SPV_lower.values)
	var_ninio_WPV = np.mean(hgt_s_WPV.z.values[index_ninio_WPV.values, :, :], axis=0)
	var_normal_WPV = np.mean(hgt_s_WPV.z.values[index_normal_WPV.values, :, :], axis=0)
	var_ninia_WPV = np.mean(hgt_s_WPV.z.values[index_ninia_WPV.values, :, :], axis=0)	
	var_ninio_SPV = np.mean(hgt_s_SPV.z.values[index_ninio_SPV.values, :, :], axis=0)
	var_normal_SPV = np.mean(hgt_s_SPV.z.values[index_normal_SPV.values, :, :], axis=0)
	var_ninia_SPV = np.mean(hgt_s_SPV.z.values[index_ninia_SPV.values, :, :], axis=0)	
	var_ninio_all = np.mean(hgt_s.z.values[index_ninio_all.values, :, :], axis=0)
	var_normal_all = np.mean(hgt_s.z.values[index_normal_all.values, :, :], axis=0)
	var_ninia_all = np.mean(hgt_s.z.values[index_ninia_all.values, :, :], axis=0)	
	#save npz file to compute correlations
	np.savez(PATH_DATA_2 + 'z50_conditioned_' + seas[i] + '.npz', var1=hgt_s.z.values[index_ninio_all.values, :, :], var2=hgt_s.z.values[index_normal_all.values, :, :], var3=hgt_s.z.values[index_ninia_all.values, :, :], var4=hgt_s_WPV.z.values[index_ninio_WPV.values, :, :], var5=hgt_s_WPV.z.values[index_normal_WPV.values, :, :], var6=hgt_s_WPV.z.values[index_ninia_WPV.values, :, :], var7=hgt_s_SPV.z.values[index_ninio_SPV.values, :, :], var8=hgt_s_SPV.z.values[index_normal_SPV.values, :, :], var9=hgt_s_SPV.z.values[index_ninia_SPV.values, :, :])
	correl_ninio_WPV[i] = TestCorrelation(hgt_s.z.values[index_ninio_all.values, :, :], hgt_s.z.values[index_normal_all.values, :, :], hgt_s_WPV.z.values[index_ninio_WPV.values, :, :], hgt_s_WPV.z.values[index_normal_WPV.values, :, :])
	correl_ninia_WPV[i] = TestCorrelation(hgt_s.z.values[index_ninia_all.values, :, :], hgt_s.z.values[index_normal_all.values, :, :], hgt_s_WPV.z.values[index_ninia_WPV.values, :, :], hgt_s_WPV.z.values[index_normal_WPV.values, :, :])
	correl_ninio_SPV[i] = TestCorrelation(hgt_s.z.values[index_ninio_all.values, :, :], hgt_s.z.values[index_normal_all.values, :, :], hgt_s_SPV.z.values[index_ninio_SPV.values, :, :], hgt_s_SPV.z.values[index_normal_SPV.values, :, :])
	correl_ninia_SPV[i] = TestCorrelation(hgt_s.z.values[index_ninia_all.values, :, :], hgt_s.z.values[index_normal_all.values, :, :], hgt_s_SPV.z.values[index_ninia_SPV.values, :, :], hgt_s_SPV.z.values[index_normal_SPV.values, :, :])


ds = xr.Dataset({'correl_ninio_WPV': (['seas'], correl_ninio_WPV),
		 'correl_ninia_WPV': (['seas'], correl_ninia_WPV),
		 'correl_ninio_SPV': (['seas'], correl_ninio_SPV),
		 'correl_ninia_SPV': (['seas'], correl_ninia_SPV),},
		 coords={'seas': (['seas'], seas)})
ds.to_netcdf(PATH_DATA_2 + 'seasonal_correlations_z50_enso_SPoV_polar.nc4')



