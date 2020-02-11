#correlations of extratropical asymmetries and computation of asymmetric index
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
def ComputeAsymmetry(field, pattern):
	"""project pattern inot fiels"""
	pattern = np.ravel(pattern)
	pattern = pattern #/ np.sqrt(np.sum(pattern * pattern))
	pattern = np.tile(pattern[np.newaxis, :], (field.shape[0], 1))
	field = np.reshape(field,[field.shape[0], field.shape[1]*field.shape[2]])
	field_norm = np.reshape(np.tile((np.sum(field * field, axis=1)), (1, field.shape[1])),
				[field.shape[0], field.shape[1]])
	print(field_norm.shape)
#	field_norm = 
	index = np.squeeze(np.sum(field /field_norm * pattern, axis=1))
	return index
#================================================
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/decile_new/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
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

#enso during all years
index_ninio_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear')
index_ninia_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear'), ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear'))
#enso during weak PoV

index_ninio_WPV = np.logical_and(index_ninio_all.values, index_SPV_upper.values)
index_ninia_WPV = np.logical_and(index_ninia_all.values, index_SPV_upper.values)
index_normal_WPV = np.logical_and(index_normal_all.values, index_SPV_upper.values)

#enso during strong PoV
index_ninio_SPV = np.logical_and(index_ninio_all.values, index_SPV_lower.values)
index_ninia_SPV = np.logical_and(index_ninia_all.values, index_SPV_lower.values)
index_normal_SPV = np.logical_and(index_normal_all.values, index_SPV_lower.values)

correl_ninio_WPV = np.empty([7])
correl_ninia_WPV = np.empty([7])
correl_ninio_SPV = np.empty([7])
correl_ninia_SPV = np.empty([7])
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

for i in np.arange(0, 7):
	var_ninio_WPV = hgt.z.values[i, index_ninio_WPV, :, :]
	var_normal_WPV = hgt.z.values[i, index_normal_WPV, :, :]
	var_ninia_WPV = hgt.z.values[i, index_ninia_WPV, :, :]	
	var_ninio_SPV = hgt.z.values[i, index_ninio_SPV, :, :]
	var_normal_SPV = hgt.z.values[i, index_normal_SPV, :, :]
	var_ninia_SPV = hgt.z.values[i, index_ninia_SPV, :, :]	
	var_ninio_all = hgt.z.values[i, index_ninio_all.values, :, :]
	var_normal_all = hgt.z.values[i, index_normal_all.values, :, :]
	var_ninia_all = hgt.z.values[i, index_ninia_all.values, :, :]
	np.savez(PATH_DATA_2 + 'z200_conditioned_' + month[i] + '_d_new.npz', var1=hgt.z.values[i, index_ninio_all.values, :, :], var2=hgt.z.values[i, index_normal_all.values, :, :], var3=hgt.z.values[i, index_ninia_all.values, :, :], var4=hgt.z.values[i, index_ninio_WPV, :, :], var5=hgt.z.values[i, index_normal_WPV, :, :], var6=hgt.z.values[i, index_ninia_WPV, :, :], var7=hgt.z.values[i, index_ninio_SPV, :, :], var8=hgt.z.values[i, index_normal_SPV, :, :], var9=hgt.z.values[i, index_ninia_SPV, :, :])

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
ds.to_netcdf(PATH_DATA_2 + 'monthly_correlations_enso_SPoV_polar_d_new.nc4')

correl_ninio_WPV = np.empty([5])
correl_ninia_WPV = np.empty([5])
correl_ninio_SPV = np.empty([5])
correl_ninia_SPV = np.empty([5])

for i in np.arange(0, 5):
	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
	var_ninio_WPV = hgt_s.z.values[index_ninio_WPV, :, :]
	var_normal_WPV = hgt_s.z.values[index_normal_WPV, :, :]
	var_ninia_WPV = hgt_s.z.values[index_ninia_WPV, :, :]	
	var_ninio_SPV = hgt_s.z.values[index_ninio_SPV, :, :]
	var_normal_SPV = hgt_s.z.values[index_normal_SPV, :, :]
	var_ninia_SPV = hgt_s.z.values[index_ninia_SPV, :, :]	
	var_ninio_all = hgt_s.z.values[index_ninio_all.values, :, :]
	var_normal_all = hgt_s.z.values[index_normal_all.values, :, :]
	var_ninia_all = hgt_s.z.values[index_ninia_all.values, :, :]

	#save npz file to compute correlations
	np.savez(PATH_DATA_2 + 'z200_conditioned_' + seas[i] + '_d_new.npz', var1=hgt_s.z.values[index_ninio_all.values, :, :], var2=hgt_s.z.values[index_normal_all.values, :, :], var3=hgt_s.z.values[index_ninia_all.values, :, :], var4=hgt_s.z.values[index_ninio_WPV, :, :], var5=hgt_s.z.values[index_normal_WPV, :, :], var6=hgt_s.z.values[index_ninia_WPV, :, :], var7=hgt_s.z.values[index_ninio_SPV, :, :], var8=hgt_s.z.values[index_normal_SPV, :, :], var9=hgt_s.z.values[index_ninia_SPV, :, :])
	#test correlation
	correl_ninio_WPV[i] = TestCorrelation(var_ninio_all, var_normal_all, var_ninio_WPV, var_normal_WPV)
	correl_ninia_WPV[i] = TestCorrelation(var_ninia_all, var_normal_all, var_ninia_WPV, var_normal_WPV)
	correl_ninio_SPV[i] = TestCorrelation(var_ninio_all, var_normal_all, var_ninio_SPV, var_normal_SPV)
	correl_ninia_SPV[i] = TestCorrelation(var_ninia_all, var_normal_all, var_ninia_SPV, var_normal_SPV)


ds = xr.Dataset({'correl_ninio_WPV': (['seas'], correl_ninio_WPV),
		 'correl_ninia_WPV': (['seas'], correl_ninia_WPV),
		 'correl_ninio_SPV': (['seas'], correl_ninio_SPV),
		 'correl_ninia_SPV': (['seas'], correl_ninia_SPV),},
		 coords={'seas': (['seas'], seas)})
ds.to_netcdf(PATH_DATA_2 + 'seasonal_correlations_enso_SPoV_polar_d_new.nc4')



