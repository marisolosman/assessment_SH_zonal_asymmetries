"""Compute asymmetric index of composites of z200 in PV events conditioned on ENSO strength
Compute correlations between composites
"""
import numpy as np
import xarray as xr
import os

def TestCorrelation(var1, var2, var3, var4):
	"""compute correlations"""
	varx = np.mean(var1, axis=0) - np.mean(var2, axis=0)
	vary = np.mean(var3, axis=0) - np.mean(var4, axis=0)
	varx = np.ravel(varx)
	vary = np.ravel(vary)
	correlacion = np.corrcoef(varx, vary)[0, 1]
	return correlacion
def ComputeAsymmetry(field, pattern):
	"""project pattern into fields"""
	pattern = np.ravel(pattern)
	pattern = pattern 
	pattern = np.tile(pattern[np.newaxis, :], (field.shape[0], 1))
	field = np.reshape(field,[field.shape[0], field.shape[1]*field.shape[2]])
	field_norm = np.reshape(np.tile(np.sqrt(np.sum(field * field, axis=1)), (1, field.shape[1])),
				[field.shape[0], field.shape[1]])
	index = np.squeeze(np.sum(field /field_norm * pattern, axis=1))
	return index
#================================================
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
hgt = hgt.sel(**{'latitude':slice(-45, -90)})
hgt = hgt - hgt.mean(dim='longitude')
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for EN years 
index_EN = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear')
#search for LN years
index_LN = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear')

#PV during all phases
index_SPV_all = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.10, dim='dim_0', interpolation='linear')
index_WPV_all = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.90, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(PV_index.SPV_index < PV_index.SPV_index.quantile(0.90, dim='dim_0', interpolation='linear'), PV_index.SPV_index > PV_index.SPV_index.quantile(0.10, dim='dim_0', interpolation='linear'))


#PV  during EN
index_SPV_EN = np.logical_and(index_SPV_all.values, index_EN.values)
index_WPV_EN = np.logical_and(index_WPV_all.values, index_EN.values)
index_normal_EN = np.logical_and(index_normal_all.values, index_EN.values)

#PV  during EN
index_SPV_LN = np.logical_and(index_SPV_all.values, index_LN.values)
index_WPV_LN = np.logical_and(index_WPV_all.values, index_LN.values)
index_normal_LN = np.logical_and(index_normal_all.values, index_LN.values)

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

correl_SPV_EN = np.empty([7])
correl_SPV_LN = np.empty([7])
correl_WPV_EN = np.empty([7])
correl_WPV_LN = np.empty([7])

for i in np.arange(0, 7):
	var_WPV_EN = hgt.z.values[i, index_WPV_EN, :, :]
	var_WPV_LN = hgt.z.values[i, index_WPV_LN, :, :]
	var_normal_EN = hgt.z.values[i, index_normal_EN, :, :]
	var_normal_LN = hgt.z.values[i, index_normal_LN, :, :]
	var_normal_all = hgt.z.values[i, index_normal_all.values, :, :]
	var_SPV_EN = hgt.z.values[i, index_SPV_EN, :, :]
	var_SPV_LN = hgt.z.values[i, index_SPV_LN, :, :]
	var_WPV_all = hgt.z.values[i, index_WPV_all.values, :, :]
	var_SPV_all = hgt.z.values[i, index_SPV_all.values, :, :]

	np.savez(PATH_DATA_2 + 'z200_PV_conditioned_ENSO_' + month[i] + '_sep_d_new.npz', var1=hgt.z.values[i, index_SPV_all.values, :, :], var2=hgt.z.values[i, index_WPV_all.values, :, :], var3=hgt.z.values[i, index_normal_all.values, :, :], var4=hgt.z.values[i, index_SPV_EN, :, :], var5=hgt.z.values[i, index_WPV_EN, :, :], var6=hgt.z.values[i, index_normal_EN, :, :], var7=hgt.z.values[i, index_SPV_LN, :, :], var8=hgt.z.values[i, index_WPV_LN, :, :], var9=hgt.z.values[i, index_normal_LN, :, :])

	#test correlation
	correl_SPV_EN[i] = TestCorrelation(var_SPV_all, var_normal_all, var_SPV_EN, var_normal_EN)
	correl_WPV_EN[i] = TestCorrelation(var_WPV_all, var_normal_all, var_WPV_EN, var_normal_EN)
	correl_SPV_LN[i] = TestCorrelation(var_SPV_all, var_normal_all, var_SPV_LN, var_normal_LN)
	correl_WPV_LN[i] = TestCorrelation(var_WPV_all, var_normal_all, var_WPV_LN, var_normal_LN)

ds = xr.Dataset({'correl_SPV_EN': (['month'], correl_SPV_EN),
		 'correl_SPV_LN': (['month'], correl_SPV_LN),
		 'correl_WPV_EN': (['month'], correl_WPV_EN),
		 'correl_WPV_LN': (['month'], correl_WPV_LN)},
		 coords={'month': (['month'], month)})
ds.to_netcdf(PATH_DATA_2 + 'monthly_correlations_SPoV_enso_polar_sep_d_new.nc4')

correl_SPV_EN = np.empty([5])
correl_SPV_LN = np.empty([5])
correl_WPV_EN = np.empty([5])
correl_WPV_LN = np.empty([5])
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
for i in np.arange(0, 5):
	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
	var_WPV_EN = hgt_s.z.values[index_WPV_EN, :, :]
	var_WPV_LN = hgt_s.z.values[index_WPV_LN, :, :]	
	var_SPV_EN = hgt_s.z.values[index_SPV_EN, :, :]
	var_SPV_LN = hgt_s.z.values[index_SPV_LN, :, :]
	var_normal_EN = hgt_s.z.values[index_normal_EN, :, :]
	var_normal_LN = hgt_s.z.values[index_normal_LN, :, :]
	var_WPV_all = hgt_s.z.values[index_WPV_all.values, :, :]
	var_SPV_all = hgt_s.z.values[index_SPV_all.values, :, :]
	var_normal_all = hgt_s.z.values[index_normal_all.values, :, :]
	#test correlation
	correl_SPV_EN[i] = TestCorrelation(var_SPV_all, var_normal_all, var_SPV_EN, var_normal_EN)
	correl_WPV_EN[i] = TestCorrelation(var_WPV_all, var_normal_all, var_WPV_EN, var_normal_EN)
	correl_SPV_LN[i] = TestCorrelation(var_SPV_all, var_normal_all, var_SPV_LN, var_normal_LN)
	correl_WPV_LN[i] = TestCorrelation(var_WPV_all, var_normal_all, var_WPV_LN, var_normal_LN)

	#save npz file to compute correlations
	np.savez(PATH_DATA_2 + 'z200_PV_conditioned_ENSO_' + seas[i] + '_sep_d_new.npz',  var1=hgt_s.z.values[index_SPV_all.values, :, :], var2=hgt_s.z.values[index_WPV_all.values, :, :], var3=hgt_s.z.values[index_normal_all.values, :, :],  var4=hgt_s.z.values[index_SPV_EN, :, :], var5=hgt_s.z.values[index_WPV_EN, :, :], var6=hgt_s.z.values[index_normal_EN, :, :],  var7=hgt_s.z.values[index_SPV_LN, :, :], var8=hgt_s.z.values[index_WPV_LN, :, :], var9=hgt_s.z.values[index_normal_LN, :, :])

ds = xr.Dataset({'correl_SPV_EN': (['seas'], correl_SPV_EN),
		 'correl_SPV_LN': (['seas'], correl_SPV_LN),
		 'correl_WPV_EN': (['seas'], correl_WPV_EN),
		 'correl_WPV_LN': (['seas'], correl_WPV_LN)},
		 coords={'seas': (['seas'], seas)})
ds.to_netcdf(PATH_DATA_2 + 'seasonal_correlations_SPoV_enso_polar_sep_d_new.nc4')



