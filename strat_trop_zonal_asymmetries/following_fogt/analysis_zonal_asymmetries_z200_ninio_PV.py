#Compute asymmetric index and correlations between compostes of ENSO conditioned or  not to the intensity of the Strat Polar Vortex
import numpy as np
import xarray as xr
import os
import plots
import random

def TestCorrelation(var1, var2):
	var1 = np.ravel(var1)
	var2 = np.ravel(var2)
	correlacion = np.corrcoef(var1, var2)[0, 1]
	distrib = np.empty([10000])
	for i in range(10000):
		j = [random.randint(0, var1.shape[0] - 1) for k in range(var1.shape[0])]
		distrib[i] = np.corrcoef(var1[j], var2[j])[0, 1]
	distrib = np.sort(distrib)
	return correlacion, distrib[499], distrib[9499]
def ComputeAsymmetry(field, pattern):
	pattern = np.ravel(pattern)
	pattern = np.tile(pattern[np.newaxis, :], (field.shape[0], 1))
	field = np.reshape(field,[field.shape[0], field.shape[1]*field.shape[2]])
	index = np.squeeze(np.sum(field * pattern, axis=1))
	return index
#================================================
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '~/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
hgt = hgt - hgt.mean(dim='longitude')
hgt = hgt.sel(**{'latitude':slice(-20, -90)})
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)


#search for years with EN
index_ENSO_upper = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.90, dim='dim_0', interpolation='linear')
#search for years with LN
index_ENSO_lower = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.10, dim='dim_0', interpolation='linear')

# compute EN-LA composites conditioned on PV anomalies
PoV_ninio = PV_index.sel(dim_0 = index_ENSO_upper.values)
hgt_ninio = hgt.sel(realiz = index_ENSO_upper.values)
PoV_ninia = PV_index.sel(dim_0 = index_ENSO_lower.values)
hgt_ninia = hgt.sel(realiz = index_ENSO_lower.values)

#enso during all years
WPoV_all = PV_index.PV_index >= PV_index.PV_index.quantile(0.90, dim='dim_0', interpolation='linear')
SPoV_all = PV_index.PV_index <= PV_index.PV_index.quantile(0.10, dim='dim_0', interpolation='linear')
NPoV_all = np.logical_and(PV_index.PV_index < PV_index.PV_index.quantile(0.90, dim='dim_0', interpolation='linear'), PV_index.PV_index > PV_index.PV_index.quantile(0.10, dim='dim_0', interpolation='linear'))
#enso during weak PoV
WPoV_ninio = PoV_ninio.PV_index >= PoV_ninio.PV_index.quantile(0.90, dim='dim_0', interpolation='linear')
SPoV_ninio = PoV_ninio.PV_index <= PoV_ninio.PV_index.quantile(0.10, dim='dim_0', interpolation='linear')
NPoV_ninio = np.logical_and(PoV_ninio.PV_index < PoV_ninio.PV_index.quantile(0.90, dim='dim_0', interpolation='linear'), PoV_ninio.PV_index > PoV_ninio.PV_index.quantile(0.10, dim='dim_0', interpolation='linear'))
#enso during strong PoV
WPoV_ninia = PoV_ninia.PV_index >= PoV_ninia.PV_index.quantile(0.90, dim='dim_0', interpolation='linear')
SPoV_ninia = PoV_ninia.PV_index <= PoV_ninia.PV_index.quantile(0.10, dim='dim_0', interpolation='linear')
NPoV_ninia = np.logical_and(PoV_ninia.PV_index < PoV_ninia.PV_index.quantile(0.90, dim='dim_0', interpolation='linear'), PoV_ninia.PV_index > PoV_ninia.PV_index.quantile(0.10, dim='dim_0', interpolation='linear'))

asimmetry_WPoV_ninio = np.empty([7, np.shape(hgt.realiz)[0]])
asimmetry_SPoV_ninio = np.empty([7, np.shape(hgt.realiz)[0]])
asimmetry_WPoV_ninia = np.empty([7, np.shape(hgt.realiz)[0]])
asimmetry_SPoV_ninia = np.empty([7, np.shape(hgt.realiz)[0]])
asimmetry_WPoV_all = np.empty([7, np.shape(hgt.realiz)[0]])
asimmetry_SPoV_all = np.empty([7, np.shape(hgt.realiz)[0]])

correl_WPoV_ninio = np.empty([7])
correl_SPoV_ninio = np.empty([7])
correl_WPoV_ninia = np.empty([7])
correl_SPoV_ninia = np.empty([7])
quant1_WPoV_ninio = np.empty([7])
quant1_SPoV_ninio = np.empty([7])
quant1_WPoV_ninia = np.empty([7])
quant1_SPoV_ninia = np.empty([7])
quant2_WPoV_ninio = np.empty([7])
quant2_SPoV_ninio = np.empty([7])
quant2_WPoV_ninia = np.empty([7])
quant2_SPoV_ninia = np.empty([7])

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

for i in np.arange(0, 7):
	var_WPoV_ninio = np.mean(hgt_WPV.z.values[i, index_WPoV_ninio.values, :, :], axis=0)
	var_NPoV_ninio = np.mean(hgt_WPV.z.values[i, index_NPoV_ninio.values, :, :], axis=0)
	var_SPoV_ninio = np.mean(hgt_WPV.z.values[i, index_SPoV_ninio.values, :, :], axis=0)	
	var_WPoV_ninia = np.mean(hgt_SPV.z.values[i, index_WPoV_ninia.values, :, :], axis=0)
	var_NPoV_ninia = np.mean(hgt_SPV.z.values[i, index_NPoV_ninia.values, :, :], axis=0)
	var_SPoV_ninia = np.mean(hgt_SPV.z.values[i, index_SPoV_ninia.values, :, :], axis=0)	
	var_WPoV_all = np.mean(hgt.z.values[i, index_WPoV_all.values, :, :], axis=0)
	var_NPoV_all = np.mean(hgt.z.values[i, index_NPoV_all.values, :, :], axis=0)
	var_SPoV_all = np.mean(hgt.z.values[i, index_SPoV_all.values, :, :], axis=0)	
	#test correlation
	[correl_ninio_WPV[i], quant1_ninio_WPV[i], quant2_ninio_WPV[i]] = TestCorrelation(var_ninio_all-var_normal_all, var_ninio_WPV-var_normal_WPV)
	[correl_ninia_WPV[i], quant1_ninia_WPV[i], quant2_ninia_WPV[i]] = TestCorrelation(var_ninia_all-var_normal_all, var_ninia_WPV-var_normal_WPV)
	#calcular indice de asimetria proyectando cada anio en la anomalia para obtener un valor 
	asimmetry_ninio_WPV[i, :] = ComputeAsymmetry(hgt.z.values[i, :, :, :], var_ninio_WPV - var_normal_WPV)
	asimmetry_ninia_WPV[i, :] = ComputeAsymmetry(hgt.z.values[i, :, :, :], var_ninia_WPV - var_normal_WPV)
	asimmetry_ninio_SPV[i, :] = ComputeAsymmetry(hgt.z.values[i, :, :, :], var_ninio_SPV - var_normal_SPV)
	asimmetry_ninia_SPV[i, :] = ComputeAsymmetry(hgt.z.values[i, :, :, :], var_ninia_SPV - var_normal_SPV)
	asimmetry_ninio_all[i, :] = ComputeAsymmetry(hgt.z.values[i, :, :, :], var_ninio_all - var_normal_all)
	asimmetry_ninia_all[i, :] = ComputeAsymmetry(hgt.z.values[i, :, :, :], var_ninia_all - var_normal_all)

ds = xr.Dataset({'asimm_ninio_WPV': (['month', 'realiz'], asimmetry_ninio_WPV),
		 'asimm_ninia_WPV': (['month', 'realiz'], asimmetry_ninia_WPV),
		 'asimm_ninio_SPV': (['month', 'realiz'], asimmetry_ninio_SPV),
		 'asimm_ninia_SPV': (['month', 'realiz'], asimmetry_ninia_SPV),
		 'asimm_ninio_all': (['month', 'realiz'], asimmetry_ninio_all),
		 'asimm_ninia_all': (['month', 'realiz'], asimmetry_ninia_all)},
		 coords={'month': (['month'], month),
		 	 'realiz': (['realiz'], hgt.realiz.values)})
ds.to_netcdf(PATH_DATA_2 + 'monthly_asimmetric_index_enso_SPoV.nc4')
ds = xr.Dataset({'correl_ninio_WPV': (['month'], correl_ninio_WPV),
		 'correl_ninia_WPV': (['month'], correl_ninia_WPV),
		 'correl_ninio_SPV': (['month'], correl_ninio_SPV),
		 'correl_ninia_SPV': (['month'], correl_ninia_SPV),
		 'quant1_ninio_WPV': (['month'], quant1_ninio_WPV),
		 'quant1_ninia_WPV': (['month'], quant1_ninia_WPV),
		 'quant1_ninio_SPV': (['month'], quant1_ninio_SPV),
		 'quant1_ninia_SPV': (['month'], quant1_ninia_SPV),
		 'quant2_ninio_WPV': (['month'], quant2_ninio_WPV),
		 'quant2_ninia_WPV': (['month'], quant2_ninia_WPV),
		 'quant2_ninio_SPV': (['month'], quant2_ninio_SPV),
		 'quant2_ninia_SPV': (['month'], quant2_ninia_SPV),},
		 coords={'month': (['month'], month)})
ds.to_netcdf(PATH_DATA_2 + 'monthly_correlations_enso_SPoV.nc4')

asimmetry_ninio_WPV = np.empty([5, np.shape(hgt.realiz)[0]])
asimmetry_ninia_WPV = np.empty([5, np.shape(hgt.realiz)[0]])
asimmetry_ninio_SPV = np.empty([5, np.shape(hgt.realiz)[0]])
asimmetry_ninia_SPV = np.empty([5, np.shape(hgt.realiz)[0]])
asimmetry_ninio_all = np.empty([5, np.shape(hgt.realiz)[0]])
asimmetry_ninia_all = np.empty([5, np.shape(hgt.realiz)[0]])
correl_ninio_WPV = np.empty([5])
correl_ninia_WPV = np.empty([5])
correl_ninio_SPV = np.empty([5])
correl_ninia_SPV = np.empty([5])
quant1_ninio_WPV = np.empty([5])
quant1_ninia_WPV = np.empty([5])
quant1_ninio_SPV = np.empty([5])
quant1_ninia_SPV = np.empty([5])
quant2_ninio_WPV = np.empty([5])
quant2_ninia_WPV = np.empty([5])
quant2_ninio_SPV = np.empty([5])
quant2_ninia_SPV = np.empty([5])


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

	asimmetry_ninio_WPV[i, :] = ComputeAsymmetry(hgt_s.z.values, var_ninio_WPV - var_normal_WPV)
	asimmetry_ninia_WPV[i, :] = ComputeAsymmetry(hgt_s.z.values, var_ninia_WPV - var_normal_WPV)
	asimmetry_ninio_SPV[i, :] = ComputeAsymmetry(hgt_s.z.values, var_ninio_SPV - var_normal_SPV)
	asimmetry_ninia_SPV[i, :] = ComputeAsymmetry(hgt_s.z.values, var_ninia_SPV - var_normal_SPV)
	asimmetry_ninio_all[i, :] = ComputeAsymmetry(hgt_s.z.values, var_ninio_all - var_normal_all)
	asimmetry_ninia_all[i, :] = ComputeAsymmetry(hgt_s.z.values, var_ninia_all - var_normal_all)
	[correl_ninio_WPV[i], quant1_ninio_WPV[i], quant2_ninio_WPV[i]] = TestCorrelation(var_ninio_all-var_normal_all, var_ninio_WPV-var_normal_WPV)
	[correl_ninia_WPV[i], quant1_ninia_WPV[i], quant2_ninia_WPV[i]] = TestCorrelation(var_ninia_all-var_normal_all, var_ninia_WPV-var_normal_WPV)
	[correl_ninio_SPV[i], quant1_ninio_SPV[i], quant2_ninio_SPV[i]] = TestCorrelation(var_ninio_all-var_normal_all, var_ninio_SPV-var_normal_SPV)
	[correl_ninia_SPV[i], quant1_ninia_SPV[i], quant2_ninia_SPV[i]] = TestCorrelation(var_ninia_all-var_normal_all, var_ninia_SPV-var_normal_SPV)



ds = xr.Dataset({'asimm_ninio_WPV': (['seas', 'realiz'], asimmetry_ninio_WPV),
		 'asimm_ninia_WPV': (['seas', 'realiz'], asimmetry_ninia_WPV),
		 'asimm_ninio_SPV': (['seas', 'realiz'], asimmetry_ninio_SPV),
		 'asimm_ninia_SPV': (['seas', 'realiz'], asimmetry_ninia_SPV),
		 'asimm_ninio_all': (['seas', 'realiz'], asimmetry_ninio_all),
		 'asimm_ninia_all': (['seas', 'realiz'], asimmetry_ninia_all)},
		 coords={'seas': (['seas'], seas),
		 	 'realiz': (['realiz'], hgt.realiz.values)})

ds.to_netcdf(PATH_DATA_2 + 'seasonal_asimmetric_index_enso_SPoV.nc4')

ds = xr.Dataset({'correl_ninio_WPV': (['seas'], correl_ninio_WPV),
		 'correl_ninia_WPV': (['seas'], correl_ninia_WPV),
		 'correl_ninio_SPV': (['seas'], correl_ninio_SPV),
		 'correl_ninia_SPV': (['seas'], correl_ninia_SPV),
		 'quant1_ninio_WPV': (['seas'], quant1_ninio_WPV),
		 'quant1_ninia_WPV': (['seas'], quant1_ninia_WPV),
		 'quant1_ninio_SPV': (['seas'], quant1_ninio_SPV),
		 'quant1_ninia_SPV': (['seas'], quant1_ninia_SPV),
		 'quant2_ninio_WPV': (['seas'], quant2_ninio_WPV),
		 'quant2_ninia_WPV': (['seas'], quant2_ninia_WPV),
		 'quant2_ninio_SPV': (['seas'], quant2_ninio_SPV),
		 'quant2_ninia_SPV': (['seas'], quant2_ninia_SPV),},
		 coords={'seas': (['seas'], seas)})
ds.to_netcdf(PATH_DATA_2 + 'seasonal_correlations_enso_SPoV.nc4')




