#determines location of ASL during EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
from metpy import calc

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/decile_new/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
hgt = hgt - hgt.mean(dim='longitude')
#retain hgt over ASL region
hgt = hgt.sel(**{'latitude': slice(-50, -80), 'longitude': slice(165, 355)})
#get gradient
grad = calc.gradient(hgt.z, axes=[2, 3], coordinates=[hgt.latitude.values, hgt.longitude.values])
#Lap_hgt = calc.laplacian(hgt.z, axes=[2, 3], coordinates=[hgt.latitude.values, hgt.longitude.values])

ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for years with weak PV
index_SPV_upper = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_SPV_normal = np.logical_and(PV_index.SPV_index > PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear'), PV_index.SPV_index < PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear'))


# compute EN-LA composites conditioned on PV anomalies
ninio34_WPV = ninio34.sel(dim_0 = index_SPV_upper.values)
hgt_WPV = hgt.sel(realiz = index_SPV_upper.values)
ninio34_SPV = ninio34.sel(dim_0 = index_SPV_lower.values)
hgt_SPV = hgt.sel(realiz = index_SPV_lower.values)

#enso during all years
index_ninio_all = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
index_ninia_all = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')
index_normal_all = np.logical_and(ninio34.ninio34_index < ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear'), ninio34.ninio34_index > ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear'))
#enso during weak PoV

index_ninio_WPV = np.logical_and(index_ninio_all.values, index_SPV_upper.values)
index_ninia_WPV = np.logical_and(index_ninia_all.values, index_SPV_upper.values)
index_normal_WPV = np.logical_and(index_normal_all.values, index_SPV_upper.values)

#enso during strong PoV
index_ninio_SPV = np.logical_and(index_ninio_all.values, index_SPV_lower.values)
index_ninia_SPV = np.logical_and(index_ninia_all.values, index_SPV_lower.values
index_normal_SPV = np.logical_and(index_normal_all.values, index_SPV_lower.values)
index_normal = np.logical_and(index_normal_all.values, index_SPV_normal.values)

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
lat_min = np.empty([7, len(hgt.realiz.values)])
lon_min = np.empty([7, len(hgt.realiz.values)])
lat_max = np.empty([7, len(hgt.realiz.values)])
lon_max = np.empty([7, len(hgt.realiz.values)])
val_min = np.empty([7, len(hgt.realiz.values)])
val_max = np.empty([7, len(hgt.realiz.values)])

for i in np.arange(0, 7):
	for j in np.arange(0, len(hgt.realiz.values)):
		aux = np.where(hgt.z.values[i, j, :, :] == np.min(hgt.z.values[i, j, :, :]))
		if (len(aux[0][:]) > 1) print(i, j):
		lat_min[i, j] = hgt.latitude.values[aux[0][0]]
		lon_min[i, j] = hgt.longitude.values[aux[1][0]]
		val_min[i, j] = hgt.z.values[i, j, aux[0][0], aux[1][0]]
		aux = np.where(hgt.z.values[i, j, :, :] == np.max(hgt.z.values[i, j, :, :]))
		lat_max[i, j] = hgt.latitude.values[aux[0][0]]
		lon_max[i, j] = hgt.longitude.values[aux[1][0]]
		val_max[i, j] = hgt.z.values[i, j, aux[0][0], aux[1][0]]
		if (len(aux[0][:]) > 1): print(i, j)
	var_ninio_WPV = np.array([val_min[i, index_ninio_WPV], lat_min[i, index_ninio_WPV],
				  lon_min[i, index_ninio_WPV], val_max[i, index_ninio_WPV],
				  lat_max[i, index_ninio_WPV], lon_max[i, index_ninio_WPV]])
	var_normal_WPV = np.array([val_min[i, index_normal_WPV], lat_min[i, index_normal_WPV],
				  lon_min[i, index_normal_WPV], val_max[i, index_normal_WPV],
				  lat_max[i, index_normal_WPV], lon_max[i, index_normal_WPV]])
	var_ninia_WPV =  np.array([val_min[i, index_ninia_WPV], lat_min[i, index_ninia_WPV],
				  lon_min[i, index_ninia_WPV], val_max[i, index_ninia_WPV],
				  lat_max[i, index_ninia_WPV], lon_max[i, index_ninia_WPV]])
	var_ninio_SPV = np.array([val_min[i, index_ninio_SPV], lat_min[i, index_ninio_SPV],
				  lon_min[i, index_ninio_SPV], val_max[i, index_ninio_SPV],
				  lat_max[i, index_ninio_SPV], lon_max[i, index_ninio_SPV]])
	var_normal_SPV = np.array([val_min[i, index_normal_SPV], lat_min[i, index_normal_SPV],
				  lon_min[i, index_normal_SPV], val_max[i, index_normal_SPV],
				  lat_max[i, index_normal_SPV], lon_max[i, index_normal_SPV]])
	var_ninia_SPV =  np.array([val_min[i, index_ninia_SPV], lat_min[i, index_ninia_SPV],
				  lon_min[i, index_ninia_SPV], val_max[i, index_ninia_SPV],
				  lat_max[i, index_ninia_SPV], lon_max[i, index_ninia_SPV]])
	var_ninio_all = np.array([val_min[i, index_ninio_all], lat_min[i, index_ninio_all],
				  lon_min[i, index_ninio_all], val_max[i, index_ninio_all],
				  lat_max[i, index_ninio_all], lon_max[i, index_ninio_all]])
	var_normal_all = np.array([val_min[i, index_normal_all], lat_min[i, index_normal_all],
				  lon_min[i, index_normal_all], val_max[i, index_normal_all],
				  lat_max[i, index_normal_all], lon_max[i, index_normal_all]])
	var_ninia_all =  np.array([val_min[i, index_ninia_all], lat_min[i, index_ninia_all],
				  lon_min[i, index_ninia_all], val_max[i, index_ninia_all],
				  lat_max[i, index_ninia_all], lon_max[i, index_ninia_all]])
	np.savez(PATH_DATA_2 + 'ASL_coords' + month[i] + '.npz', var1=var_ninio_WPV,
		 var2=var_normal_WPV, var3=var_ninia_WPV, var4=var_ninio_SPV, var5=var_normal_SPV,
		 var6=var_ninia_SPV, var7=var_ninio_all, var8=var_normal_all, var9=var_ninia_all)

lat_min = np.empty([5, len(hgt.realiz.values)])
lon_min = np.empty([5, len(hgt.realiz.values)])
lat_max = np.empty([5, len(hgt.realiz.values)])
lon_max = np.empty([5, len(hgt.realiz.values)])
val_min = np.empty([5, len(hgt.realiz.values)])
val_max = np.empty([5, len(hgt.realiz.values)])
for i in np.arange(0, 5):
	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
	#get laplacian
	#Lap_hgt = calc.laplacian(hgt_s.z, axes=[1, 2],
#				 coordinates=[hgt.latitude.values, hgt.longitude.values])
	for j in np.arange(0, len(hgt.realiz.values)):
		aux = np.where(hgt_s.z.values[j, :, :] == np.min(hgt_s.z.values[j, :, :]))
		lat_min[i, j] = hgt.latitude.values[aux[0][0]]
		lon_min[i, j] = hgt.longitude.values[aux[1][0]]
		val_min[i, j] = hgt_s.z.values[j, aux[0][0], aux[1][0]]
		if (len(aux[0][:]) > 1): print(i, j)
		aux = np.where(hgt_s.z.values[j, :, :] == np.max(hgt_s.z.values[j, :, :]))
		lat_max[i, j] = hgt.latitude.values[aux[0][0]]
		lon_max[i, j] = hgt.longitude.values[aux[1][0]]
		val_max[i, j] = hgt_s.z.values[j, aux[0][0], aux[1][0]]
		if (len(aux[0][:]) > 1): print(i, j)
	var_ninio_WPV = np.array([val_min[i, index_ninio_WPV], lat_min[i, index_ninio_WPV],
				  lon_min[i, index_ninio_WPV], val_max[i, index_ninio_WPV],
				  lat_max[i, index_ninio_WPV], lon_max[i, index_ninio_WPV]])
	var_normal_WPV = np.array([val_min[i, index_normal_WPV], lat_min[i, index_normal_WPV],
				  lon_min[i, index_normal_WPV], val_max[i, index_normal_WPV],
				  lat_max[i, index_normal_WPV], lon_max[i, index_normal_WPV]])
	var_ninia_WPV =  np.array([val_min[i, index_ninia_WPV], lat_min[i, index_ninia_WPV],
				  lon_min[i, index_ninia_WPV], val_max[i, index_ninia_WPV],
				  lat_max[i, index_ninia_WPV], lon_max[i, index_ninia_WPV]])
	var_ninio_SPV = np.array([val_min[i, index_ninio_SPV], lat_min[i, index_ninio_SPV],
				  lon_min[i, index_ninio_SPV], val_max[i, index_ninio_SPV],
				  lat_max[i, index_ninio_SPV], lon_max[i, index_ninio_SPV]])
	var_normal_SPV = np.array([val_min[i, index_normal_SPV], lat_min[i, index_normal_SPV],
				  lon_min[i, index_normal_SPV], val_max[i, index_normal_SPV],
				  lat_max[i, index_normal_SPV], lon_max[i, index_normal_SPV]])
	var_ninia_SPV =  np.array([val_min[i, index_ninia_SPV], lat_min[i, index_ninia_SPV],
				  lon_min[i, index_ninia_SPV], val_max[i, index_ninia_SPV],
				  lat_max[i, index_ninia_SPV], lon_max[i, index_ninia_SPV]])
	var_ninio_all = np.array([val_min[i, index_ninio_all], lat_min[i, index_ninio_all],
				  lon_min[i, index_ninio_all], val_max[i, index_ninio_all],
				  lat_max[i, index_ninio_all], lon_max[i, index_ninio_all]])
	var_normal_all = np.array([val_min[i, index_normal_all], lat_min[i, index_normal_all],
				  lon_min[i, index_normal_all], val_max[i, index_normal_all],
				  lat_max[i, index_normal_all], lon_max[i, index_normal_all]])
	var_ninia_all =  np.array([val_min[i, index_ninia_all], lat_min[i, index_ninia_all],
				  lon_min[i, index_ninia_all], val_max[i, index_ninia_all],
				  lat_max[i, index_ninia_all], lon_max[i, index_ninia_all]])
	np.savez(PATH_DATA_2 + 'ASL_coords' + seas[i] + '.npz', var1=var_ninio_WPV,
		 var2=var_normal_WPV, var3=var_ninia_WPV, var4=var_ninio_SPV, var5=var_normal_SPV,
		 var6=var_ninia_SPV, var7=var_ninio_all, var8=var_normal_all, var9=var_ninia_all)


