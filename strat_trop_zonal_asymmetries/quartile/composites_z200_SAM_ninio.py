#composites of PV events conditioned on ENSO phase
import numpy as np
import xarray as xr
import os
import plots

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/strat_trop_zonal_asymmetries/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_SAM_S4 = 'fogt/SAM_monthly_index_s4.nc4'
hgt = xr.open_dataset(PATH_DATA + FILE_HGT_S4)
hgt = hgt - hgt.mean(dim='longitude')
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
SAM_index =  xr.open_dataset(PATH_DATA + FILE_SAM_S4)

#search for EN years 
index_EN = ninio34.ninio34_index >= ninio34.ninio34_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for LN years
index_LN = ninio34.ninio34_index <= ninio34.ninio34_index.quantile(0.25, dim='dim_0', interpolation='linear')

# compute SAM composites conditioned on ENSO phase
SAM_index_EN = SAM_index.SAM_index.sel(realiz = index_EN.values)
hgt_EN = hgt.sel(realiz = index_EN.values)
SAM_index_LN = SAM_index.SAM_index.sel(realiz = index_LN.values)
hgt_LN = hgt.sel(realiz = index_LN.values)

#PV during all phases
index_SSAM_all = SAM_index.SAM_index <= SAM_index.SAM_index.quantile(0.75, dim='realiz', interpolation='linear')
index_WSAM_all = SAM_index.SAM_index >= SAM_index.SAM_index.quantile(0.25, dim='realiz', interpolation='linear')
#PV  during EN
index_SSAM_EN = SAM_index_EN <= SAM_index_EN.quantile(0.75, dim='realiz', interpolation='linear')
index_WSAM_EN = SAM_index_EN >= SAM_index_EN.quantile(0.25, dim='realiz', interpolation='linear')
#PV  during LN
index_SSAM_LN = SAM_index_LN <= SAM_index_LN.quantile(0.75, dim='realiz', interpolation='linear')
index_WSAM_LN = SAM_index_LN >= SAM_index_LN.quantile(0.25, dim='realiz', interpolation='linear')

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

for i in np.arange(0, 7):
	var_WSAM_EN = np.mean(hgt_EN.z.values[i, index_WSAM_EN.values[i, :], :, :], axis=0)
	var_WSAM_LN = np.mean(hgt_LN.z.values[i, index_WSAM_LN.values[i, :], :, :], axis=0)	
	var_SSAM_EN = np.mean(hgt_EN.z.values[i, index_SSAM_EN.values[i, :], :, :], axis=0)
	var_SSAM_LN = np.mean(hgt_LN.z.values[i, index_SSAM_LN.values[i, :], :, :], axis=0)	
	var_WSAM_all = np.mean(hgt.z.values[i, index_WSAM_all.values[i, :], :, :], axis=0)
	var_SSAM_all = np.mean(hgt.z.values[i, index_SSAM_all.values[i, :], :, :], axis=0)
	tit = 'Composites S4 Z* 200hPa Strong SAM - Weak SAM Conditioned - ENSO - ' + month[i]
	filename = FIG_PATH + 'z200_composites_SAM_' + month[i] +'_ENSO_q.png'
	plots.PlotPoVCompositesDiffENSO(var_SSAM_all-var_WSAM_all, var_SSAM_EN - var_WSAM_EN,
				    var_SSAM_LN - var_WSAM_LN, hgt.latitude, hgt.longitude,
				    tit, filename)
FILE_SAM_S4 = 'fogt/SAM_index_s4.nc4'
SAM_index =  xr.open_dataset(PATH_DATA + FILE_SAM_S4)
# compute SAM composites conditioned on ENSO phase
SAM_index_EN = SAM_index.SAM_index.sel(realiz = index_EN.values)
hgt_EN = hgt.sel(realiz = index_EN.values)
SAM_index_LN = SAM_index.SAM_index.sel(realiz = index_LN.values)
hgt_LN = hgt.sel(realiz = index_LN.values)

#PV during all phases
index_SSAM_all = SAM_index.SAM_index <= SAM_index.SAM_index.quantile(0.75, dim='realiz', interpolation='linear')
index_WSAM_all = SAM_index.SAM_index >= SAM_index.SAM_index.quantile(0.25, dim='realiz', interpolation='linear')
#PV  during EN
index_SSAM_EN = SAM_index_EN <= SAM_index_EN.quantile(0.75, dim='realiz', interpolation='linear')
index_WSAM_EN = SAM_index_EN >= SAM_index_EN.quantile(0.25, dim='realiz', interpolation='linear')
#PV  during LN
index_SSAM_LN = SAM_index_LN <= SAM_index_LN.quantile(0.75, dim='realiz', interpolation='linear')
index_WSAM_LN = SAM_index_LN >= SAM_index_LN.quantile(0.25, dim='realiz', interpolation='linear')



for i in np.arange(0, 5):
	hgt_s = hgt.isel(month=range(i, i+3)).mean(dim='month')
	hgt_s_EN = hgt_s.sel(realiz=index_EN.values)
	hgt_s_LN = hgt_s.sel(realiz=index_LN.values)
	var_WSAM_EN = np.mean(hgt_s_EN.z.values[index_WSAM_EN.values[i, :], :, :], axis=0)
	var_WSAM_LN = np.mean(hgt_s_LN.z.values[index_WSAM_LN.values[i, :], :, :], axis=0)	
	var_SSAM_EN = np.mean(hgt_s_EN.z.values[index_SSAM_EN.values[i, :], :, :], axis=0)
	var_SSAM_LN = np.mean(hgt_s_LN.z.values[index_SSAM_LN.values[i, :], :, :], axis=0)	
	var_WSAM_all = np.mean(hgt_s.z.values[index_WSAM_all.values[i, :], :, :], axis=0)
	var_SSAM_all = np.mean(hgt_s.z.values[index_SSAM_all.values[i, :], :, :], axis=0)
	tit = 'Composites S4 Z* 200hPa Strong SAM - Weak SAM Conditioned - ENSO - ' + seas[i]
	filename = FIG_PATH + 'z200_composites_SAM_' + seas[i] +'_ENSO_q.png'
	plots.PlotPoVCompositesDiffENSO(var_SSAM_all-var_WSAM_all, var_SSAM_EN - var_WSAM_EN,
				    var_SSAM_LN - var_WSAM_LN, hgt.latitude, hgt.longitude,
				    tit, filename)


