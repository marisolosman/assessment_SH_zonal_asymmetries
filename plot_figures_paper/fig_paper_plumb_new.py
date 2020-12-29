#composites of z50 and plumb fluxes for EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import plots_paper as plots
import plumb_flux
#================================================
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '~/datos/data/fogt/'
FIG_PATH = '/storage/silver/acrcc/vg140344/figures/strat_trop_zonal_asymmetries/quartile_new/'
FILE_HGT200_S4 = 'monthly_hgt200_aug_feb.nc4'
FILE_HGT50_S4 = 'monthly_hgt50_aug_feb.nc4'

FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
FILE_WINDS = 'fogt/monthly_winds200_aug_feb.nc4'
FILE_WINDS50 = 'fogt/monthly_winds50_aug_feb.nc4'

ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)

#search for years with weak PV
index_SPV_upper = PV_index.SPV_index >= PV_index.SPV_index.quantile(0.75, dim='dim_0', interpolation='linear')
#search for years with strong PV
index_SPV_lower = PV_index.SPV_index <= PV_index.SPV_index.quantile(0.25, dim='dim_0', interpolation='linear')
#search for years with strong PV
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

index_ninio_neutral = np.logical_and(index_ninio_all.values, index_SPV_normal)
index_ninia_neutral = np.logical_and(index_ninia_all.values, index_SPV_normal)

index_SPV_neutral = np.logical_and(index_SPV_lower, index_normal_all)
index_WPV_neutral = np.logical_and(index_SPV_upper, index_normal_all)

index_normal = np.logical_and(index_SPV_normal, index_normal_all)

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
if ~os.path.isfile('pf_oct_50_2.npz'):

	hgt50 = xr.open_dataset(PATH_DATA + FILE_HGT50_S4, chunks={'latitude':10})
	hgt50 = hgt50 - hgt50.mean(dim='longitude')
	hgt50 = hgt50.sel(latitude=slice(10, -90)).compute()
	winds = xr.open_dataset(PATH_DATA + FILE_WINDS50, chunks={'latitude':10})
	winds = winds.transpose('month', 'realiz', 'latitude', 'longitude')
	winds = winds.sel(latitude=slice(10, -90)).compute()
	winds_clm = winds.mean(dim='realiz')
	i = 2
	var50_normal = np.mean(hgt50.z.values[i, :, :, :], axis=0)

	var50_ninio_WPV = hgt50.z.values[i, index_ninio_WPV, :, :]
	pf50_ninio_WPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var50_ninio_WPV, var50_normal, hgt50.latitude.values, hgt50.longitude.values, 5000)

	var50_ninia_WPV = hgt50.z.values[i, index_ninia_WPV, :, :]
	pf50_ninia_WPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var50_ninia_WPV, var50_normal, hgt50.latitude.values, hgt50.longitude.values, 5000)

	var50_ninio_SPV = hgt50.z.values[i, index_ninio_SPV, :, :]
	pf50_ninio_SPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var50_ninio_SPV, var50_normal, hgt50.latitude.values, hgt50.longitude.values, 5000)

	var50_ninia_SPV = hgt50.z.values[i, index_ninia_SPV, :, :]
	pf50_ninia_SPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var50_ninia_SPV, var50_normal, hgt50.latitude.values, hgt50.longitude.values, 5000)

	var50_WPV_neutral = hgt50.z.values[i, index_WPV_neutral, :, :]
	pf50_WPV_neutral = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var50_WPV_neutral, var50_normal, hgt50.latitude.values, hgt50.longitude.values, 5000)

	var50_SPV_neutral = hgt50.z.values[i, index_SPV_neutral, :, :]
	pf50_SPV_neutral = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var50_SPV_neutral, var50_normal, hgt50.latitude.values, hgt50.longitude.values, 5000)

	var50_ninio = hgt50.z.values[i, index_ninio_neutral, :, :]
	pf50_ninio = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var50_ninio, var50_normal, hgt50.latitude.values, hgt50.longitude.values, 5000)

	var50_ninia = hgt50.z.values[i, index_ninia_neutral, :, :]
	pf50_ninia = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var50_ninia, var50_normal, hgt50.latitude.values, hgt50.longitude.values, 5000)

	var_l_50 = {'z1': pf50_WPV_neutral['divl'],
	      'px1': pf50_WPV_neutral['pxl'], 'py1': pf50_WPV_neutral['pyl'],
	      'z2': pf50_SPV_neutral['divl'],
	      'px2': pf50_SPV_neutral['pxl'], 'py2': pf50_SPV_neutral['pyl'],
	      'z3': pf50_ninio['divl'],
	      'px3': pf50_ninio['pxl'], 'py3': pf50_ninio['pyl'],
	      'z4': pf50_ninia['divl'],
	      'px4': pf50_ninia['pxl'], 'py4': pf50_ninia['pyl'],
	      'z5': pf50_ninio_WPV['divl'],
	      'px5': pf50_ninio_WPV['pxl'], 'py5': pf50_ninio_WPV['pyl'],
	      'z6': pf50_ninia_SPV['divl'],
	      'px6': pf50_ninia_SPV['pxl'], 'py6': pf50_ninia_SPV['pyl'],
	      'z7': pf50_ninio_SPV['divl'], 
	      'px7': pf50_ninio_SPV['pxl'], 'py7': pf50_ninio_SPV['pyl'],
	      'z8': pf50_ninia_WPV['divl'],
	      'px8': pf50_ninia_WPV['pxl'], 'py8': pf50_ninia_WPV['pyl']}

	var_nl_50 = {'z1': pf50_WPV_neutral['divnl'],
	      'px1': pf50_WPV_neutral['pxnl'], 'py1': pf50_WPV_neutral['pynl'],
	      'z2': pf50_SPV_neutral['divnl'],
	      'px2': pf50_SPV_neutral['pxnl'], 'py2': pf50_SPV_neutral['pynl'],
	      'z3': pf50_ninio['divnl'],
	      'px3': pf50_ninio['pxnl'], 'py3': pf50_ninio['pynl'],
	      'z4': pf50_ninia['divnl'],
	      'px4': pf50_ninia['pxnl'], 'py4': pf50_ninia['pynl'],
	      'z5': pf50_ninio_WPV['divnl'],
	      'px5': pf50_ninio_WPV['pxnl'], 'py5': pf50_ninio_WPV['pynl'],
	      'z6': pf50_ninia_SPV['divnl'],
	      'px6': pf50_ninia_SPV['pxnl'], 'py6': pf50_ninia_SPV['pynl'],
	      'z7': pf50_ninio_SPV['divnl'], 
	      'px7': pf50_ninio_SPV['pxnl'], 'py7': pf50_ninio_SPV['pynl'],
	      'z8': pf50_ninia_WPV['divnl'],
	      'px8': pf50_ninia_WPV['pxnl'], 'py8': pf50_ninia_WPV['pynl']}

	var_em_50 = {'z1': pf50_WPV_neutral['divem'],
	      'px1': pf50_WPV_neutral['pxem'], 'py1': pf50_WPV_neutral['pyem'],
	      'z2': pf50_SPV_neutral['divem'],
	      'px2': pf50_SPV_neutral['pxem'], 'py2': pf50_SPV_neutral['pyem'],
	      'z3': pf50_ninio['divem'],
	      'px3': pf50_ninio['pxem'], 'py3': pf50_ninio['pyem'],
	      'z4': pf50_ninia['divem'],
	      'px4': pf50_ninia['pxem'], 'py4': pf50_ninia['pyem'],
	      'z5': pf50_ninio_WPV['divem'],
	      'px5': pf50_ninio_WPV['pxem'], 'py5': pf50_ninio_WPV['pyem'],
	      'z6': pf50_ninia_SPV['divem'],
	      'px6': pf50_ninia_SPV['pxem'], 'py6': pf50_ninia_SPV['pyem'],
	      'z7': pf50_ninio_SPV['divem'], 
	      'px7': pf50_ninio_SPV['pxem'], 'py7': pf50_ninio_SPV['pyem'],
	      'z8': pf50_ninia_WPV['divem'],
	      'px8': pf50_ninia_WPV['pxem'], 'py8': pf50_ninia_WPV['pyem']}

	np.savez('pf_oct_50_2.npz', var_l_50, var_nl_50, var_em_50)

	del hgt50
else:
	f = np.load('pf_oct_50_2.npz')
	print(f.keys())
	var_l_50 = f['var1']
	var_nl_50 = f['var2']
	var_em_50 = f['var3']

#seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
hgt200 = xr.open_dataset(PATH_DATA + FILE_HGT200_S4, chunks={'latitude':10})
hgt200 = hgt200 - hgt200.mean(dim='longitude')
hgt200 = hgt200.sel(latitude=slice(10, -90)).compute()
winds = xr.open_dataset(PATH_DATA + FILE_WINDS, chunks={'latitude':10})
winds = winds.transpose('month', 'realiz', 'latitude', 'longitude')
winds = winds.sel(latitude=slice(10, -90)).compute()
winds_clm = winds.mean(dim='realiz')
i = 2
if ~os.path.isfile('pf_oct_200_2.npz'):
		
	var200_normal = np.mean(hgt200.z.values[i, :, :, :], axis=0)
	var200_ninio_WPV = hgt200.z.values[i, index_ninio_WPV, :, :]
	pf200_ninio_WPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_ninio_WPV, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)
	var200_ninia_WPV = hgt200.z.values[i, index_ninia_WPV, :, :]
	pf200_ninia_WPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_ninia_WPV, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)
	var200_ninio_SPV = hgt200.z.values[i, index_ninio_SPV, :, :]
	pf200_ninio_SPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_ninio_SPV, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)
	var200_ninia_SPV = hgt200.z.values[i, index_ninia_SPV, :, :]
	pf200_ninia_SPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_ninia_SPV, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)

	var200_WPV_neutral = hgt200.z.values[i, index_WPV_neutral, :, :]
	pf200_WPV_neutral = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_WPV_neutral, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)

	var200_SPV_neutral = hgt200.z.values[i, index_SPV_neutral, :, :]
	pf200_SPV_neutral = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_SPV_neutral, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)

	var200_ninio = hgt200.z.values[i, index_ninio_neutral, :, :]
	pf200_ninio = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_ninio, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)

	var200_ninia = hgt200.z.values[i, index_ninia_neutral, :, :]
	pf200_ninia = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_ninia, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)

	var_l_200 = {'z1': pf200_WPV_neutral['divl'],
	      'px1': pf200_WPV_neutral['pxl'], 'py1': pf200_WPV_neutral['pyl'],
	      'z2': pf200_SPV_neutral['divl'],
	      'px2': pf200_SPV_neutral['pxl'], 'py2': pf200_SPV_neutral['pyl'],
	      'z3': pf200_ninio['divl'],
	      'px3': pf200_ninio['pxl'], 'py3': pf200_ninio['pyl'],
	      'z4': pf200_ninia['divl'],
	      'px4': pf200_ninia['pxl'], 'py4': pf200_ninia['pyl'],
	      'z5': pf200_ninio_WPV['divl'],
	      'px5': pf200_ninio_WPV['pxl'], 'py5': pf200_ninio_WPV['pyl'],
	      'z6': pf200_ninia_SPV['divl'],
	      'px6': pf200_ninia_SPV['pxl'], 'py6': pf200_ninia_SPV['pyl'],
	      'z7': pf200_ninio_SPV['divl'], 
	      'px7': pf200_ninio_SPV['pxl'], 'py7': pf200_ninio_SPV['pyl'],
	      'z8': pf200_ninia_WPV['divl'],
	      'px8': pf200_ninia_WPV['pxl'], 'py8': pf200_ninia_WPV['pyl']}

	var_nl_200 = {'z1': pf200_WPV_neutral['divnl'],
	      'px1': pf200_WPV_neutral['pxnl'], 'py1': pf200_WPV_neutral['pynl'],
	      'z2': pf200_SPV_neutral['divnl'],
	      'px2': pf200_SPV_neutral['pxnl'], 'py2': pf200_SPV_neutral['pynl'],
	      'z3': pf200_ninio['divnl'],
	      'px3': pf200_ninio['pxnl'], 'py3': pf200_ninio['pynl'],
	      'z4': pf200_ninia['divnl'],
	      'px4': pf200_ninia['pxnl'], 'py4': pf200_ninia['pynl'],
	      'z5': pf200_ninio_WPV['divnl'],
	      'px5': pf200_ninio_WPV['pxnl'], 'py5': pf200_ninio_WPV['pynl'],
	      'z6': pf200_ninia_SPV['divnl'],
	      'px6': pf200_ninia_SPV['pxnl'], 'py6': pf200_ninia_SPV['pynl'],
	      'z7': pf200_ninio_SPV['divnl'], 
	      'px7': pf200_ninio_SPV['pxnl'], 'py7': pf200_ninio_SPV['pynl'],
	      'z8': pf200_ninia_WPV['divnl'],
	      'px8': pf200_ninia_WPV['pxnl'], 'py8': pf200_ninia_WPV['pynl']}

	var_em_200 = {'z1': pf200_WPV_neutral['divem'],
	      'px1': pf200_WPV_neutral['pxem'], 'py1': pf200_WPV_neutral['pyem'],
	      'z2': pf200_SPV_neutral['divem'],
	      'px2': pf200_SPV_neutral['pxem'], 'py2': pf200_SPV_neutral['pyem'],
	      'z3': pf200_ninio['divem'],
	      'px3': pf200_ninio['pxem'], 'py3': pf200_ninio['pyem'],
	      'z4': pf200_ninia['divem'],
	      'px4': pf200_ninia['pxem'], 'py4': pf200_ninia['pyem'],
	      'z5': pf200_ninio_WPV['divem'],
	      'px5': pf200_ninio_WPV['pxem'], 'py5': pf200_ninio_WPV['pyem'],
	      'z6': pf200_ninia_SPV['divem'],
	      'px6': pf200_ninia_SPV['pxem'], 'py6': pf200_ninia_SPV['pyem'],
	      'z7': pf200_ninio_SPV['divem'], 
	      'px7': pf200_ninio_SPV['pxem'], 'py7': pf200_ninio_SPV['pyem'],
	      'z8': pf200_ninia_WPV['divem'],
	      'px8': pf200_ninia_WPV['pxem'], 'py8': pf200_ninia_WPV['pyem']}

	np.savez('pf_oct_200_2.npz', var_l_200, var_nl_200, var_em_200)
else:
	f = np.load('pf_oct_200_2.npz')
	print(f.keys())
	var_l_200 = f['var1']
	var_nl_200 = f['var2']
	var_em_200 = f['var3']


tit = 'Composites EM Plumb fluxes and its divergence - ' + month[i]
filename = FIG_PATH + 'fig_pf_composites_oct_em_new8.eps'
plots.PlotCompositesdivPF(var_em_50, var_em_200, hgt200.latitude, hgt200.longitude, tit, filename)
tit = 'Composites EM linear term \n Plumb fluxes and its divergence - ' + month[i]
filename = FIG_PATH + 'fig_pf_composites_oct_lt_new8.eps'
plots.PlotCompositesdivPF(var_l_50, var_l_200, hgt200.latitude, hgt200.longitude, tit, filename)
tit = 'Composites EM non linear term \n Plumb fluxes and its divergence - ' + month[i]
filename = FIG_PATH + 'fig_pf_composites_oct_nlt_new8.eps'
#plots.PlotCompositesdivPF(var_nl_50, var_nl_200, hgt200.latitude, hgt200.longitude, tit, filename)
filename =  FIG_PATH + 'fig_pf_composites_oct_nlt_new8.eps'
plots.PlotCompositesdivPF(var_nl_50, var_nl_200, hgt200.latitude, hgt200.longitude, tit, filename)

