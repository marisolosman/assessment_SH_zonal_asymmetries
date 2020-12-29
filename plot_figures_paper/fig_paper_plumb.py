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

index_normal = np.logical_and(index_SPV_normal, index_normal_all)

month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
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
var50_ninio_all = hgt50.z.values[i, index_ninio_all, :, :]
pf50_ninio_all = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var50_ninio_all, var50_normal, hgt50.latitude.values, hgt50.longitude.values, 5000)
var50_ninia_all = hgt50.z.values[i, index_ninia_all, :, :]
pf50_ninia_all = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var50_ninia_all, var50_normal, hgt50.latitude.values, hgt50.longitude.values, 5000)
var50_SPV_all = hgt50.z.values[i, index_SPV_lower, :, :]
pf50_SPV_all = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var50_SPV_all, var50_normal, hgt50.latitude.values, hgt50.longitude.values, 5000)
var50_WPV_all = hgt50.z.values[i, index_SPV_upper, :, :]
pf50_WPV_all = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var50_WPV_all, var50_normal, hgt50.latitude.values, hgt50.longitude.values, 5000)

var_l_50 = {'z1': pf50_WPV_all['divl'],
      'px1': pf50_WPV_all['pxl'], 'py1': pf50_WPV_all['pyl'],
      'z2': pf50_SPV_all['divl'],
      'px2': pf50_SPV_all['pxl'], 'py2': pf50_SPV_all['pyl'],
      'z3': pf50_ninio_all['divl'],
      'px3': pf50_ninio_all['pxl'], 'py3': pf50_ninio_all['pyl'],
      'z4': pf50_ninia_all['divl'],
      'px4': pf50_ninia_all['pxl'], 'py4': pf50_ninia_all['pyl'],
      'z5': pf50_ninio_WPV['divl'],
      'px5': pf50_ninio_WPV['pxl'], 'py5': pf50_ninio_WPV['pyl'],
      'z6': pf50_ninia_SPV['divl'],
      'px6': pf50_ninia_SPV['pxl'], 'py6': pf50_ninia_SPV['pyl'],
      'z7': pf50_ninio_SPV['divl'], 
      'px7': pf50_ninio_SPV['pxl'], 'py7': pf50_ninio_SPV['pyl'],
      'z8': pf50_ninia_WPV['divl'],
      'px8': pf50_ninia_WPV['pxl'], 'py8': pf50_ninia_WPV['pyl']}
var_nl_50 = {'z1': pf50_WPV_all['divnl'],
      'px1': pf50_WPV_all['pxnl'], 'py1': pf50_WPV_all['pynl'],
      'z2': pf50_SPV_all['divnl'],
      'px2': pf50_SPV_all['pxnl'], 'py2': pf50_SPV_all['pynl'],
      'z3': pf50_ninio_all['divnl'],
      'px3': pf50_ninio_all['pxnl'], 'py3': pf50_ninio_all['pynl'],
      'z4': pf50_ninia_all['divnl'],
      'px4': pf50_ninia_all['pxnl'], 'py4': pf50_ninia_all['pynl'],
      'z5': pf50_ninio_WPV['divnl'],
      'px5': pf50_ninio_WPV['pxnl'], 'py5': pf50_ninio_WPV['pynl'],
      'z6': pf50_ninia_SPV['divnl'],
      'px6': pf50_ninia_SPV['pxnl'], 'py6': pf50_ninia_SPV['pynl'],
      'z7': pf50_ninio_SPV['divnl'], 
      'px7': pf50_ninio_SPV['pxnl'], 'py7': pf50_ninio_SPV['pynl'],
      'z8': pf50_ninia_WPV['divnl'],
      'px8': pf50_ninia_WPV['pxnl'], 'py8': pf50_ninia_WPV['pynl']}
var_em_50 = {'z1': pf50_WPV_all['divem'],
      'px1': pf50_WPV_all['pxem'], 'py1': pf50_WPV_all['pyem'],
      'z2': pf50_SPV_all['divem'],
      'px2': pf50_SPV_all['pxem'], 'py2': pf50_SPV_all['pyem'],
      'z3': pf50_ninio_all['divem'],
      'px3': pf50_ninio_all['pxem'], 'py3': pf50_ninio_all['pyem'],
      'z4': pf50_ninia_all['divem'],
      'px4': pf50_ninia_all['pxem'], 'py4': pf50_ninia_all['pyem'],
      'z5': pf50_ninio_WPV['divem'],
      'px5': pf50_ninio_WPV['pxem'], 'py5': pf50_ninio_WPV['pyem'],
      'z6': pf50_ninia_SPV['divem'],
      'px6': pf50_ninia_SPV['pxem'], 'py6': pf50_ninia_SPV['pyem'],
      'z7': pf50_ninio_SPV['divem'], 
      'px7': pf50_ninio_SPV['pxem'], 'py7': pf50_ninio_SPV['pyem'],
      'z8': pf50_ninia_WPV['divem'],
      'px8': pf50_ninia_WPV['pxem'], 'py8': pf50_ninia_WPV['pyem']}

del hgt50

#seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
hgt200 = xr.open_dataset(PATH_DATA + FILE_HGT200_S4, chunks={'latitude':10})
hgt200 = hgt200 - hgt200.mean(dim='longitude')
hgt200 = hgt200.sel(latitude=slice(10, -90)).compute()
winds = xr.open_dataset(PATH_DATA + FILE_WINDS, chunks={'latitude':10})
winds = winds.transpose('month', 'realiz', 'latitude', 'longitude')
winds = winds.sel(latitude=slice(10, -90)).compute()
winds_clm = winds.mean(dim='realiz')
i = 2

var200_normal = np.mean(hgt200.z.values[i, :, :, :], axis=0)
var200_ninio_WPV = hgt200.z.values[i, index_ninio_WPV, :, :]
pf200_ninio_WPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_ninio_WPV, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)
var200_ninia_WPV = hgt200.z.values[i, index_ninia_WPV, :, :]
pf200_ninia_WPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_ninia_WPV, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)
var200_ninio_SPV = hgt200.z.values[i, index_ninio_SPV, :, :]
pf200_ninio_SPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_ninio_SPV, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)
var200_ninia_SPV = hgt200.z.values[i, index_ninia_SPV, :, :]
pf200_ninia_SPV = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_ninia_SPV, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)
var200_ninio_all = hgt200.z.values[i, index_ninio_all, :, :]
pf200_ninio_all = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_ninio_all, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)
var200_ninia_all = hgt200.z.values[i, index_ninia_all, :, :]
pf200_ninia_all = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_ninia_all, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)
var200_SPV_all = hgt200.z.values[i, index_SPV_lower, :, :]
pf200_SPV_all = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_SPV_all, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)
var200_WPV_all = hgt200.z.values[i, index_SPV_upper, :, :]
pf200_WPV_all = plumb_flux.ComputePlumbFluxesDecomp(winds_clm.u.values[i, :, :], winds_clm.v.values[i, :, :], var200_WPV_all, var200_normal, hgt200.latitude.values, hgt200.longitude.values, 20000)
var_l_200 ={'z1': pf200_WPV_all['divl'],
      'px1': pf200_WPV_all['pxl'], 'py1': pf200_WPV_all['pyl'],
      'z2': pf200_SPV_all['divl'],
      'px2': pf200_SPV_all['pxl'], 'py2': pf200_SPV_all['pyl'],
      'z3': pf200_ninio_all['divl'],
      'px3': pf200_ninio_all['pxl'], 'py3': pf200_ninio_all['pyl'],
      'z4': pf200_ninia_all['divl'],
      'px4': pf200_ninia_all['pxl'], 'py4': pf200_ninia_all['pyl'],
      'z5': pf200_ninio_WPV['divl'],
      'px5': pf200_ninio_WPV['pxl'], 'py5': pf200_ninio_WPV['pyl'],
      'z6': pf200_ninia_SPV['divl'],
      'px6': pf200_ninia_SPV['pxl'], 'py6': pf200_ninia_SPV['pyl'],
      'z7': pf200_ninio_SPV['divl'], 
      'px7': pf200_ninio_SPV['pxl'], 'py7': pf200_ninio_SPV['pyl'],
      'z8': pf200_ninia_WPV['divl'],
      'px8': pf200_ninia_WPV['pxl'], 'py8': pf200_ninia_WPV['pyl']}
var_nl_200 ={'z1': pf200_WPV_all['divnl'],
      'px1': pf200_WPV_all['pxnl'], 'py1': pf200_WPV_all['pynl'],
      'z2': pf200_SPV_all['divnl'],
      'px2': pf200_SPV_all['pxnl'], 'py2': pf200_SPV_all['pynl'],
      'z3': pf200_ninio_all['divnl'],
      'px3': pf200_ninio_all['pxnl'], 'py3': pf200_ninio_all['pynl'],
      'z4': pf200_ninia_all['divnl'],
      'px4': pf200_ninia_all['pxnl'], 'py4': pf200_ninia_all['pynl'],
      'z5': pf200_ninio_WPV['divnl'],
      'px5': pf200_ninio_WPV['pxnl'], 'py5': pf200_ninio_WPV['pynl'],
      'z6': pf200_ninia_SPV['divnl'],
      'px6': pf200_ninia_SPV['pxnl'], 'py6': pf200_ninia_SPV['pynl'],
      'z7': pf200_ninio_SPV['divnl'], 
      'px7': pf200_ninio_SPV['pxnl'], 'py7': pf200_ninio_SPV['pynl'],
      'z8': pf200_ninia_WPV['divnl'],
      'px8': pf200_ninia_WPV['pxnl'], 'py8': pf200_ninia_WPV['pynl']}
var_em_200 ={'z1': pf200_WPV_all['divem'],
      'px1': pf200_WPV_all['pxem'], 'py1': pf200_WPV_all['pyem'],
      'z2': pf200_SPV_all['divem'],
      'px2': pf200_SPV_all['pxem'], 'py2': pf200_SPV_all['pyem'],
      'z3': pf200_ninio_all['divem'],
      'px3': pf200_ninio_all['pxem'], 'py3': pf200_ninio_all['pyem'],
      'z4': pf200_ninia_all['divem'],
      'px4': pf200_ninia_all['pxem'], 'py4': pf200_ninia_all['pyem'],
      'z5': pf200_ninio_WPV['divem'],
      'px5': pf200_ninio_WPV['pxem'], 'py5': pf200_ninio_WPV['pyem'],
      'z6': pf200_ninia_SPV['divem'],
      'px6': pf200_ninia_SPV['pxem'], 'py6': pf200_ninia_SPV['pyem'],
      'z7': pf200_ninio_SPV['divem'], 
      'px7': pf200_ninio_SPV['pxem'], 'py7': pf200_ninio_SPV['pyem'],
      'z8': pf200_ninia_WPV['divem'],
      'px8': pf200_ninia_WPV['pxem'], 'py8': pf200_ninia_WPV['pyem']}
tit = 'Composites EM Plumb fluxes and its divergence - ' + month[i]
filename = FIG_PATH + 'fig_pf_composites_oct_em.eps'
plots.PlotCompositesdivPF(var_em_50, var_em_200, hgt200.latitude, hgt200.longitude, tit, filename)
tit = 'Composites EM linear term \n Plumb fluxes and its divergence - ' + month[i]
filename = FIG_PATH + 'fig_pf_composites_oct_lt.eps'
plots.PlotCompositesdivPF(var_l_50, var_l_200, hgt200.latitude, hgt200.longitude, tit, filename)
tit = 'Composites EM non linear term \n Plumb fluxes and its divergence - ' + month[i]
filename = FIG_PATH + 'fig_pf_composites_oct_nlt.eps'
#plots.PlotCompositesdivPF(var_nl_50, var_nl_200, hgt200.latitude, hgt200.longitude, tit, filename)
filename =  FIG_PATH + 'fig_pf_composites_oct_nlt2.eps'
plots.PlotCompositesdivPF(var_nl_50, var_nl_200, hgt200.latitude, hgt200.longitude, tit, filename)

