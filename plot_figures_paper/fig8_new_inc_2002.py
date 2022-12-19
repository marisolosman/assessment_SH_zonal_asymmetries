#composites of EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
from numpy import ma
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
#from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.path as mpath
from cartopy.util import add_cyclic_point
from matplotlib.transforms import offset_copy
from scipy import stats

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/fogt/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures_paper/'
FILE_HGT_S4 = 'hgt50_monthly_harmonics.nc4'
FILE_HGT_CLIMA = 'hgt50_monthly_climatology_harmonics.nc4'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
hgt = xr.open_dataset(PATH_DATA_2 + FILE_HGT_S4)
hgt_clim = xr.open_dataset(PATH_DATA_2 + FILE_HGT_CLIMA)
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

nn_ninio_all = sum(index_ninio_all.values)
nn_ninia_all = sum(index_ninia_all.values)
nn_ninio_WPV = sum(index_ninio_WPV)
nn_ninia_WPV = sum(index_ninia_WPV)
nn_ninio_SPV = sum(index_ninio_SPV)
nn_ninia_SPV = sum(index_ninia_SPV)
nn_index_normal = sum(index_normal)
nn_SPV_upper = sum(index_SPV_upper.values)
nn_SPV_lower = sum(index_SPV_lower.values)

month = ['Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
var_ninio_WPV = np.mean(hgt.Harmonic1.values[:, index_ninio_WPV, :, :], axis=1)
SS_ninio_WPV = np.var(hgt.Harmonic1.values[:, index_ninio_WPV, :, :], axis=1) / nn_ninio_WPV
var_ninia_WPV = np.mean(hgt.Harmonic1.values[:, index_ninia_WPV, :, :], axis=1)	
SS_ninia_WPV = np.var(hgt.Harmonic1.values[:, index_ninia_WPV, :, :], axis=1) / nn_ninia_WPV
var_ninio_SPV = np.mean(hgt.Harmonic1.values[:, index_ninio_SPV, :, :], axis=1)
SS_ninio_SPV = np.var(hgt.Harmonic1.values[:, index_ninio_SPV, :, :], axis=1) / nn_ninio_SPV
var_ninia_SPV = np.mean(hgt.Harmonic1.values[:, index_ninia_SPV, :, :], axis=1)	
SS_ninia_SPV = np.var(hgt.Harmonic1.values[:, index_ninia_SPV, :, :], axis=1) / nn_ninia_SPV
var_ninio_all = np.mean(hgt.Harmonic1.values[:, index_ninio_all.values, :, :], axis=1)
SS_ninio_all = np.var(hgt.Harmonic1.values[:, index_ninio_all.values, :, :], axis=1) / nn_ninio_all
var_normal = np.mean(hgt.Harmonic1.values[:, index_normal, :, :], axis=1)
SS_normal = np.var(hgt.Harmonic1.values[:, index_normal, :, :], axis=1) / nn_index_normal
var_SPV_all = np.mean(hgt.Harmonic1.values[:, index_SPV_lower.values, :, :], axis=1)
SS_SPV_all = np.var(hgt.Harmonic1.values[:, index_SPV_lower.values, :, :], axis=1) / nn_SPV_lower
var_WPV_all = np.mean(hgt.Harmonic1.values[:, index_SPV_upper.values, :, :], axis=1)
SS_WPV_all = np.var(hgt.Harmonic1.values[:, index_SPV_upper.values, :, :], axis=1) / nn_SPV_upper
var_ninia_all = np.mean(hgt.Harmonic1.values[:, index_ninia_all.values, :, :], axis=1)
SS_ninia_all = np.var(hgt.Harmonic1.values[:, index_ninia_all.values, :, :], axis=1) / nn_ninia_all


tit = 'Composites S4 Z* 50hPa Wave-1'
filename = FIG_PATH + 'z50_W1_composites_ENSO_SPoV_2.eps'

titulos = ['Ninio & Weak SPV', 'Ninia & Strong SPV',
	   'Ninio & Strong SPV', 'Ninia & Weak SPV']
harmonic = 1
lon = hgt.longitude.values
lat = hgt.latitude.values
proj = ccrs.Stereographic(central_longitude=-60, central_latitude=-90)
clevs2 = np.arange(-150, 200, 50)
clevs = np.arange(-60, 70, 10)
barra = plt.cm.RdBu_r
theta = np.linspace(0, 2 * np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
fig, ax = plt.subplots(nrows=4, ncols=4, figsize=(7, 8.5), subplot_kw={'projection': proj})
for i in range(4):
	ax[i, 0].set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
	ax[i, 0].set_boundary(circle, transform=ax[i, 0].transAxes)
	ax[i, 0].coastlines()
	variable = var_ninio_WPV[i + 1, :, :] - var_normal[i + 1, :, :]
	variable_std = variable / np.sqrt(SS_ninio_WPV[i + 1, :, :] + SS_normal[i + 1, :, :])
	cota = stats.t.ppf(1 - 0.025, nn_ninio_WPV + nn_index_normal - 2)
	M = np.ma.array(variable, mask= np.abs(variable_std) < cota)
	ax[i, 0].contourf(lon, lat, M, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	ax[i, 0].contour(lon, lat, hgt_clim.Harmonic1.values[i + 1, :, :], clevs2,
			 transform=ccrs.PlateCarree(), linewidths=0.4, colors='black')
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
#
	ax[i, 1].set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
	ax[i, 1].set_boundary(circle, transform=ax[i, 1].transAxes)
	ax[i, 1].coastlines()
	variable = var_ninia_SPV[i + 1, :, :] - var_normal[i + 1, :, :]
	variable_std = variable / np.sqrt(SS_ninia_SPV[i + 1, :, :] + SS_normal[i + 1, :, :])
	cota = stats.t.ppf(1 - 0.025, nn_ninia_SPV + nn_index_normal - 2)
	M = np.ma.array(variable, mask= np.abs(variable_std) < cota)
	ax[i, 1].contourf(lon, lat, M, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	ax[i, 1].contour(lon,lat, hgt_clim.Harmonic1.values[i + 1, :, :], clevs2,
			 transform=ccrs.PlateCarree(), linewidths=0.4, colors='black')
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	
	ax[i, 2].set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
	ax[i, 2].set_boundary(circle, transform=ax[i, 2].transAxes)
	ax[i, 2].coastlines()
	variable = var_ninio_SPV[i + 1, :, :] - var_normal[i + 1, :, :]
	variable_std = variable / np.sqrt(SS_ninio_SPV[i + 1, :, :] + SS_normal[i + 1, :, :])
	cota = stats.t.ppf(1 - 0.025, nn_ninio_SPV + nn_index_normal - 2)
	M = np.ma.array(variable, mask= np.abs(variable_std) < cota)
	ax[i, 2].contourf(lon, lat, M, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	ax[i, 2].contour(lon, lat, hgt_clim.Harmonic1.values[i + 1, :, :], clevs2,
			 transform=ccrs.PlateCarree(), linewidths=0.4, colors='black')
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))

	ax[i, 3].set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
	ax[i, 3].set_boundary(circle, transform=ax[i, 3].transAxes)
	ax[i, 3].coastlines()
	variable = var_ninia_WPV[i + 1, :, :] - var_normal[i + 1, :, :]
	variable_std = variable / np.sqrt(SS_ninia_WPV[i + 1, :, :] + SS_normal[i + 1, :, :])
	cota = stats.t.ppf(1 - 0.025, nn_ninia_WPV + nn_index_normal - 2)
	M = np.ma.array(variable, mask= np.abs(variable_std) < cota)
	im = ax[i, 3].contourf(lon, lat, M, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	ax[i, 3].contour(lon, lat, hgt_clim.Harmonic1.values[i + 1, :, :], clevs2, transform=ccrs.PlateCarree(),
			      linewidths=0.4, colors='black')
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
pad = 5
for bx, col in zip(ax[0], titulos):
	bx.annotate(col, xy=(0.5, 1), xytext=(0, pad),
		   xycoords='axes fraction', textcoords='offset points',
		   size='large', ha='center', va='baseline', fontsize=8)
	#bx.set_title(col)
for bx, row in zip(ax[:, 0], month):
	bx.annotate(row, xy=(0, 0.5), xytext=(-bx.yaxis.labelpad-pad, 0),
		   xycoords='axes fraction', rotation=90, textcoords='offset points',
		   size='large', ha='right', va='center', fontsize=10)
#	bx.set_ylabel(row, size='large')
plt.suptitle(tit, fontsize=12, x=0.55, y=0.9)
fig.tight_layout()
fig.subplots_adjust(left=0.15, bottom=0.08, top=0.84, hspace=0.1, wspace=0.1)
cbar_ax = fig.add_axes([0.32, 0.05, 0.5, 0.02])
fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
plt.savefig(filename, dpi=500, bbox_inches='tight', papertype='A4')
plt.clf()
plt.cla()
plt.close()


tit = 'Composites S4 Z* 50hPa Wave-1'
filename = FIG_PATH + 'z50_W1_composites_ENSO_SPoV_3.eps'

titulos = ['Weak-Strong SPV \n cond Ninio', 'Weak-Strong SPV \n cond Ninia',
	   'Ninio-Ninia \n cond Weak SPV', 'Ninio-Ninia \n cond Strong SPV']
harmonic = 1
lon = hgt.longitude.values
lat = hgt.latitude.values
proj = ccrs.Stereographic(central_longitude=-60, central_latitude=-90)
clevs2 = np.arange(-150, 200, 50)
clevs = np.arange(-60, 70, 10)
barra = plt.cm.RdBu_r
theta = np.linspace(0, 2 * np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
fig, ax = plt.subplots(nrows=4, ncols=4, figsize=(7, 8.5), subplot_kw={'projection': proj})
for i in range(4):
	ax[i, 0].set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
	ax[i, 0].set_boundary(circle, transform=ax[i, 0].transAxes)
	ax[i, 0].coastlines()
	variable = var_ninio_WPV[i + 1, :, :] - var_ninio_SPV[i + 1, :, :]
	variable_std = variable / np.sqrt(SS_ninio_WPV[i + 1, :, :] + SS_ninio_SPV[i + 1, :, :])
	cota = stats.t.ppf(1 - 0.025, nn_ninio_WPV + nn_ninio_SPV - 2)
	M = np.ma.array(variable, mask= np.abs(variable_std) < cota)
	ax[i, 0].contourf(lon, lat, M, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	ax[i, 0].contour(lon, lat, hgt_clim.Harmonic1.values[i + 1, :, :], clevs2,
			 transform=ccrs.PlateCarree(), linewidths=0.4, colors='black')
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
#
	ax[i, 1].set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
	ax[i, 1].set_boundary(circle, transform=ax[i, 1].transAxes)
	ax[i, 1].coastlines()
	variable = var_ninia_WPV[i + 1, :, :] - var_ninia_SPV[i + 1, :, :]
	variable_std = variable / np.sqrt(SS_ninia_WPV[i + 1, :, :] + SS_ninia_SPV[i + 1, :, :])
	cota = stats.t.ppf(1 - 0.025, nn_ninia_WPV + nn_ninia_SPV - 2)
	M = np.ma.array(variable, mask= np.abs(variable_std) < cota)
	ax[i, 1].contourf(lon, lat, M, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	ax[i, 1].contour(lon,lat, hgt_clim.Harmonic1.values[i + 1, :, :], clevs2,
			 transform=ccrs.PlateCarree(), linewidths=0.4, colors='black')
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	
	ax[i, 2].set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
	ax[i, 2].set_boundary(circle, transform=ax[i, 2].transAxes)
	ax[i, 2].coastlines()
	variable = var_ninio_WPV[i + 1, :, :] - var_ninia_WPV[i + 1, :, :]
	variable_std = variable / np.sqrt(SS_ninio_WPV[i + 1, :, :] + SS_ninia_WPV[i + 1, :, :])
	cota = stats.t.ppf(1 - 0.025, nn_ninio_WPV + nn_ninia_WPV - 2)
	M = np.ma.array(variable, mask= np.abs(variable_std) < cota)
	ax[i, 2].contourf(lon, lat, M, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	ax[i, 2].contour(lon, lat, hgt_clim.Harmonic1.values[i + 1, :, :], clevs2,
			 transform=ccrs.PlateCarree(), linewidths=0.4, colors='black')
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))

	ax[i, 3].set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
	ax[i, 3].set_boundary(circle, transform=ax[i, 3].transAxes)
	ax[i, 3].coastlines()
	variable = var_ninio_SPV[i + 1, :, :] - var_ninia_SPV[i + 1, :, :]
	variable_std = variable / np.sqrt(SS_ninio_SPV[i + 1, :, :] + SS_ninia_SPV[i + 1, :, :])
	cota = stats.t.ppf(1 - 0.025, nn_ninio_SPV + nn_ninia_SPV - 2)
	M = np.ma.array(variable, mask= np.abs(variable_std) < cota)
	im = ax[i, 3].contourf(lon, lat, M, clevs, transform=ccrs.PlateCarree(),
			 cmap=barra, extend='both', vmin=-60, vmax=60)
	ax[i, 3].contour(lon, lat, hgt_clim.Harmonic1.values[i + 1, :, :], clevs2, transform=ccrs.PlateCarree(),
			      linewidths=0.4, colors='black')
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
pad = 5
for bx, col in zip(ax[0], titulos):
	bx.annotate(col, xy=(0.5, 1), xytext=(0, pad),
		   xycoords='axes fraction', textcoords='offset points',
		   size='large', ha='center', va='baseline', fontsize=8)
	#bx.set_title(col)
for bx, row in zip(ax[:, 0], month):
	bx.annotate(row, xy=(0, 0.5), xytext=(-bx.yaxis.labelpad-pad, 0),
		   xycoords='axes fraction', rotation=90, textcoords='offset points',
		   size='large', ha='right', va='center', fontsize=10)
#	bx.set_ylabel(row, size='large')
plt.suptitle(tit, fontsize=12, x=0.55, y=0.9)
fig.tight_layout()
fig.subplots_adjust(left=0.15, bottom=0.08, top=0.84, hspace=0.1, wspace=0.1)
cbar_ax = fig.add_axes([0.32, 0.05, 0.5, 0.02])
fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
plt.savefig(filename, dpi=500, bbox_inches='tight', papertype='A4')
plt.clf()
plt.cla()
plt.close()

