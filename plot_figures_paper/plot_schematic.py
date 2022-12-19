#plot schematic of ENSO and SPV
from __future__ import unicode_literals
import numpy as np
import xarray as xr
import numpy.ma as ma
import os
from numpy import ma
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
import matplotlib.path as mpath
from cartopy.util import add_cyclic_point
from matplotlib.transforms import offset_copy

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '/home/users/vg140344/datos/data/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures_paper/'
FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc4'

FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
hgt = xr.open_dataset(PATH_DATA_2 + FILE_HGT_S4)
ninio34 =  xr.open_dataset(PATH_DATA + FILE_NINIO_S4)
PV_index =  xr.open_dataset(PATH_DATA + FILE_PV_S4)
hgt = hgt - hgt.mean(dim='longitude')
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

month = ['Oct']

var_ninio_WPV = np.mean(hgt.z.values[:, index_ninio_WPV, :, :], axis=1)
SS_ninio_WPV = np.var(hgt.z.values[:, index_ninio_WPV, :, :], axis=1) / nn_ninio_WPV
var_ninia_WPV = np.mean(hgt.z.values[:, index_ninia_WPV, :, :], axis=1)
SS_ninia_WPV = np.var(hgt.z.values[:, index_ninia_WPV, :, :], axis=1) / nn_ninia_WPV
var_ninio_SPV = np.mean(hgt.z.values[:, index_ninio_SPV, :, :], axis=1)
SS_ninio_SPV = np.var(hgt.z.values[:, index_ninio_SPV, :, :], axis=1) / nn_ninio_SPV
var_ninia_SPV = np.mean(hgt.z.values[:, index_ninia_SPV, :, :], axis=1)
SS_ninia_SPV = np.var(hgt.z.values[:, index_ninia_SPV, :, :], axis=1) / nn_ninia_SPV
var_ninio_all = np.mean(hgt.z.values[:, index_ninio_all.values, :, :], axis=1)
SS_ninio_all = np.var(hgt.z.values[:, index_ninio_all.values, :, :], axis=1) / nn_ninio_all
var_normal = np.mean(hgt.z.values[:, index_normal, :, :], axis=1)
SS_normal = np.var(hgt.z.values[:, index_normal, :, :], axis=1) / nn_index_normal
var_SPV_all = np.mean(hgt.z.values[:, index_SPV_lower.values, :, :], axis=1)
SS_SPV_all = np.var(hgt.z.values[:, index_SPV_lower.values, :, :], axis=1) / nn_SPV_lower
var_WPV_all = np.mean(hgt.z.values[:, index_SPV_upper.values, :, :], axis=1)
SS_WPV_all = np.var(hgt.z.values[:, index_SPV_upper.values, :, :], axis=1) / nn_SPV_upper
var_ninia_all = np.mean(hgt.z.values[:, index_ninia_all.values, :, :], axis=1)
SS_ninia_all = np.var(hgt.z.values[:, index_ninia_all.values, :, :], axis=1) / nn_ninia_all

tit = 'Z* response to ENSO and SPV'
filename = FIG_PATH + 'z200_schematic.eps'

titulos = ['a) Niño - Niña cond SPV', 'b) Weak - Strong SPV cond ENSO',
	   'c) In-Phase Events', 'd) Out-of-Phase Events']
lon = hgt.longitude.values
lat = hgt.latitude.values
proj = ccrs.Stereographic(central_longitude=-60, central_latitude=-90)
barra = plt.cm.RdBu_r
theta = np.linspace(0, 2 * np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(7, 7.5), subplot_kw={'projection': proj})
for i in [2]:
	# ENSO
	ax[0, 0].set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
	ax[0, 0].set_boundary(circle, transform=ax[0, 0].transAxes)
	ax[0, 0].coastlines(linewidths=0.4, color='gray')
	variable = var_ninio_WPV[i, :, :] - var_ninia_WPV[i, :, :]
	clevs = [-16, 16]
	ax[0, 0].contour(lon, lat, variable, levels=[clevs[0]], transform=ccrs.PlateCarree(),
			 linewidths=1, linestyles='solid', colors='blue')
	ax[0, 0].contour(lon, lat, variable, levels=[clevs[1]], transform=ccrs.PlateCarree(),
			 linewidths=1, linestyles='solid', colors='red')
	ax[0,0].title.set_text(titulos[0])
	# SPV
	ax[0, 1].set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
	ax[0, 1].set_boundary(circle, transform=ax[0, 1].transAxes)
	ax[0, 1].coastlines(linewidths=0.4, color='gray')
	variable = var_ninia_WPV[i, :, :] - var_ninia_SPV[i, :, :]
	clevs = [-18, 18]
	ax[0, 1].contour(lon, lat, variable, levels=[clevs[0]], transform=ccrs.PlateCarree(),
			 linewidths=1, linestyles='dotted', colors='blue')
	ax[0, 1].contour(lon, lat, variable, levels=[clevs[1]], transform=ccrs.PlateCarree(),
			 linewidths=1, linestyles='dotted', colors='red')
	ax[0, 1].title.set_text(titulos[1])
	# in-phase
	ax[1, 0].set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
	ax[1, 0].set_boundary(circle, transform=ax[1, 0].transAxes)
	ax[1, 0].coastlines(linewidths=0.4, color='gray')
	variable = var_ninio_WPV[i, :, :] - var_ninia_WPV[i, :, :]
	constructive = np.logical_or(np.logical_and((var_ninio_WPV[i, :, :] - var_ninia_WPV[i, :, :])> 16,
			(var_ninia_WPV[i, :, :] - var_ninia_SPV[i, :, :]) > 18), 
			np.logical_and((var_ninio_WPV[i, :, :] - var_ninia_WPV[i, :, :])< -16,
			(var_ninia_WPV[i, :, :] - var_ninia_SPV[i, :, :]) < -18))

	destructive_mask = np.logical_or(np.logical_and((var_ninio_WPV[i, :, :] - var_ninia_WPV[i, :, :])> 16,
			(var_ninia_WPV[i, :, :] - var_ninia_SPV[i, :, :]) < -18), 
			np.logical_and((var_ninio_WPV[i, :, :] - var_ninia_WPV[i, :, :]) < -16,
			(var_ninia_WPV[i, :, :] - var_ninia_SPV[i, :, :]) > 18))
	constructive = ma.masked_array(variable, mask=np.logical_not(constructive))
	destructive = ma.masked_array(variable, mask=np.logical_not(destructive_mask))
	clevs = [-16, 16]
	ax[1, 0].contour(lon, lat, variable, levels=[clevs[1]], transform=ccrs.PlateCarree(),
			 linewidths=1, linestyles='solid', colors='red')
	ax[1, 0].contour(lon, lat, variable, levels=[clevs[0]], transform=ccrs.PlateCarree(),
			 linewidths=1, linestyles='solid', colors='blue')
	clevs = [-18, 18]
	variable = var_ninia_WPV[i, :, :] - var_ninia_SPV[i, :, :]
	ax[1, 0].contour(lon, lat, variable, levels=[clevs[0]], transform=ccrs.PlateCarree(),
			 linewidths=1, linestyles='dotted', colors='blue')
	ax[1, 0].contour(lon, lat, variable, levels=[clevs[1]], transform=ccrs.PlateCarree(),
			 linewidths=1, linestyles='dotted', colors='red')
	ax[1, 0].contourf(lon, lat, constructive, levels=[-100, -16], colors='lightblue',
				  transform=ccrs.PlateCarree(), alpha=0.3)
	ax[1, 0].contourf(lon, lat, constructive, levels=[16, 100], colors='lightsalmon', alpha=0.3,
				  transform=ccrs.PlateCarree())
	try:
		clevs = [-16, 16]
		ax[1, 0].contourf(lon, lat, destructive, levels=clevs, extend='both', colors='lightgray',
				 transform=ccrs.PlateCarree())
	except:
		pass
	ax[1, 0].title.set_text(titulos[2])
	ax[1, 1].set_extent([0, 359, -90, -30], crs=ccrs.PlateCarree())
	ax[1, 1].set_boundary(circle, transform=ax[1, 1].transAxes)
	ax[1, 1].coastlines(linewidths=0.4, color='gray')
	variable = var_ninio_WPV[i, :, :] - var_ninia_WPV[i, :, :]
	constructive = np.logical_or(np.logical_and((var_ninio_WPV[i, :, :] - var_ninia_WPV[i, :, :])>= 16,
			(-1 * var_ninia_WPV[i, :, :] + var_ninia_SPV[i, :, :]) >= 18), 
			np.logical_and((var_ninio_WPV[i, :, :] - var_ninia_WPV[i, :, :])<= -16,
			(-1 * var_ninia_WPV[i, :, :] + var_ninia_SPV[i, :, :]) <= -18))
	destructive = np.logical_or(np.logical_and((var_ninio_WPV[i, :, :] - var_ninia_WPV[i, :, :])>= 16,
			(-1 * var_ninia_WPV[i, :, :] + var_ninia_SPV[i, :, :]) <= -18), 
			np.logical_and((var_ninio_WPV[i, :, :] - var_ninia_WPV[i, :, :]) <= -16,
			(-1 * var_ninia_WPV[i, :, :] + var_ninia_SPV[i, :, :]) >= 18))
	constructive = ma.masked_array(variable, mask=np.logical_not(constructive))
	destructive = ma.masked_array(variable, mask=np.logical_not(destructive))
	ax[1, 1].contour(lon, lat, variable, levels=[clevs[0]], transform=ccrs.PlateCarree(),
			 linewidths=1, linestyles='solid', colors='blue')
	ax[1, 1].contour(lon, lat, variable, levels=[clevs[1]], transform=ccrs.PlateCarree(),
			 linewidths=1, linestyles='solid', colors='red')
	variable = var_ninia_WPV[i, :, :] - var_ninia_SPV[i, :, :]
	clevs = [-18, 18]
	ax[1, 1].contour(lon, lat, -1 * variable, levels=[clevs[0]], transform=ccrs.PlateCarree(),
			 linewidths=1, linestyles='dotted', colors='blue')
	ax[1, 1].contour(lon, lat, -1 * variable, levels=[clevs[1]], transform=ccrs.PlateCarree(),
			 linewidths=1, linestyles='dotted', colors='red')
	ax[1, 1].contourf(lon, lat, constructive, levels=[-100, -16], colors='lightblue',
				  transform=ccrs.PlateCarree(), alpha=0.3)
	try:
		ax[1, 1].contourf(lon, lat, constructive, levels=[16, 100], colors='lightsalmon', alpha=0.3,
				  transform=ccrs.PlateCarree())
	except:
		pass
	try:
		clevs = [-16, 16]
		ax[1, 1].contourf(lon, lat, destructive, levels=clevs, extend='both', colors='lightgray',
				  transform=ccrs.PlateCarree())
	except:
		pass
	ax[1, 1].title.set_text(titulos[3])

plt.suptitle(tit, fontsize=12, x=0.5, y=0.92)
fig.subplots_adjust(top=0.85, hspace=0.08)
plt.savefig(filename, dpi=500, bbox_inches='tight', papertype='A4')
plt.clf()
plt.cla()
plt.close()


