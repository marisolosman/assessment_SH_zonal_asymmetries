#composites of z50 and plumb fluxes for EN events conditioned on PV strength
import numpy as np
import xarray as xr
import os
import plots_paper as plots
import matplotlib.pyplot as plt
#================================================
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
PATH_DATA = '~/datos/data/'
PATH_DATA_2 = '~/datos/data/fogt/'
FIG_PATH = '/storage/silver/acrcc/vg140344/figures/strat_trop_zonal_asymmetries/quartile_new/'
FILE_HGT = 'HGT_S4_Mlevels_60S.nc4'

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
hgt = xr.open_dataset(PATH_DATA + FILE_HGT)
#hgt = hgt.groupby('realiz.month')
hgt = hgt - hgt.mean(dim='longitude')
hgt = hgt.sel(realiz=hgt['realiz.month']==12)
#hgt = hgt.transpose('realiz', 'number', 'isobaricInhPa', 'longitude')
hgt = hgt.stack(experiment=['realiz', 'number'])
#hgt = hgt.transpose(['experiment', 'isobaricInhPa', 'longitude'])
hgt_normal = np.mean(hgt.z.values, axis=2)
var_ninio_WPV = np.mean(hgt.z.values[:, :, index_ninio_WPV], axis=2)
var_ninia_WPV = np.mean(hgt.z.values[:, :, index_ninia_WPV], axis=2)
var_ninio_SPV = np.mean(hgt.z.values[:, :, index_ninio_SPV], axis=2)
var_ninia_SPV = np.mean(hgt.z.values[:, :, index_ninia_SPV], axis=2)
var_ninio_all = np.mean(hgt.z.values[:, :, index_ninio_all], axis=2)
var_ninia_all = np.mean(hgt.z.values[:, :, index_ninia_all], axis=2)
var_SPV_all = np.mean(hgt.z.values[:, :, index_SPV_lower], axis=2)
var_WPV_all = np.mean(hgt.z.values[:, :, index_SPV_upper], axis=2)
var = {'z1': var_WPV_all, 'z2': var_SPV_all,
	'z3': var_ninio_all,  'z4': var_ninia_all,
	'z5': var_ninio_WPV, 'z6': var_ninia_SPV,
	'z7': var_ninio_SPV, 'z8': var_ninia_WPV}

titulo = 'Composites z* 60S Nov'
filename = FIG_PATH + 'fig_z60S_composites_nov.eps'
var_key = [*var]
tit = ['Weak SPoV (All ENSO)', 'Strong SpoV (All ENSO)',
	'Ninio (All PoV)', 'Ninia (All PoV)',
	'Ninio & Weak PoV', 'Ninia & Strong PoV',
	'Ninio & Strong PoV', 'Ninia & Weak PoV']
clevs = np.arange(-120, 140, 20)
clevs2 = np.arange(-300, 400, 100)
barra = plt.cm.RdBu_r
fig = plt.figure(1, (15, 15), 300)
for i in range(8):
	ax = plt.subplot(4, 2, i + 1)
	im = ax.contourf(hgt.longitude.values, hgt.isobaricInhPa, (var[var_key[i]] - hgt_normal) / 10,
		       clevs, cmap=barra, extend='both', vmin=clevs[0], vmax=clevs[-1])
	ax.set_yscale('log')
	ax.invert_yaxis()
	ax.set_yticks([1000, 100, 10])
	ax.set_yticklabels(['1000', '100', '10'])
	barra.set_under(barra(0))
	barra.set_over(barra(barra.N-1))
	ax.contour(hgt.longitude.values, hgt.isobaricInhPa, hgt_normal / 10, clevs2,
		   colors='black')
	plt.title(tit[i])
plt.suptitle(titulo, fontsize=12, x=0.47, y=0.9)
fig.subplots_adjust(right=0.8)
fig.subplots_adjust(bottom=0.2, top=0.85, hspace=0.25, wspace=0.15)
cbar_ax = fig.add_axes([0.34, 0.1, 0.25, 0.05])
#plt.quiverkey(im2, 0.5, -0.8, 1, "1 m/s", coordinates='axes', color='k')
fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
plt.savefig(filename, dpi=300, bbox_inches='tight', papertype='A4')
plt.clf()
plt.cla()
plt.close()

