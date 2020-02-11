import numpy as np
import datetime
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs	
import cartopy.feature 	
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

RUTA = '~/datos/data/'
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/model_assessment/'
#FILE_HGT_S4 = 'monthly_hgt200_aug_feb.nc'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
FILE_NINIO_ERAI = 'fogt/ninio34_erai_index.nc4'
FILE_PV_ERAI = 'fogt/SPV_index_erai.nc4'
ninio34_erai =  xr.open_dataset(RUTA + FILE_NINIO_ERAI)
PV_erai =  xr.open_dataset(RUTA + FILE_PV_ERAI)
#PV_erai.z.values = PV_erai.z.values / 10
#PV_erai.time.values = PV_erai.valid_time.values
#PV_erai = PV_erai.sel({'time':slice('1981-08-01', '2018-02-28')})
#select ninio, ninia and neutral
index_ninio_erai_upper = ninio34_erai.ninio34_index >= ninio34_erai.ninio34_index.quantile(0.75, dim='dim_0')
index_ninio_erai_lower = ninio34_erai.ninio34_index <= ninio34_erai.ninio34_index.quantile(0.25, dim='dim_0')
index_ninio_erai_normal = np.logical_and(ninio34_erai.ninio34_index < ninio34_erai.ninio34_index.quantile(0.75, dim='dim_0'), ninio34_erai.ninio34_index > ninio34_erai.ninio34_index.quantile(0.25, dim='dim_0'))
index_sninio_erai_upper = ninio34_erai.ninio34_index >= ninio34_erai.ninio34_index.quantile(0.90, dim='dim_0')
index_sninio_erai_lower = ninio34_erai.ninio34_index <= ninio34_erai.ninio34_index.quantile(0.10, dim='dim_0')
index_sninio_erai_normal = np.logical_and(ninio34_erai.ninio34_index < ninio34_erai.ninio34_index.quantile(0.90, dim='dim_0'), ninio34_erai.ninio34_index > ninio34_erai.ninio34_index.quantile(0.10, dim='dim_0'))

#idem s4
#abrir archivo de sst y PV del erai
ninio34_s4 =  xr.open_dataset(RUTA + FILE_NINIO_S4)
PV_s4 =  xr.open_dataset(RUTA + FILE_PV_S4)
#PV_s4 = PV_s4.stack(realiz=['year', 'number'])
#select ninio, ninia and neutral
index_ninio_s4_upper = ninio34_s4.ninio34_index >= ninio34_s4.ninio34_index.quantile(0.75, dim='dim_0')
index_ninio_s4_lower = ninio34_s4.ninio34_index <= ninio34_s4.ninio34_index.quantile(0.25, dim='dim_0')
index_ninio_s4_normal = np.logical_and(ninio34_s4.ninio34_index < ninio34_s4.ninio34_index.quantile(0.75, dim='dim_0'), ninio34_s4.ninio34_index > ninio34_s4.ninio34_index.quantile(0.25, dim='dim_0'))
index_sninio_s4_upper = ninio34_s4.ninio34_index >= ninio34_s4.ninio34_index.quantile(0.90, dim='dim_0')
index_sninio_s4_lower = ninio34_s4.ninio34_index <= ninio34_s4.ninio34_index.quantile(0.10, dim='dim_0')
index_sninio_s4_normal = np.logical_and(ninio34_s4.ninio34_index < ninio34_s4.ninio34_index.quantile(0.90, dim='dim_0'), ninio34_s4.ninio34_index > ninio34_s4.ninio34_index.quantile(0.10, dim='dim_0'))

season = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']
lmonth = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec']
#loop over seasons, selec data and performs boxplot
PV_erai_ninio = np.zeros([5, index_ninio_erai_upper.values.sum()])
PV_erai_sninio = np.zeros([5, index_sninio_erai_upper.values.sum()])
PV_erai_ninia = np.zeros([5, index_ninio_erai_lower.values.sum()])
PV_erai_sninia = np.zeros([5, index_sninio_erai_lower.values.sum()])
PV_s4_ninio = np.zeros([5, index_ninio_s4_upper.values.sum()])
PV_s4_sninio = np.zeros([5, index_sninio_s4_upper.values.sum()])
PV_s4_ninia = np.zeros([5, index_ninio_s4_lower.values.sum()])
PV_s4_sninia = np.zeros([5, index_sninio_s4_lower.values.sum()])

for i in np.arange(0, 5):
	aux = PV_erai['SPV_mon'].resample(time='QS-' + lmonth[i]).mean(dim='time',skipna=True)
	mes = datetime.datetime.strptime(lmonth[i], '%b').month
	aux = aux.sel(time= np.logical_and(aux['time.month'] == mes, aux['time.year']!=2002))
	PV_erai_ninio[i, :] = (aux.values[index_ninio_erai_upper] - aux.values.mean()) / aux.values.std()
	PV_erai_ninia[i, :] = (aux.values[index_ninio_erai_lower] - aux.values.mean()) / aux.values.std()
	PV_erai_sninio[i, :] = (aux.values[index_sninio_erai_upper] - aux.values.mean()) / aux.values.std()
	PV_erai_sninia[i, :] = (aux.values[index_sninio_erai_lower] - aux.values.mean()) / aux.values.std()
	PV_s4_smean = np.nanmean(PV_s4.SPV_index.values[i:i + 3, :], axis=0)	
	PV_s4_ninio[i, :] = (PV_s4_smean[index_ninio_s4_upper] - PV_s4_smean.mean()) / PV_s4_smean.std()
	PV_s4_ninia[i, :] = (PV_s4_smean[index_ninio_s4_lower] - PV_s4_smean.mean()) / PV_s4_smean.std()
	PV_s4_sninio[i, :] = (PV_s4_smean[index_sninio_s4_upper] - PV_s4_smean.mean()) / PV_s4_smean.std()
	PV_s4_sninia[i, :] = (PV_s4_smean[index_sninio_s4_lower] - PV_s4_smean.mean()) / PV_s4_smean.std()

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
# plot violin plot
axes[0].boxplot(PV_erai_ninio.T, positions=np.arange(1.20, 8.8,1.6))
axes[0].boxplot(PV_erai_sninio.T, positions=np.arange(1.8,9.3,1.6))

axes[0].set_title('Strat PoV index - ERAI')

# plot box plot
axes[1].violinplot(PV_s4_ninio.T,positions=np.arange(1.20, 8.8, 1.6),
		   showmeans=False,
		   showmedians=True)
axes[1].violinplot(PV_s4_sninio.T, positions=np.arange(1.8, 9.3, 1.6),
		   showmeans=False, showmedians=True)

axes[1].set_title('Strat PoV index - S4')
for ax in axes:
	# adding horizontal grid lines
	ax.yaxis.grid(True)
	#ax.set_xticks(range(5))
	#ax.set_xtickslabel(season)
	ax.set_xlim([0.5, 9])
	ax.set_ylim([-4.5, 4.5])
	ax.set_xlabel('Seasons')
	ax.set_ylabel('Polar Vortex')
# add x-tick labels
plt.setp(axes, xticks=np.arange(1.5, 9, 1.6), xticklabels=season)

# add x-tick labels
plt.savefig(FIG_PATH + 'boxplot_PV_ninio.png', res=400)


fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
# plot violin plot
b=axes[0].boxplot(PV_erai_ninio.T, positions=np.arange(1.20, 8.8,1.6))
axes[0].boxplot(PV_erai_sninio.T, positions=np.arange(1.8,9.3,1.6), patch_artist=True)
axes[0].set_title('Strat PoV index - ERAI')

# plot box plot
axes[1].boxplot(PV_s4_ninio.T,positions=np.arange(1.20, 8.8, 1.6))
axes[1].boxplot(PV_s4_sninio.T, positions=np.arange(1.8, 9.3, 1.6), patch_artist=True)

axes[1].set_title('Strat PoV index - S4')
for ax in axes:
	# adding horizontal grid lines
	ax.yaxis.grid(True)
	#ax.set_xticks(range(5))
	#ax.set_xtickslabel(season)
	ax.set_xlim([0.5, 9])
	ax.set_ylim([-4.5, 4.5])
	ax.set_xlabel('Seasons')
	ax.set_ylabel('Polar Vortex')
# add x-tick labels
plt.setp(axes, xticks=np.arange(1.5, 9, 1.6), xticklabels=season)

# add x-tick labels
plt.savefig(FIG_PATH + 'boxplot_PV_ninio_2.png', res=400)



