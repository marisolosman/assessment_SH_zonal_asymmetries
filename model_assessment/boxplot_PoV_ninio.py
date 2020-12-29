# compute boxplot of spv for ninios
# compute regression between spv and enso
import numpy as np
import datetime
import pandas as pd
import xarray as xr
from eofs.standard import Eof
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import matplotlib.pyplot as plt
import cartopy.crs as ccrs	
import cartopy.feature 
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

RUTA = '~/datos/data/' #este spanglish no te lo robo amiga
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures/model_assessment/'
FILE_NINIO_S4 = 'fogt/ninio34_monthly.nc4'
FILE_PV_S4 = 'fogt/SPV_index.nc4'
FILE_NINIO_ERAI = 'fogt/ninio34_erai_index.nc4'
FILE_PV_ERAI = 'fogt/SPV_index_erai.nc4'
ninio34_erai =  xr.open_dataset(RUTA + FILE_NINIO_ERAI)
#select ninio, ninia and neutral
index_ninio_erai_upper = ninio34_erai.ninio34_index >= ninio34_erai.ninio34_index.quantile(0.75, dim='dim_0')
index_ninio_erai_lower = ninio34_erai.ninio34_index <= ninio34_erai.ninio34_index.quantile(0.25, dim='dim_0')
index_ninio_erai_normal = np.logical_and(ninio34_erai.ninio34_index < ninio34_erai.ninio34_index.quantile(0.75, dim='dim_0'), ninio34_erai.ninio34_index > ninio34_erai.ninio34_index.quantile(0.25, dim='dim_0'))
index_sninio_erai_upper = ninio34_erai.ninio34_index >= ninio34_erai.ninio34_index.quantile(0.90, dim='dim_0')
index_sninio_erai_lower = ninio34_erai.ninio34_index <= ninio34_erai.ninio34_index.quantile(0.10, dim='dim_0')
index_sninio_erai_normal = np.logical_and(ninio34_erai.ninio34_index < ninio34_erai.ninio34_index.quantile(0.90, dim='dim_0'), ninio34_erai.ninio34_index > ninio34_erai.ninio34_index.quantile(0.10, dim='dim_0'))

# idem s4
ninio34_s4 =  xr.open_dataset(RUTA + FILE_NINIO_S4)
# select ninio, ninia and neutral
index_ninio_s4_upper = ninio34_s4.ninio34_index >= ninio34_s4.ninio34_index.quantile(0.75, dim='dim_0')
index_ninio_s4_lower = ninio34_s4.ninio34_index <= ninio34_s4.ninio34_index.quantile(0.25, dim='dim_0')
index_ninio_s4_normal = np.logical_and(ninio34_s4.ninio34_index < ninio34_s4.ninio34_index.quantile(0.75, dim='dim_0'), ninio34_s4.ninio34_index > ninio34_s4.ninio34_index.quantile(0.25, dim='dim_0'))
index_sninio_s4_upper = ninio34_s4.ninio34_index >= ninio34_s4.ninio34_index.quantile(0.90, dim='dim_0')
index_sninio_s4_lower = ninio34_s4.ninio34_index <= ninio34_s4.ninio34_index.quantile(0.10, dim='dim_0')
index_sninio_s4_normal = np.logical_and(ninio34_s4.ninio34_index < ninio34_s4.ninio34_index.quantile(0.90, dim='dim_0'), ninio34_s4.ninio34_index > ninio34_s4.ninio34_index.quantile(0.10, dim='dim_0'))

PoV_erai = xr.open_dataset(RUTA + FILE_PV_ERAI)
PoV_s4 =  xr.open_dataset(RUTA + FILE_PV_S4)

PoV_erai_ninio = PoV_erai.SPV_mon.values[index_ninio_erai_upper.values]
PoV_erai_sninio = PoV_erai.SPV_mon.values[index_sninio_erai_upper.values]
PoV_erai_ninia = PoV_erai.SPV_mon.values[index_ninio_erai_lower.values]
PoV_erai_sninia = PoV_erai.SPV_mon.values[index_sninio_erai_lower.values]
PoV_s4_ninio = PoV_s4.SPV_index.values[index_ninio_s4_upper.values]
PoV_s4_sninio = PoV_s4.SPV_index.values[index_sninio_s4_upper.values]
PoV_s4_ninia = PoV_s4.SPV_index.values[index_ninio_s4_lower.values]
PoV_s4_sninia = PoV_s4.SPV_index.values[index_sninio_s4_lower.values]

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
# plot violin plot
axes[0].boxplot(PoV_erai_ninio.T, positions=[1.20])
axes[0].boxplot(PoV_erai_sninio.T, positions=[1.8])

axes[0].set_title('Strat. Polar Vortex index - ERAI')

# plot box plot
axes[1].violinplot(PoV_s4_ninio,positions=[1.20],
		   showmeans=False,
		   showmedians=True)
axes[1].violinplot(PoV_s4_sninio, positions=[1.8],
		   showmeans=False, showmedians=True)

axes[1].set_title('Strat. Polar Vortex index - S4')
for ax in axes:
	# adding horizontal grid lines
	ax.yaxis.grid(True)
	#ax.set_xticks(range(5))
	#ax.set_xtickslabel(season)
	ax.set_xlim([0.5, 2.5])
	#ax.set_ylim([-4.8, 4.8])
	#ax.set_xlabel('Seasons')
	ax.set_ylabel('SAM index')
# add x-tick labels
#plt.setp(axes, xticks=np.arange(1.5, 9, 1.6), xticklabels=season)

plt.savefig(FIG_PATH + 'boxplot_PoV_NINIOS.png', res=400)

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
# plot violin plot
axes[0].boxplot(PoV_erai.SPV_mon.values, positions=[1.20])
axes[0].boxplot(PoV_erai_ninio, positions=[1.8])

axes[0].set_title('Strat Polar Vortex index - ERAI')

# plot box plot
axes[1].violinplot(PoV_s4.SPV_index.values,positions=[1.20],
		   showmeans=False,
		   showmedians=True)
axes[1].violinplot(PoV_s4_ninio, positions=[1.8],
		   showmeans=False, showmedians=True)

axes[1].set_title('Strat Polar Vortex index - S4')
for ax in axes:
	# adding horizontal grid lines
	ax.yaxis.grid(True)
	#ax.set_xticks(range(5))
	#ax.set_xtickslabel(season)
	ax.set_xlim([0.5, 2.5])
	#ax.set_ylim([-4.8, 4.8])
	#ax.set_xlabel('Seasons')
	ax.set_ylabel('SAM index')
# add x-tick labels
#plt.setp(axes, xticks=np.arange(1.5, 9, 1.6), xticklabels=season)
plt.savefig(FIG_PATH + 'boxplot_PoV.png', res=400)

fig = plt.figure(figsize=(12, 12), dpi=300)
results = sm.OLS(PoV_erai.SPV_mon.values, sm.add_constant(ninio34_erai.ninio34_index.values)).fit()
x_pred = np.linspace(-5, 5, 1000)
y_pred = x_pred * results.params[1] + results.params[0]
print(results.summary())
x_pred2 = sm.add_constant(x_pred)
sdev, lower, upper = wls_prediction_std(results, exog=x_pred2, alpha=0.05)
interv = results.conf_int(alpha=0.05)

plt.scatter(ninio34_erai.ninio34_index.values, PoV_erai.SPV_mon.values)
plt.plot(x_pred, y_pred)
plt.fill_between(x_pred, lower, upper, color='#888888', alpha=0.1)
plt.xlim([-6, 6])
#plt.ylim([-4, 4])
plt.xlabel('Ninio 3.4')
plt.ylabel('PoV index')
plt.title('Regression PoV vs ENSO - ERAI - correlation: ' + '{:03.2f}'.format(np.corrcoef(ninio34_erai.ninio34_index.values, PoV_erai.SPV_mon.values) [0,1]))
plt.savefig(FIG_PATH + 'regression_pov_enso.png', res=400, box_inches='tight')
limite = np.max(np.abs(ninio34_s4.ninio34_index.values))
#scatter plot between ENSO and SAM
fig = plt.figure(figsize=(14, 12), dpi=300)
results = sm.OLS(PoV_s4.SPV_index.values, sm.add_constant(ninio34_s4.ninio34_index.values)).fit()
x_pred = np.linspace(-limite, limite, 1000)
y_pred = x_pred * results.params[1] + results.params[0]
print(results.summary())
x_pred2 = sm.add_constant(x_pred)
sdev, lower, upper = wls_prediction_std(results, exog=x_pred2, alpha=0.05)
interv = results.conf_int(alpha=0.05)

plt.scatter(ninio34_s4.ninio34_index.values, PoV_s4.SPV_index.values)
plt.plot(x_pred, y_pred)
plt.fill_between(x_pred, lower, upper, color='#888888', alpha=0.1)
plt.xlim([-limite - 0.5 , limite + 0.5])
#plt.ylim([-4, 4])
plt.xlabel('Ninio 3.4')
plt.ylabel('PoV index')
plt.suptitle('Regression PoV vs ENSO - S4 - correlation: ' + '{:03.2f}'.format(np.corrcoef(ninio34_s4.ninio34_index.values, PoV_s4.SPV_index.values)[0,1]))
plt.savefig(FIG_PATH + 'regression_pov_enso_s4.png', res=400, bbox_inches='tight' )


