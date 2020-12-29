#compute boxplot of spv for ninios
#compute regression between spv and enso
import numpy as np
import datetime
import pandas as pd
import xarray as xr
import seaborn as sns
from eofs.standard import Eof
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import matplotlib.pyplot as plt
import cartopy.crs as ccrs	
import cartopy.feature 
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

RUTA = '/home/users/vg140344/datos/data/fogt/' #este spanglish no te lo robo amiga
FIG_PATH = '/home/users/vg140344/assessment_SH_zonal_asymmetries/figures_paper/'
month = ['Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']
seas = ['ASO', 'SON', 'OND', 'NDJ', 'DJF']

fig = plt.figure(figsize=(14, 12), dpi=500)
for i in np.arange(0, 7):
	f = np.load(RUTA + 'ASL_coords' + month[i] + '.npz')
	ax = plt.subplot(4, 2, i + 1)
	plt.scatter(f['var9'][2], f['var9'][1], marker="o", c='tab:green',
		    alpha=0.3, label='All PoV' )
	plt.scatter(f['var6'][2], f['var6'][1], marker="o", c='tab:red',
		    alpha=0.3,  label='Strong PoV')
	plt.scatter(f['var3'][2], f['var3'][1], marker="o", c='tab:blue',
		    alpha=0.3,  label='Weak PoV')
	#plt.xlim([180, 300])
	#plt.ylim([-75, -55])
	plt.xlabel('Longitude')
	plt.ylabel('Latitude')
	plt.title(month[i])
plt.legend()
plt.suptitle('Location Min HGT - Ninia')
plt.savefig(FIG_PATH + 'Min_monthly_location_ninia.png', res=500,bbox_inches='tight' )
plt.close('all')

fig = plt.figure(figsize=(14, 12), dpi=500)
for i in np.arange(0, 7):
	f = np.load(RUTA + 'ASL_coords' + month[i] + '.npz')
	ax = plt.subplot(4, 2, i + 1)
	plt.scatter(f['var9'][5], f['var9'][4], marker="o", c='tab:green',
		    alpha=0.3, label='All PoV' )
	plt.scatter(f['var6'][5], f['var6'][4], marker="o", c='tab:red',
		    alpha=0.3,  label='Strong PoV')
	plt.scatter(f['var3'][5], f['var3'][4], marker="o", c='tab:blue',
		    alpha=0.3,  label='Weak PoV')
	#plt.xlim([180, 300])
	#plt.ylim([-75, -55])
	plt.xlabel('Longitude')
	plt.ylabel('Latitude')
	plt.title(month[i])
plt.legend()
plt.suptitle('Location Max HGT - Ninia')
plt.savefig(FIG_PATH + 'Max_monthly_location_ninia.png', res=500,bbox_inches='tight' )
plt.close('all')

fig = plt.figure(figsize=(14, 12), dpi=500)
for i in np.arange(0, 5):
	f = np.load(RUTA + 'ASL_coords' + seas[i] + '.npz')
	ax = plt.subplot(3, 2, i + 1)
	plt.scatter(f['var9'][5, :], f['var9'][4, :], marker="o", c='tab:green',
		    alpha=0.3, label='All PoV' )
	plt.scatter(f['var6'][5, :], f['var6'][4, :], marker="o", c='tab:red',
		    alpha=0.3, label='Strong PoV')
	plt.scatter(f['var3'][5, :], f['var3'][4, :], marker="o", c='tab:blue',
		    alpha=0.3, label='Weak PoV')
	#plt.xlim([180, 300])
	#plt.ylim([-75, -55])
	plt.xlabel('Longitude')
	plt.ylabel('Latitude')
	plt.title(seas[i])
plt.legend()
plt.suptitle('Location Max HGT - Ninia')
plt.savefig(FIG_PATH + 'Max_seasonal_location_ninia.png', res=500,bbox_inches='tight' )
plt.close('all')
#
fig = plt.figure(figsize=(14, 12), dpi=500)
for i in np.arange(0, 5):
	f = np.load(RUTA + 'ASL_coords' + seas[i] + '.npz')
	ax = plt.subplot(3, 2, i + 1)
	plt.scatter(f['var9'][2, :], f['var9'][1, :], marker="o", c='tab:green',
		    alpha=0.3, label='All PoV' )
	plt.scatter(f['var6'][2, :], f['var6'][1, :], marker="o", c='tab:red',
		    alpha=0.3, label='Strong PoV')
	plt.scatter(f['var3'][2, :], f['var3'][1, :], marker="o", c='tab:blue',
		    alpha=0.3, label='Weak PoV')
	#plt.xlim([180, 300])
	#plt.ylim([-75, -55])
	plt.xlabel('Longitude')
	plt.ylabel('Latitude')
	plt.title(seas[i])
plt.legend()
plt.suptitle('Location Min HGT - Ninia')
plt.savefig(FIG_PATH + 'Min_seasonal_location_ninia.png', res=500,bbox_inches='tight' )
plt.close('all')

#=================================================================================
df_max_all = pd.DataFrame()
df_min_all = pd.DataFrame()
df_max_SPV = pd.DataFrame()
df_min_SPV = pd.DataFrame()
df_max_WPV = pd.DataFrame()
df_min_WPV = pd.DataFrame()

#Intensity
for i in np.arange(0, 7):
	f = np.load(RUTA + 'ASL_coords' + month[i] + '.npz')
	df = pd.DataFrame({month[i]: f['var9'][3, :]})
	df_max_all = pd.concat([df_max_all, df], axis=0)
	df = pd.DataFrame({month[i]: f['var6'][3, :]})
	df_max_SPV = pd.concat([df_max_SPV, df], axis=0)
	df = pd.DataFrame({month[i]: f['var3'][3, :]})
	df_max_WPV = pd.concat([df_max_WPV, df], axis=0)
	df = pd.DataFrame({month[i]: f['var9'][0, :]})
	df_min_all = pd.concat([df_min_all, df], axis=0)
	df = pd.DataFrame({month[i]: f['var6'][0, :]})
	df_min_SPV = pd.concat([df_min_SPV, df], axis=0)
	df = pd.DataFrame({month[i]: f['var3'][0, :]})
	df_min_WPV = pd.concat([df_min_WPV, df], axis=0)
df_max_all = pd.melt(df_max_all, value_vars=month)
df_max_SPV = pd.melt(df_max_SPV, value_vars=month)
df_max_WPV = pd.melt(df_max_WPV, value_vars=month)
df_min_all = pd.melt(df_min_all, value_vars=month)
df_min_SPV = pd.melt(df_min_SPV, value_vars=month)
df_min_WPV = pd.melt(df_min_WPV, value_vars=month)

df_max_all['Vortex'] = 'All'
df_max_SPV['Vortex'] = 'SPoV'
df_max_WPV['Vortex'] = 'WPoV'
df_min_all['Vortex'] = 'All'
df_min_SPV['Vortex'] = 'SPoV'
df_min_WPV['Vortex'] = 'WPoV'

dfMax = pd.concat([df_max_all, df_max_SPV, df_max_WPV], axis=0)
dfMin = pd.concat([df_min_all, df_min_SPV, df_min_WPV], axis=0)
dfMax = dfMax.reindex()
print(dfMax)
fig = plt.figure(figsize=(14, 12), dpi=500)
ax = plt.subplot(1, 2, 1)
sns.boxplot(x='variable', y='value', data=dfMax, hue='Vortex')
plt.title('Maximum HGT Value')
ax = plt.subplot(1, 2, 2)
sns.boxplot(x='variable', y='value', data=dfMin, hue='Vortex')
plt.title('Minimum HGT Value')
plt.savefig(FIG_PATH + 'monthly_intesity_HGT_ninia.png', res=500)
#=========================================================================
df_max_all = pd.DataFrame()
df_min_all = pd.DataFrame()
df_max_SPV = pd.DataFrame()
df_min_SPV = pd.DataFrame()
df_max_WPV = pd.DataFrame()
df_min_WPV = pd.DataFrame()

#Intensity
for i in np.arange(0, 5):
	f = np.load(RUTA + 'ASL_coords' + seas[i] + '.npz')
	df = pd.DataFrame({seas[i]: f['var9'][3, :]})
	df_max_all = pd.concat([df_max_all, df], axis=0)
	df = pd.DataFrame({seas[i]: f['var6'][3, :]})
	df_max_SPV = pd.concat([df_max_SPV, df], axis=0)
	df = pd.DataFrame({seas[i]: f['var3'][3, :]})
	df_max_WPV = pd.concat([df_max_WPV, df], axis=0)
	df = pd.DataFrame({seas[i]: f['var9'][0, :]})
	df_min_all = pd.concat([df_min_all, df], axis=0)
	df = pd.DataFrame({seas[i]: f['var6'][0, :]})
	df_min_SPV = pd.concat([df_min_SPV, df], axis=0)
	df = pd.DataFrame({seas[i]: f['var3'][0, :]})
	df_min_WPV = pd.concat([df_min_WPV, df], axis=0)
df_max_all = pd.melt(df_max_all, value_vars=seas)
df_max_SPV = pd.melt(df_max_SPV, value_vars=seas)
df_max_WPV = pd.melt(df_max_WPV, value_vars=seas)
df_min_all = pd.melt(df_min_all, value_vars=seas)
df_min_SPV = pd.melt(df_min_SPV, value_vars=seas)
df_min_WPV = pd.melt(df_min_WPV, value_vars=seas)

df_max_all['Vortex'] = 'All'
df_max_SPV['Vortex'] = 'SPoV'
df_max_WPV['Vortex'] = 'WPoV'
df_min_all['Vortex'] = 'All'
df_min_SPV['Vortex'] = 'SPoV'
df_min_WPV['Vortex'] = 'WPoV'

dfMax = pd.concat([df_max_all, df_max_SPV, df_max_WPV], axis=0)
dfMin = pd.concat([df_min_all, df_min_SPV, df_min_WPV], axis=0)
dfMax = dfMax.reindex()
print(dfMax)
fig = plt.figure(figsize=(14, 12), dpi=500)
ax = plt.subplot(1, 2, 1)
sns.boxplot(x='variable', y='value', data=dfMax, hue='Vortex')
plt.title('Maximum HGT Value')
ax = plt.subplot(1, 2, 2)
sns.boxplot(x='variable', y='value', data=dfMin, hue='Vortex')
plt.title('Minimum HGT Value')
plt.savefig(FIG_PATH + 'seasonal_intesity_HGT_ninia.png', res=500)

##========================================================================

fig = plt.figure(figsize=(14, 12), dpi=500)
for i in np.arange(0, 7):
	f = np.load(RUTA + 'ASL_coords' + month[i] + '.npz')
	ax = plt.subplot(4, 2, i + 1)
	plt.scatter(f['var7'][2, :], f['var7'][1, :], marker="o", c='tab:green',
		    alpha=0.3, label='All PoV' )
	plt.scatter(f['var4'][2, :], f['var4'][1, :], marker="o", c='tab:red',
		    alpha=0.3, label='Strong PoV')
	plt.scatter(f['var1'][2, :], f['var1'][1, :], marker="o", c='tab:blue',
		    alpha=0.3, label='Weak PoV')
	#plt.xlim([180, 300])
	#plt.ylim([-75, -55])
	plt.xlabel('Longitude')
	plt.ylabel('Latitude')
	plt.title(month[i])
plt.legend()
plt.suptitle('Location Mim HGT - Ninio')
plt.savefig(FIG_PATH + 'Min_monthly_location_ninio.png', res=500,bbox_inches='tight' )
plt.close('all')

fig = plt.figure(figsize=(14, 12), dpi=500)
for i in np.arange(0, 5):
	f = np.load(RUTA + 'ASL_coords' + seas[i] + '.npz')
	ax = plt.subplot(4, 2, i + 1)
	plt.scatter(f['var7'][5, :], f['var7'][4, :], marker="o", c='tab:green',
		    alpha=0.3, label='All PoV' )
	plt.scatter(f['var4'][5, :], f['var4'][4, :], marker="o", c='tab:red',
		    alpha=0.3, label='Strong PoV')
	plt.scatter(f['var1'][5, :], f['var1'][4, :], marker="o", c='tab:blue',
		    alpha=0.3, label='Weak PoV')
	#plt.xlim([180, 300])
	#plt.ylim([-75, -55])
	plt.xlabel('Longitude')
	plt.ylabel('Latitude')
	plt.title(seas[i])
plt.legend()
plt.suptitle('Location Max HGT - Ninio')
plt.savefig(FIG_PATH + 'Max_seasonal_location_ninio.png', res=500,bbox_inches='tight' )
plt.close('all')

fig = plt.figure(figsize=(14, 12), dpi=500)
for i in np.arange(0, 7):
	f = np.load(RUTA + 'ASL_coords' + month[i] + '.npz')
	ax = plt.subplot(4, 2, i + 1)
	plt.scatter(f['var7'][5, :], f['var7'][4, :], marker="o", c='tab:green',
		    alpha=0.3, label='All PoV' )
	plt.scatter(f['var4'][5, :], f['var4'][4, :], marker="o", c='tab:red',
		    alpha=0.3, label='Strong PoV')
	plt.scatter(f['var1'][5, :], f['var1'][4, :], marker="o", c='tab:blue',
		    alpha=0.3, label='Weak PoV')
	#plt.xlim([180, 300])
	#plt.ylim([-75, -55])
	plt.xlabel('Longitude')
	plt.ylabel('Latitude')
	plt.title(month[i])
plt.legend()
plt.suptitle('Location Max HGT - Ninio')
plt.savefig(FIG_PATH + 'Max_monthly_location_ninio.png', res=500,bbox_inches='tight' )
plt.close('all')

fig = plt.figure(figsize=(14, 12), dpi=500)
for i in np.arange(0, 5):
	f = np.load(RUTA + 'ASL_coords' + seas[i] + '.npz')
	ax = plt.subplot(4, 2, i + 1)
	plt.scatter(f['var7'][2, :], f['var7'][1, :], marker="o", c='tab:green',
		    alpha=0.3, label='All PoV' )
	plt.scatter(f['var4'][2, :], f['var4'][1, :], marker="o", c='tab:red',
		    alpha=0.3, label='Strong PoV')
	plt.scatter(f['var1'][2, :], f['var1'][1, :], marker="o", c='tab:blue',
		    alpha=0.3, label='Weak PoV')
	#plt.xlim([180, 300])
	#plt.ylim([-75, -55])
	plt.xlabel('Longitude')
	plt.ylabel('Latitude')
	plt.title(seas[i])
plt.legend()
plt.suptitle('Location Min HGT - Ninio')
plt.savefig(FIG_PATH + 'Min_seasonal_location_ninio.png', res=500,bbox_inches='tight' )
plt.close('all')
#=================================================================================
df_max_all = pd.DataFrame()
df_min_all = pd.DataFrame()
df_max_SPV = pd.DataFrame()
df_min_SPV = pd.DataFrame()
df_max_WPV = pd.DataFrame()
df_min_WPV = pd.DataFrame()

#Intensity
for i in np.arange(0, 7):
	f = np.load(RUTA + 'ASL_coords' + month[i] + '.npz')
	df = pd.DataFrame({month[i]: f['var7'][3, :]})
	df_max_all = pd.concat([df_max_all, df], axis=0)
	df = pd.DataFrame({month[i]: f['var4'][3, :]})
	df_max_SPV = pd.concat([df_max_SPV, df], axis=0)
	df = pd.DataFrame({month[i]: f['var1'][3, :]})
	df_max_WPV = pd.concat([df_max_WPV, df], axis=0)
	df = pd.DataFrame({month[i]: f['var7'][0, :]})
	df_min_all = pd.concat([df_min_all, df], axis=0)
	df = pd.DataFrame({month[i]: f['var4'][0, :]})
	df_min_SPV = pd.concat([df_min_SPV, df], axis=0)
	df = pd.DataFrame({month[i]: f['var1'][0, :]})
	df_min_WPV = pd.concat([df_min_WPV, df], axis=0)
df_max_all = pd.melt(df_max_all, value_vars=month)
df_max_SPV = pd.melt(df_max_SPV, value_vars=month)
df_max_WPV = pd.melt(df_max_WPV, value_vars=month)
df_min_all = pd.melt(df_min_all, value_vars=month)
df_min_SPV = pd.melt(df_min_SPV, value_vars=month)
df_min_WPV = pd.melt(df_min_WPV, value_vars=month)

df_max_all['Vortex'] = 'All'
df_max_SPV['Vortex'] = 'SPoV'
df_max_WPV['Vortex'] = 'WPoV'
df_min_all['Vortex'] = 'All'
df_min_SPV['Vortex'] = 'SPoV'
df_min_WPV['Vortex'] = 'WPoV'

dfMax = pd.concat([df_max_all, df_max_SPV, df_max_WPV], axis=0)
dfMin = pd.concat([df_min_all, df_min_SPV, df_min_WPV], axis=0)
dfMax = dfMax.reindex()
print(dfMax)
fig = plt.figure(figsize=(14, 12), dpi=500)
ax = plt.subplot(1, 2, 1)
sns.boxplot(x='variable', y='value', data=dfMax, hue='Vortex')
plt.title('Maximum HGT Value')
ax = plt.subplot(1, 2, 2)
sns.boxplot(x='variable', y='value', data=dfMin, hue='Vortex')
plt.title('Minimum HGT Value')
plt.savefig(FIG_PATH + 'monthly_intesity_HGT_ninio.png', res=500)
#=========================================================================
df_max_all = pd.DataFrame()
df_min_all = pd.DataFrame()
df_max_SPV = pd.DataFrame()
df_min_SPV = pd.DataFrame()
df_max_WPV = pd.DataFrame()
df_min_WPV = pd.DataFrame()

#Intensity
for i in np.arange(0, 5):
	f = np.load(RUTA + 'ASL_coords' + seas[i] + '.npz')
	df = pd.DataFrame({seas[i]: f['var7'][3, :]})
	df_max_all = pd.concat([df_max_all, df], axis=0)
	df = pd.DataFrame({seas[i]: f['var4'][3, :]})
	df_max_SPV = pd.concat([df_max_SPV, df], axis=0)
	df = pd.DataFrame({seas[i]: f['var1'][3, :]})
	df_max_WPV = pd.concat([df_max_WPV, df], axis=0)
	df = pd.DataFrame({seas[i]: f['var7'][0, :]})
	df_min_all = pd.concat([df_min_all, df], axis=0)
	df = pd.DataFrame({seas[i]: f['var4'][0, :]})
	df_min_SPV = pd.concat([df_min_SPV, df], axis=0)
	df = pd.DataFrame({seas[i]: f['var1'][0, :]})
	df_min_WPV = pd.concat([df_min_WPV, df], axis=0)
df_max_all = pd.melt(df_max_all, value_vars=seas)
df_max_SPV = pd.melt(df_max_SPV, value_vars=seas)
df_max_WPV = pd.melt(df_max_WPV, value_vars=seas)
df_min_all = pd.melt(df_min_all, value_vars=seas)
df_min_SPV = pd.melt(df_min_SPV, value_vars=seas)
df_min_WPV = pd.melt(df_min_WPV, value_vars=seas)

df_max_all['Vortex'] = 'All'
df_max_SPV['Vortex'] = 'SPoV'
df_max_WPV['Vortex'] = 'WPoV'
df_min_all['Vortex'] = 'All'
df_min_SPV['Vortex'] = 'SPoV'
df_min_WPV['Vortex'] = 'WPoV'

dfMax = pd.concat([df_max_all, df_max_SPV, df_max_WPV], axis=0)
dfMin = pd.concat([df_min_all, df_min_SPV, df_min_WPV], axis=0)
dfMax = dfMax.reindex()
print(dfMax)
fig = plt.figure(figsize=(14, 12), dpi=500)
ax = plt.subplot(1, 2, 1)
sns.boxplot(x='variable', y='value', data=dfMax, hue='Vortex')
plt.title('Maximum HGT Value')
ax = plt.subplot(1, 2, 2)
sns.boxplot(x='variable', y='value', data=dfMin, hue='Vortex')
plt.title('Minimum HGT Value')
plt.savefig(FIG_PATH + 'seasonal_intesity_HGT_ninio.png', res=500)


