#sort el ninio events
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import mpl_toolkits.basemap as bm
RUTA='~/data/'
#abro el archivo de soi y clasifico las realizaciones
soi_monthly = xr.open_dataset(RUTA + 'soi_monthly.nc')
PV_monthly = xr.open_dataset(RUTA + 'PV_monthly2.nc')
#
##abro el archivo con los cuartiles
soi_monthly_upper = xr.open_dataset(RUTA + 'soi_monthly_upper.nc')
soi_monthly_lower = xr.open_dataset(RUTA + 'soi_monthly_lower.nc')
#abro el archivo de geopotencial y junto la coordenada year y numbre
ds = xr.open_dataset(RUTA + 'hgt200.nc')
ds = ds.stack(realiz = ['year', 'number'])
ds = ds.transpose('month', 'realiz', 'latitude', 'longitude')
ds_aso = ds.sel(**{'month':slice(8,10)}).mean(dim='month')
ds_son = ds.sel(**{'month':slice(9,11)}).mean(dim='month')

#
##composites monthly values
index_soi_upper = soi_monthly.msl[0,:] > soi_monthly_upper.msl[0]
index_soi_lower = soi_monthly.msl[0,:] < soi_monthly_lower.msl[0]

PV_upper_soi_upper = PV_monthly.z[:, index_soi_upper.values].quantile(0.75, dim='realiz', interpolation='linear')
PV_lower_soi_upper = PV_monthly.z[:, index_soi_upper.values].quantile(0.25, dim='realiz', interpolation='linear')
index_PV_lower_soi_upper = PV_monthly.z < PV_lower_soi_upper
index_PV_upper_soi_upper = PV_monthly.z < PV_upper_soi_upper

PV_upper_soi_lower = PV_monthly.z[:, index_soi_lower.values].quantile(0.75, dim='realiz', interpolation='linear')
PV_lower_soi_lower = PV_monthly.z[:, index_soi_lower.values].quantile(0.25, dim='realiz', interpolation='linear')
index_PV_lower_soi_lower = PV_monthly.z < PV_lower_soi_lower
index_PV_upper_soi_lower = PV_monthly.z < PV_upper_soi_lower

PV_aso = PV_monthly.sel(**{'month':slice(8,10)}).mean(dim='month')
PV_son = PV_monthly.sel(**{'month':slice(9,11)}).mean(dim='month')

PV_aso_upper_soi_upper = PV_aso.z[index_soi_upper.values].quantile(0.75, dim='realiz', interpolation='linear')
PV_aso_lower_soi_upper = PV_aso.z[index_soi_upper.values].quantile(0.25, dim='realiz', interpolation='linear')
index_aso_PV_lower_soi_upper = PV_aso.z < PV_aso_lower_soi_upper
index_aso_PV_upper_soi_upper = PV_aso.z < PV_aso_upper_soi_upper

PV_aso_upper_soi_lower = PV_aso.z[index_soi_lower.values].quantile(0.75, dim='realiz', interpolation='linear')
PV_aso_lower_soi_lower = PV_aso.z[index_soi_lower.values].quantile(0.25, dim='realiz', interpolation='linear')
index_aso_PV_lower_soi_lower = PV_aso.z < PV_aso_lower_soi_lower
index_aso_PV_upper_soi_lower = PV_aso.z < PV_aso_upper_soi_lower
#====
PV_son_upper_soi_upper = PV_son.z[index_soi_upper.values].quantile(0.75, dim='realiz', interpolation='linear')
PV_son_lower_soi_upper = PV_son.z[index_soi_upper.values].quantile(0.25, dim='realiz', interpolation='linear')
index_son_PV_lower_soi_upper = PV_son.z < PV_son_lower_soi_upper
index_son_PV_upper_soi_upper = PV_son.z < PV_son_upper_soi_upper

PV_son_upper_soi_lower = PV_son.z[index_soi_lower.values].quantile(0.75, dim='realiz', interpolation='linear')
PV_son_lower_soi_lower = PV_son.z[index_soi_lower.values].quantile(0.25, dim='realiz', interpolation='linear')
index_son_PV_lower_soi_lower = PV_son.z < PV_son_lower_soi_lower
index_son_PV_upper_soi_lower = PV_son.z < PV_son_upper_soi_lower

month = ['Aug', 'Sep', 'Oct', 'Nov']
for i in np.arange(0,4):
    fig = plt.figure(i + 1, figsize=(10,6.7), dpi=300)  #fig size in inches
    plt.subplot(2, 1, 1)
    mapproj = bm.Basemap(projection='cyl',
                         llcrnrlat=-88.0, llcrnrlon=0.0,
                         urcrnrlat=10, urcrnrlon=359.5)    #projection and map limits
    
    mapproj.drawcoastlines()          # coast
    mapproj.drawparallels(np.array([-75, -60, -45, -30,-15, 0]), labels=[1,0,0,0])    #draw parallels
    mapproj.drawmeridians(np.array([-180, -90, 0, 90, 180]), labels=[0,0,0,1])     #draw meridians
    mapproj.drawcountries()
    lonall, latall = np.meshgrid(ds.longitude.values, ds.latitude.values)          #array of grid
    lonproj, latproj = mapproj(lonall, latall)      #poject grid
    
    clevs = np.arange(-70, 80, 10)
    barra = plt.cm.RdBu_r #colorbar
    var = np.mean(ds.z[i, np.logical_and(index_soi_upper, index_PV_lower_soi_upper[i, :]), :, :], axis=0) -\
	  np.mean(ds.z[i, np.logical_and(index_soi_upper, index_PV_upper_soi_upper[i, :]), :, :], axis=0)
    CS1 = mapproj.contourf(lonproj, latproj, var, clevs, cmap=barra, extend='both') 
    barra.set_under(barra(0))
    barra.set_over(barra(barra.N-1))
    plt.title('Composites S-PV - W-PV cond +SOI YEARS')
#
    plt.subplot(2, 1, 2)
    mapproj = bm.Basemap(projection='cyl',
                         llcrnrlat=-88.0, llcrnrlon=0.0,
                         urcrnrlat=10, urcrnrlon=359.5)    #projection and map limits
#    
    mapproj.drawcoastlines()          # coast
    mapproj.drawparallels(np.array([-75, -60, -45, -30,-15, 0]), labels=[1,0,0,0])    #draw parallels
    mapproj.drawmeridians(np.array([-180, -90, 0, 90, 180]), labels=[0,0,0,1])     #draw meridians
    mapproj.drawcountries()
    lonall, latall = np.meshgrid(ds.longitude.values, ds.latitude.values)          #array of grid
    lonproj, latproj = mapproj(lonall, latall)      #poject grid
    clevs = np.arange(-70, 80, 10)
    barra = plt.cm.RdBu_r #colorbar
    var = np.mean(ds.z[i, np.logical_and(index_soi_lower, index_PV_lower_soi_lower[i, :]), :, :], axis=0) -\
	  np.mean(ds.z[i,np.logical_and(index_soi_lower, index_PV_upper_soi_lower[i, :]), :, :], axis=0)
    CS1 = mapproj.contourf(lonproj, latproj, var, clevs, cmap=barra, extend='both') 
    barra.set_under(barra(0))
    barra.set_over(barra(barra.N-1))
    plt.title('Composites S-PV - W-PV cond -SOI YEARS')
    plt.suptitle('Composites S4 HGT 200hPa - ' + month[i], fontsize=12, x=0.52,
	         y=0.9)
    #fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(bottom=0.2, top=0.82, hspace=0.25, wspace=0.05)
    cbar_ax = fig.add_axes([0.39, 0.1, 0.25, 0.05])
    fig.colorbar(CS1, cax=cbar_ax, orientation='horizontal')
    #cbar_ax.set_xticklabels([0, 40, 80, 120, 160], size=9)
    plt.savefig('hgt_200_composites_cond2_SOI_PV_' + month[i] +'.png', dpi=300,
		bbox_inches='tight', orientation='landscape',
		papertype='A4')
seas = ['ASO', 'SON']
for i in np.arange(0,2):
    fig = plt.figure(8 + i + 1, figsize=(10,6.7), dpi=300)  #fig size in inches
    plt.subplot(2, 1, 1)
    mapproj = bm.Basemap(projection='cyl',
                         llcrnrlat=-88.0, llcrnrlon=0.0,
                         urcrnrlat=10, urcrnrlon=359.5)    #projection and map limits
    
    mapproj.drawcoastlines()          # coast
    mapproj.drawparallels(np.array([-75, -60, -45, -30,-15, 0]), labels=[1,0,0,0])    #draw parallels
    mapproj.drawmeridians(np.array([-180, -90, 0, 90, 180]), labels=[0,0,0,1])     #draw meridians
    mapproj.drawcountries()
    lonall, latall = np.meshgrid(ds.longitude.values, ds.latitude.values)          #array of grid
    lonproj, latproj = mapproj(lonall, latall)      #poject grid
    
    clevs = np.arange(-70, 80, 10)
    barra = plt.cm.RdBu_r #colorbar
    if i == 0:
        var = np.mean(ds_aso.z[np.logical_and(index_soi_upper, index_aso_PV_lower_soi_upper), :, :], axis=0) -\
	      np.mean(ds_aso.z[np.logical_and(index_soi_upper, index_aso_PV_upper_soi_upper), :, :], axis=0)
    else:
         var = np.mean(ds_son.z[np.logical_and(index_soi_upper, index_son_PV_lower_soi_upper), :, :], axis=0) -\
	       np.mean(ds_son.z[np.logical_and(index_soi_upper, index_son_PV_upper_soi_upper), :, :], axis=0)
    CS1 = mapproj.contourf(lonproj, latproj, var, clevs, cmap=barra, extend='both') 
    barra.set_under(barra(0))
    barra.set_over(barra(barra.N-1))
    plt.title('Composites S-PV - W-PV cond +SOI YEARS')
#
    plt.subplot(2, 1, 2)
    mapproj = bm.Basemap(projection='cyl',
                         llcrnrlat=-88.0, llcrnrlon=0.0,
                         urcrnrlat=10, urcrnrlon=359.5)    #projection and map limits
#    
    mapproj.drawcoastlines()          # coast
    mapproj.drawparallels(np.array([-75, -60, -45, -30,-15, 0]), labels=[1,0,0,0])    #draw parallels
    mapproj.drawmeridians(np.array([-180, -90, 0, 90, 180]), labels=[0,0,0,1])     #draw meridians
    mapproj.drawcountries()
    lonall, latall = np.meshgrid(ds.longitude.values, ds.latitude.values)          #array of grid
    lonproj, latproj = mapproj(lonall, latall)      #poject grid
#    
    clevs = np.arange(-120, 140, 20)
    barra = plt.cm.RdBu_r #colorbar
    if i == 0:
        var = np.mean(ds_aso.z[np.logical_and(index_soi_lower, index_aso_PV_lower_soi_lower), :, :], axis=0) -\
	      np.mean(ds_aso.z[np.logical_and(index_soi_lower, index_aso_PV_upper_soi_lower), :, :], axis=0)
    else:
         var = np.mean(ds_son.z[np.logical_and(index_soi_upper, index_son_PV_lower_soi_lower), :, :], axis=0) -\
	       np.mean(ds_son.z[np.logical_and(index_soi_upper, index_son_PV_upper_soi_lower), :, :], axis=0)

    CS1 = mapproj.contourf(lonproj, latproj, var, clevs, cmap=barra, extend='both') 
    barra.set_under(barra(0))
    barra.set_over(barra(barra.N-1))
    plt.title('Composites S-PV - W-PV cond -SOI YEARS')
    plt.suptitle('Composites S4 HGT 200hPa - ' + seas[i], fontsize=12, x=0.52,
	         y=0.9)
    #fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(bottom=0.2, top=0.82, hspace=0.25, wspace=0.05)
    cbar_ax = fig.add_axes([0.39, 0.1, 0.25, 0.05])
    fig.colorbar(CS1, cax=cbar_ax, orientation='horizontal')
    #cbar_ax.set_xticklabels([0, 40, 80, 120, 160], size=9)
    plt.savefig('hgt_200_composites_cond2_SOI_PV_' + seas[i] +'.png', dpi=300,
		bbox_inches='tight', orientation='landscape',
		papertype='A4')

