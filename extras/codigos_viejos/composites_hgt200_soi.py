#sort el ninio events
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import mpl_toolkits.basemap as bm
RUTA='~/data/'
#abro el archivo de soi y clasifico las realizaciones
soi_monthly = xr.open_dataset(RUTA + 'soi_monthly.nc')
#
##abro el archivo con los cuartiles
soi_monthly_upper = xr.open_dataset(RUTA + 'soi_monthly_upper.nc')
soi_monthly_lower = xr.open_dataset(RUTA + 'soi_monthly_lower.nc')
#
#abro el archivo de geopotencial y junto la coordenada year y numbre
ds = xr.open_dataset(RUTA + 'hgt200.nc')
ds = ds.stack(realiz = ['year', 'number'])
ds = ds.transpose('month', 'realiz', 'latitude', 'longitude')
#
##composites monthly values
index_monthly_upper = soi_monthly > soi_monthly_upper
#print(np.sum(index_monthly_upper,axis=0))
index_monthly_lower = soi_monthly < soi_monthly_lower
#print(np.sum(index_monthly_lower,axis=0))

index_monthly_normal = np.logical_and(soi_monthly > soi_monthly_lower,
				      soi_monthly < soi_monthly_upper)
#print(np.sum(index_monthly_normal,axis=0))

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
    var = np.mean(ds.z[i, index_monthly_upper.msl[0, :], :, :], axis=0) -\
	  np.mean(ds.z[i, index_monthly_normal.msl[0, :], :, :], axis=0)
    CS1 = mapproj.contourf(lonproj, latproj, var, clevs, cmap=barra, extend='both') 
    barra.set_under(barra(0))
    barra.set_over(barra(barra.N-1))
    plt.title('+SOI YEARS')

    plt.subplot(2, 1, 2)
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
    var = np.mean(ds.z[i, index_monthly_lower.msl[0, :], :, :], axis=0) -\
	  np.mean(ds.z[i, index_monthly_normal.msl[0, :], :, :], axis=0)
    CS1 = mapproj.contourf(lonproj, latproj, var, clevs, cmap=barra, extend='both') 
    barra.set_under(barra(0))
    barra.set_over(barra(barra.N-1))
    plt.title('-SOI YEARS')
    plt.suptitle('Composites S4 HGT 200hPa - ' + month[i], fontsize=12, x=0.52,
	         y=0.9)
    #fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(bottom=0.2, top=0.82, hspace=0.2, wspace=0.05)
    cbar_ax = fig.add_axes([0.39, 0.1, 0.25, 0.05])
    fig.colorbar(CS1, cax=cbar_ax, orientation='horizontal')
    #cbar_ax.set_xticklabels([0, 40, 80, 120, 160], size=9)
    plt.savefig('hgt_200_composites_SOI_' + month[i] +'.png', dpi=300,
		bbox_inches='tight', orientation='landscape',
		papertype='A4')
for i in np.arange(0,4):
    fig = plt.figure(4 + i + 1, figsize=(10,6.7), dpi=300)  #fig size in inches
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
    var = np.mean(ds.z[i, index_monthly_lower.msl[0, :], :, :], axis=0) -\
	  np.mean(ds.z[i, index_monthly_upper.msl[0, :], :, :], axis=0)
    CS1 = mapproj.contourf(lonproj, latproj, var, clevs, cmap=barra, extend='both') 
    barra.set_under(barra(0))
    barra.set_over(barra(barra.N-1))
    plt.title('Composites diff -SOI - +SOI YEARS - ' + month[i])
    fig.colorbar(CS1, orientation='horizontal')
    plt.savefig('hgt_200_composites_+-SOI_' + month[i] +'.png', dpi=300,
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
        var = ds.sel(**{'month':slice(8,10)}).mean(dim='month')
    else:
        var = ds.sel(**{'month':slice(9,11)}).mean(dim='month')
    var = np.mean(var.z[index_monthly_upper.msl[0, :], :, :], axis=0) -\
	 	 np.mean(var.z[index_monthly_normal.msl[0, :], :, :], axis=0)	
    CS1 = mapproj.contourf(lonproj, latproj, var, clevs, cmap=barra, extend='both') 
    barra.set_under(barra(0))
    barra.set_over(barra(barra.N-1))
    plt.title('+SOI YEARS')
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
        var = ds.sel(**{'month':slice(8,10)}).mean(dim='month')
    else:
        var = ds.sel(**{'month':slice(9,11)}).mean(dim='month')
    var = np.mean(var.z[index_monthly_lower.msl[0, :], :, :], axis=0) -\
	 	 np.mean(var.z[index_monthly_normal.msl[0, :], :, :], axis=0)	
    CS1 = mapproj.contourf(lonproj, latproj, var, clevs, cmap=barra, extend='both') 
    barra.set_under(barra(0))
    barra.set_over(barra(barra.N-1))
    plt.title('-SOI YEARS')
    plt.suptitle('Composites S4 HGT 200hPa - ' + seas[i], fontsize=12, x=0.52,
	         y=0.9)
    fig.subplots_adjust(bottom=0.2, top=0.82, hspace=0.2, wspace=0.05)
    cbar_ax = fig.add_axes([0.39, 0.1, 0.25, 0.05])
    fig.colorbar(CS1, cax=cbar_ax, orientation='horizontal')
    plt.savefig('hgt_200_composites_SOI_' + seas[i] +'.png', dpi=300,
		bbox_inches='tight', orientation='landscape',
		papertype='A4')
for i in np.arange(0,2):
    fig = plt.figure(10 + i + 1, figsize=(10,6.7), dpi=300)  #fig size in inches
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
        var = ds.sel(**{'month':slice(8,10)}).mean(dim='month')
    else:
        var = ds.sel(**{'month':slice(9,11)}).mean(dim='month')
    var = np.mean(var.z[index_monthly_lower.msl[0, :], :, :], axis=0) -\
	 	 np.mean(var.z[index_monthly_upper.msl[0, :], :, :], axis=0)	
    CS1 = mapproj.contourf(lonproj, latproj, var, clevs, cmap=barra, extend='both') 
    barra.set_under(barra(0))
    barra.set_over(barra(barra.N-1))
    plt.title('Composites diff -SOI - +SOI YEARS ' + seas[i])
    fig.colorbar(CS1, orientation='horizontal')
    plt.savefig('hgt_200_composites_+-SOI_' + seas[i] +'.png', dpi=300,
		bbox_inches='tight', orientation='landscape',
		papertype='A4')

