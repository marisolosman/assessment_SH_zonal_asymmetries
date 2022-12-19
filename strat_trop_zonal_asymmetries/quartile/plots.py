import numpy as np
from numpy import ma
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


def PlotCompDiff(var, lat, lon, title, filename):
        proj = ccrs.PlateCarree(central_longitude=180)
        fig = plt.figure(1, (10, 3.7), 300)
        ax = plt.subplot(projection=proj)
        clevs = np.arange(-60, 70, 10)
        barra = plt.cm.RdBu_r
        ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
        im = ax.contourf(lon, lat, var, clevs, transform=ccrs.PlateCarree(),
                         cmap=barra, extend='both', vmin=-60, vmax=60)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        ax.coastlines()
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
        ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        plt.colorbar(im, orientation='horizontal')
        plt.title(title)
        plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
                    papertype='A4')
        plt.clf()
        plt.cla()
        plt.close()

def PlotEnsoComposites(var_pos, var_neg, lat, lon, title, filename):
        proj = ccrs.PlateCarree(central_longitude=180)
        fig = plt.figure(1, (10, 6.7), 300)
        ax = plt.subplot(2, 1, 1, projection=proj)
        clevs = np.arange(-60, 70, 10)
        barra = plt.cm.RdBu_r
        ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
        im = ax.contourf(lon, lat, var_pos, clevs, transform=ccrs.PlateCarree(),
                         cmap=barra, extend='both', vmin=-60, vmax=60)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        ax.coastlines()
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
        ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        plt.title('+NINIO3.4 YEARS')

        ax1 = plt.subplot(2, 1, 2, projection=proj)
        barra = plt.cm.RdBu_r
        ax1.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
        im = ax1.contourf(lon, lat, var_neg, clevs, transform=ccrs.PlateCarree(),
                         cmap=barra, extend='both', vmin=-60, vmax=60)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        ax1.coastlines()
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
        ax1.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax1.xaxis.set_major_formatter(lon_formatter)
        ax1.yaxis.set_major_formatter(lat_formatter)
        plt.title('-NINIO3.4 YEARS')
        plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
        fig.subplots_adjust(right=0.8)
        fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
        cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
                    papertype='A4')
        plt.clf()
        plt.cla()
        plt.close()

def PlotEnsoCompositesPoV(var1, var2, var3, var4, var5, var6, lat, lon, title, filename):
	vars = [var1, var2, var3, var4, var5, var6]
	tit = ['Ninio - Normal : All PoV', 'Ninia - Normal : All PoV',
		'Ninio - Normal : Weak PoV', 'Ninia - Normal : Weak PoV',
		'Ninio - Normal : Strong PoV', 'Ninia - Normal: Strong PoV']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
			 	cmap=barra, extend='both', vmin=-60, vmax=60)
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		ax.coastlines()
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()
def PlotKsEnsoCompositesPoV(var1, var2, var3, var4, var5, var6, lat, lon, title, filename):
	vars = [var1, var2, var3, var4, var5, var6]
	tit = ['Ninio - Normal : All PoV', 'Ninia - Normal : All PoV',
		'Ninio - Normal : Weak PoV', 'Ninia - Normal : Weak PoV',
		'Ninio - Normal : Strong PoV', 'Ninia - Normal: Strong PoV']
	proj = ccrs.PlateCarree(central_longitude=180)
	clevs = np.arange(1, 6, 1)
	barra = plt.cm.OrRd
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
			 	cmap=barra, extend='both', vmin=clevs[0], vmax=clevs[-1])
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		ax.coastlines()
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

def PlotPoVCompositesENSO(var1, var2, var3, var4, var5, var6, lat, lon, title, filename):
	vars = [var1, var2, var3, var4, var5, var6]
	tit = ['Strong SPoV - Normal SPov : All ', 'Weak SPoV - Normal SPov : All',
		'Strong SPoV - Normal SPov : Ninio ', 'Weak SPoV - Normal SPov : Ninio',
		'Strong SPoV - Normal SPov : Ninia ', 'Weak SPoV - Normal SPov : Ninia']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
			 	cmap=barra, extend='both', vmin=-60, vmax=60)
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		ax.coastlines()
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

def PlotPoVCompositesDiffENSO(var1, var2, var3, lat, lon, title, filename):
	vars = [var1, var2, var3]
	tit = ['All', 'El Ninio', 'La Ninia']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(3):
		ax = plt.subplot(3, 1, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
			 	cmap=barra, extend='both', vmin=-60, vmax=60)
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		ax.coastlines()
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

def PlotEnsoCompositesSAM(var1, var2, var3, var4, var5, var6, lat, lon, title, filename):
	vars = [var1, var2, var3, var4, var5, var6]
	tit = ['Ninio - Normal : All SAM', 'Ninia - Normal : All SAM',
		'Ninio - Normal : Weak SAM', 'Ninia - Normal : Weak SAM',
		'Ninio - Normal : Strong SAM', 'Ninia - Normal: Strong SAM']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, vars[i], clevs, transform=ccrs.PlateCarree(),
			 	cmap=barra, extend='both', vmin=-60, vmax=60)
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		ax.coastlines()
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()

def PlotEnsoCompositesPoVZPF(var, lat, lon, title, filename):
	var_keys = [*var]
	print(var_keys)
	tit = ['Ninio - Normal : All PoV', 'Ninia - Normal : All PoV',
		'Ninio - Normal : Weak PoV', 'Ninia - Normal : Weak PoV',
		'Ninio - Normal : Strong PoV', 'Ninia - Normal: Strong PoV']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(6):
		ax = plt.subplot(3, 2, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, var[var_keys[i * 3]], clevs, transform=ccrs.PlateCarree(),
			 	cmap=barra, extend='both', vmin=-60, vmax=60)
		lat1 = lat[lat < 0]
		u = var[var_keys[i * 3 + 1]][lat1 < -18, :]
		v = var[var_keys[i * 3 + 2]][lat1 < -18, :]
		lat1 = lat1[lat1 < -18]
		M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.07
		#mask array
		u_mask = ma.array(u,mask = M)
		v_mask = ma.array(v,mask = M)
		ax.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6],
			  v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(),
			  angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3,
			  headaxislength=2.8)
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		ax.coastlines()
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()
def PlotPoVCompositesDiffENSOZPF(var, lat, lon, title, filename):
	var_keys = [*var]
	#print(var_keys)
	tit = ['All', 'El Ninio', 'La Ninia']
	proj = ccrs.PlateCarree(central_longitude=180)
	fig = plt.figure(1, (10, 6.7), 300)
	for i in range(3):
		ax = plt.subplot(3, 1, i + 1, projection=proj)
		clevs = np.arange(-60, 70, 10)
		barra = plt.cm.RdBu_r
		ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
		im = ax.contourf(lon, lat, var[var_keys[i * 3]], clevs, transform=ccrs.PlateCarree(),
			 	cmap=barra, extend='both', vmin=-60, vmax=60)
		lat1 = lat[lat < 0]
		u = var[var_keys[i * 3 + 1]][lat1 < -18, :]
		v = var[var_keys[i * 3 + 2]][lat1 < -18, :]
		lat1 = lat1[lat1 < -18]
		M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.07
		#mask array
		u_mask = ma.array(u,mask = M)
		v_mask = ma.array(v,mask = M)
		ax.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6],
			  v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(),
			  angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3,
			  headaxislength=2.8)
		barra.set_under(barra(0))
		barra.set_over(barra(barra.N-1))
		ax.coastlines()
		ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
		ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
		lon_formatter = LongitudeFormatter(zero_direction_label=True)
		lat_formatter = LatitudeFormatter()
		ax.xaxis.set_major_formatter(lon_formatter)
		ax.yaxis.set_major_formatter(lat_formatter)
		plt.title(tit[i])
	plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
	fig.subplots_adjust(right=0.8)
	fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
	cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
	fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
	plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
		    papertype='A4')
	plt.clf()
	plt.cla()
	plt.close()




def PlotPVComposites(var_pos, var_neg, lat, lon, title, filename):
        proj = ccrs.PlateCarree(central_longitude=180)
        fig = plt.figure(1, (10, 6.7), 300)
        ax = plt.subplot(2, 1, 1, projection=proj)
        clevs = np.arange(-60, 70, 10)
        barra = plt.cm.RdBu_r
        ax.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
        im = ax.contourf(lon, lat, var_pos, clevs, transform=ccrs.PlateCarree(),
                         cmap=barra, extend='both', vmin=-60, vmax=60)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        ax.coastlines()
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
        ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        plt.title('STRONG PV YEARS')

        ax1 = plt.subplot(2, 1, 2, projection=proj)
        barra = plt.cm.RdBu_r
        ax1.set_extent([0, 359, -90, 20], crs=ccrs.PlateCarree())
        im = ax1.contourf(lon, lat, var_neg, clevs, transform=ccrs.PlateCarree(),
                         cmap=barra, extend='both', vmin=-60, vmax=60)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        ax1.coastlines()
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
        ax1.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax1.xaxis.set_major_formatter(lon_formatter)
        ax1.yaxis.set_major_formatter(lat_formatter)
        plt.title('WEAK PV YEARS')
        plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
        fig.subplots_adjust(right=0.8)
        fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
        cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
                    papertype='A4')
        plt.clf()
        plt.cla()
        plt.close()

def PlotCompPlumbDiff(var, u, v, lat, lon, title, filename):
        proj = ccrs.PlateCarree(central_longitude=180)
        fig = plt.figure(1, (10, 3.7), 300)
        ax = plt.subplot(projection=proj)
        clevs = np.arange(-60, 70, 10)
        barra = plt.cm.RdBu_r
        ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
        im = ax.contourf(lon, lat, var, clevs, transform=ccrs.PlateCarree(),
                         cmap=barra, extend='both', vmin=-60, vmax=60)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        lat1 = lat[lat < 0]
        u = u[lat1 < -18, :]
        v = v[lat1 < -18, :]
        lat1 = lat[lat < -18]
        Q50=np.percentile(np.sqrt(np.add(np.power(u,2),np.power(v,2))), 1) 
        M = np.sqrt(np.add(np.power(u,2),np.power(v,2))) < 0.08
        #mask array
        u_mask = ma.array(u,mask = M)
        v_mask = ma.array(v,mask = M)
        ax.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6], v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(), angles='xy', headwidth=2.5, headlength=1.8, width=1.5e-3, headaxislength=2.5)
        ax.coastlines()
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
        ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        plt.colorbar(im, orientation='horizontal')
        plt.title(title)
        plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
                    papertype='A4')
        plt.clf()
        plt.cla()
        plt.close()

def PlotCompositesPlumb(var_pos, var_neg, u_pos, u_neg, v_pos, v_neg, lat,
			lon, title, filename):
        proj = ccrs.PlateCarree(central_longitude=180)
        fig = plt.figure(1, (10, 6.7), 300)
        ax = plt.subplot(2, 1, 1, projection=proj)
        clevs = np.arange(-60, 70, 10)
        barra = plt.cm.RdBu_r
        ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
        im = ax.contourf(lon, lat, var_pos, clevs, transform=ccrs.PlateCarree(),
                         cmap=barra, extend='both', vmin=-60, vmax=60)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        lat1 = lat[lat < 0]
        u_pos = u_pos[lat1 < -18, :]
        v_pos = v_pos[lat1 < -18, :]
        lat1 = lat[lat < -18]
        Q50=np.percentile(np.sqrt(np.add(np.power(u_pos,2),np.power(v_pos,2))), 1) 
        M = np.sqrt(np.add(np.power(u_pos,2),np.power(v_pos,2))) < 0.07
        #mask array
        u_mask = ma.array(u_pos,mask = M)
        v_mask = ma.array(v_pos,mask = M)
        ax.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6], v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(), angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
        ax.coastlines()
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
        ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        plt.title('+NINIO3.4 YEARS')

        ax1 = plt.subplot(2, 1, 2, projection=proj)
        barra = plt.cm.RdBu_r
        ax1.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
        im = ax1.contourf(lon, lat, var_neg, clevs, transform=ccrs.PlateCarree(),
                         cmap=barra, extend='both', vmin=-60, vmax=60)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        lat1 = lat[lat < 0]
        u_neg = u_neg[lat1 < -18, :]
        v_neg = v_neg[lat1 < -18, :]
        lat1 = lat[lat < -18]
        Q50=np.percentile(np.sqrt(np.add(np.power(u_neg,2),np.power(v_neg,2))), 1) 
        M = np.sqrt(np.add(np.power(u_neg,2),np.power(v_neg,2))) < 0.07
        #mask array
        u_mask = ma.array(u_neg,mask = M)
        v_mask = ma.array(v_neg,mask = M)
        ax1.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6], v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(), angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
        ax1.coastlines()
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
        ax1.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax1.xaxis.set_major_formatter(lon_formatter)
        ax1.yaxis.set_major_formatter(lat_formatter)
        plt.title('-NINIO3.4 YEARS')
        plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
        fig.subplots_adjust(right=0.8)
        fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
        cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
                    papertype='A4')
        plt.clf()
        plt.cla()
        plt.close()

def PlotCompositesPlumbPV(var_pos, var_neg, u_pos, u_neg, v_pos, v_neg, lat,
			lon, title, filename):
        proj = ccrs.PlateCarree(central_longitude=180)
        fig = plt.figure(1, (10, 6.7), 300)
        ax = plt.subplot(2, 1, 1, projection=proj)
        clevs = np.arange(-60, 70, 10)
        barra = plt.cm.RdBu_r
        ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
        im = ax.contourf(lon, lat, var_pos, clevs, transform=ccrs.PlateCarree(),
                         cmap=barra, extend='both', vmin=-60, vmax=60)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        lat1 = lat[lat < 0]
        u_pos = u_pos[lat1 < -18, :]
        v_pos = v_pos[lat1 < -18, :]
        lat1 = lat[lat < -18]
        Q50=np.percentile(np.sqrt(np.add(np.power(u_pos,2),np.power(v_pos,2))), 1) 
        M = np.sqrt(np.add(np.power(u_pos,2),np.power(v_pos,2))) < 0.07
        #mask array
        u_mask = ma.array(u_pos,mask = M)
        v_mask = ma.array(v_pos,mask = M)
        ax.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6], v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(), angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
        ax.coastlines()
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
        ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        plt.title('Strong PV YEARS')

        ax1 = plt.subplot(2, 1, 2, projection=proj)
        barra = plt.cm.RdBu_r
        ax1.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
        im = ax1.contourf(lon, lat, var_neg, clevs, transform=ccrs.PlateCarree(),
                         cmap=barra, extend='both', vmin=-60, vmax=60)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        lat1 = lat[lat < 0]
        u_neg = u_neg[lat1 < -18, :]
        v_neg = v_neg[lat1 < -18, :]
        lat1 = lat[lat < -18]
        Q50=np.percentile(np.sqrt(np.add(np.power(u_neg,2),np.power(v_neg,2))), 1) 
        M = np.sqrt(np.add(np.power(u_neg,2),np.power(v_neg,2))) < 0.07
        #mask array
        u_mask = ma.array(u_neg,mask = M)
        v_mask = ma.array(v_neg,mask = M)
        ax1.quiver(lon[2:-1:6], lat1[2:-1:6], u_mask[2:-1:6, 2:-1:6], v_mask[2:-1:6, 2:-1:6], transform=ccrs.PlateCarree(), angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
        ax1.coastlines()
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
        ax1.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax1.xaxis.set_major_formatter(lon_formatter)
        ax1.yaxis.set_major_formatter(lat_formatter)
        plt.title('Weak PV YEARS')
        plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
        fig.subplots_adjust(right=0.8)
        fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
        cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
                    papertype='A4')
        plt.clf()
        plt.cla()
        plt.close()

def PlotCompDivPlumbDiff(var,u, v, lat, lon, title, filename):
        proj = ccrs.PlateCarree(central_longitude=180)
        fig = plt.figure(1, (10, 3.7), 300)
        ax = plt.subplot(projection=proj)
        clevs = np.arange(-2.4, 4, 1.6)
        barra = plt.cm.RdBu_r
        ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
        im = ax.contourf(lon, lat, var / 10, clevs,
			transform=ccrs.PlateCarree(), cmap=barra, extend='both',
			vmin=-2.4, vmax=2.4)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        M = np.sqrt(np.add(np.power(u, 2),np.power(v, 2))) < 0.5
        #mask array
        u_mask = ma.array(u, mask = M)
        v_mask = ma.array(v, mask = M)
        ax.quiver(lon[0:-1:8], lat[0:-1:8], u_mask[0:-1:8, 0:-1:8], v_mask[0:-1:8, 0:-1:8], transform=ccrs.PlateCarree(), angles='xy', headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
        ax.coastlines()
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
        ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        plt.colorbar(im, orientation='horizontal', shrink=0.25)
        plt.title(title)
        plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
                    papertype='A4')
        plt.clf()
        plt.cla()
        plt.close()

def PlotCompositesDivPlumb(var_pos, var_neg, u_pos, u_neg, v_pos, v_neg, lat, lon, title, filename):
        proj = ccrs.PlateCarree(central_longitude=180)
        fig = plt.figure(1, (10, 6.7), 300)
        ax = plt.subplot(2, 1, 1, projection=proj)
        clevs = np.arange(-1.2, 1.6, 0.8)
        barra = plt.cm.RdBu_r
        ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
        im = ax.contourf(lon, lat, var_pos / 10, clevs,
			transform=ccrs.PlateCarree(), cmap=barra, extend='both',
			vmin=-1.2, vmax=1.2)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        M = np.sqrt(np.add(np.power(u_pos, 2),np.power(v_pos, 2))) < 0.1
        #mask array
        u_mask = ma.array(u_pos, mask = M)
        v_mask = ma.array(v_pos, mask = M)
        ax.quiver(lon[0:-1:6], lat[0:-1:6], u_mask[0:-1:6, 0:-1:6], v_mask[0:-1:6, 0:-1:6], transform=ccrs.PlateCarree(), angles='xy')#, headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
        ax.coastlines()
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
        ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        plt.title('+NINIO3.4 YEARS')
        ax1 = plt.subplot(2, 1, 2, projection=proj)
        ax1.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
        im = ax1.contourf(lon, lat, var_neg / 10, clevs,
			 transform=ccrs.PlateCarree(),
                         cmap=barra, extend='both', vmin=-1.2, vmax=1.2)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        M = np.sqrt(np.add(np.power(u_neg, 2),np.power(v_neg, 2))) < 0.3
        #mask array
        u_mask = ma.array(u_neg, mask = M)
        v_mask = ma.array(v_neg, mask = M)
        ax1.quiver(lon[0:-1:6], lat[0:-1:6], u_mask[0:-1:6, 0:-1:6], v_mask[0:-1:6, 0:-1:6], transform=ccrs.PlateCarree(), angles='xy')#, headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
        ax1.coastlines()
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
        ax1.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax1.xaxis.set_major_formatter(lon_formatter)
        ax1.yaxis.set_major_formatter(lat_formatter)
        plt.title('-NINIO3.4 YEARS')
        plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
        fig.subplots_adjust(right=0.8)
        fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
        cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
                    papertype='A4')
        plt.clf()
        plt.cla()
        plt.close()
def PlotCompositesDivPlumbPV(var_pos, var_neg, u_pos, u_neg, v_pos, v_neg, lat, lon, title, filename):
        proj = ccrs.PlateCarree(central_longitude=180)
        fig = plt.figure(1, (10, 6.7), 300)
        ax = plt.subplot(2, 1, 1, projection=proj)
        clevs = np.arange(-1.2, 1.6, 0.8)
        barra = plt.cm.RdBu_r
        ax.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
        im = ax.contourf(lon, lat, var_pos / 10, clevs,
			transform=ccrs.PlateCarree(), cmap=barra, extend='both',
			vmin=-1.2, vmax=1.2)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        M = np.sqrt(np.add(np.power(u_pos, 2),np.power(v_pos, 2))) < 0.3
        #mask array
        u_mask = ma.array(u_pos, mask = M)
        v_mask = ma.array(v_pos, mask = M)
        ax.quiver(lon[0:-1:6], lat[0:-1:6], u_mask[0:-1:6, 0:-1:6], v_mask[0:-1:6, 0:-1:6], transform=ccrs.PlateCarree(), angles='xy')#, headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
        ax.coastlines()
        ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
        ax.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
        plt.title('Strong PV YEARS')
        ax1 = plt.subplot(2, 1, 2, projection=proj)
        ax1.set_extent([0, 359, -90, 0], crs=ccrs.PlateCarree())
        im = ax1.contourf(lon, lat, var_neg / 10, clevs,
			 transform=ccrs.PlateCarree(),
                         cmap=barra, extend='both', vmin=-1.2, vmax=1.2)
        barra.set_under(barra(0))
        barra.set_over(barra(barra.N-1))
        M = np.sqrt(np.add(np.power(u_neg, 2),np.power(v_neg, 2))) < 0.3
        #mask array
        u_mask = ma.array(u_neg, mask = M)
        v_mask = ma.array(v_neg, mask = M)
        ax1.quiver(lon[0:-1:6], lat[0:-1:6], u_mask[0:-1:6, 0:-1:6], v_mask[0:-1:6, 0:-1:6], transform=ccrs.PlateCarree(), angles='xy')#, headwidth=2.8, headlength=1.8, width=1.5e-3, headaxislength=2.8)
        ax1.coastlines()
        ax1.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5) #countries
        ax1.gridlines(crs=proj, linewidth=0.3, linestyle='-')
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax1.xaxis.set_major_formatter(lon_formatter)
        ax1.yaxis.set_major_formatter(lat_formatter)
        plt.title('Weak PV YEARS')
        plt.suptitle(title, fontsize=12, x=0.47, y=0.9)
        fig.subplots_adjust(right=0.8)
        fig.subplots_adjust(bottom=0.17, top=0.82, hspace=0.2, wspace=0.05)
        cbar_ax = fig.add_axes([0.33, 0.1, 0.25, 0.05])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        plt.savefig(filename, dpi=300, bbox_inches='tight', orientation='landscape',
                    papertype='A4')
        plt.clf()
        plt.cla()
        plt.close()
