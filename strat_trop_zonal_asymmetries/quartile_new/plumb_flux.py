#computation of plum fluxes. Script adapted from Kazuaki Nishii and Hisashi Nakamura
import numpy as np
import math
import xarray as xr
from metpy import calc
        
def c_diff(arr, h, dim, cyclic = False):  
	#compute derivate of array variable respect to h associated to dim
	#adapted from kuchaale script
	ndim = arr.ndim
	lst = [i for i in range(ndim)]
	lst[dim], lst[0] = lst[0], lst[dim]
	rank = lst
	arr = np.transpose(arr, tuple(rank))
	if ndim == 3:
		shp = (arr.shape[0]-2,1,1)
	elif ndim == 4:
		shp = (arr.shape[0]-2,1,1,1)
	else:
		shp = (arr.shape[0]-2,1)
	d_arr = np.copy(arr)
	if not cyclic:
		d_arr[0, ...] = (arr[1, ...] - arr[0, ...]) / (h[1] - h[0])
		d_arr[-1, ...] = (arr[-1, ...] - arr[-2, ...]) / (h[-1] - h[-2])
		d_arr[1:-1, ...] = (arr[2:, ...] - arr[0:-2, ...]) / np.reshape(h[2:] - h[0:-2], shp)
	elif cyclic:
		d_arr[0, ...] = (arr[1, ...] - arr[-1, ...]) / (h[1] - h[-1])
		d_arr[-1, ...] = (arr[0, ...] - arr[-2, ...]) / (h[0] - h[-2])
		d_arr[1:-1, ...] = (arr[2:, ...] - arr[0:-2, ...]) / np.reshape(h[2:] - h[0:-2], shp)
	d_arr = np.transpose(d_arr, tuple(rank))
	return d_arr

def ComputePlumbFluxes(uclm, vclm, z, zclm, lat, lon, lev):
	#imput: model u, v, z with multiple realizations
	#restringe domain to latitudes south to 0 (to avoid problems with sin (0°))
	z = z[:, lat < 0, :]
	zclm = zclm[lat < 0, :]
	uclm = uclm[lat < 0, :]
	vclm = vclm[lat < 0, :]
	lat = lat[lat < 0]
	#compute climatologies
	[nrealiz, nlats, nlons] = np.shape(z) #get dimensions
	#earth radius
	a = 6400000
	coslat = np.cos(lat * 3.14 / 180)
	sinlat = np.sin(lat * 3.14 / 180)
	#Coriolis parameter
	f = 2 * 7.29e-5 * sinlat
	#gravity
	g = 9.8
	# basic state (climatology): uclm vclm zclm
	# anomalies zaa
	zaa = z - zclm
	# QG stream function
#	psiaa = g / np.transpose(np.tile(f, (nlons, 1))) * zaa
	psiaa = g / np.transpose(np.tile(f, (nrealiz, nlons, 1)), (0, 2, 1)) * zaa

	# magnitude of basic state wind speed
	magU = np.sqrt(np.add(np.power(uclm, 2), np.power(vclm, 2)))
	#psi derivatives
	dpsidlon = c_diff(psiaa, lon * 3.14 /180, 2)
	ddpsidlonlon = c_diff(dpsidlon, lon * 3.14, 2)
	dpsidlat = c_diff(psiaa, lat * 3.14 / 180, 1)
	ddpsidlatlat = c_diff(dpsidlat, lat * 3.14 / 180, 1)
	ddpsidlatlon = c_diff(dpsidlat, lon * 3.14 /180, 2)
	termxu = dpsidlon * dpsidlon - psiaa * ddpsidlonlon
	termxv = dpsidlon * dpsidlat - ddpsidlatlon * psiaa
	termyv = dpsidlat * dpsidlat - psiaa * ddpsidlatlat
	# "p" is normalized by 1000hPa
	coefcos = np.transpose(np.tile(coslat, (nrealiz, nlons, 1)), (0, 2, 1))
	coeff = coefcos * (lev / 100000) / (2 * magU)
	#coeff = 1 / (2 * magU)
	#x-component
	px = coeff / (a * a * coefcos) * (uclm * termxu / coefcos + vclm * termxv)
	#px = coeff * (uclm * termxu + vclm * termxv)
	#y-component
	py = coeff / (a * a) * ( uclm / coefcos * termxv + vclm * termyv)
	#py = coeff * ( uclm * termxv + vclm * termyv)

	return px, py, lat
def ComputePlumbFluxesDecomp(uclm, vclm, z, zclm, lat, lon, lev):
	#input: model climatology of u,v, z
	#z for the event to analyze
	#restringe domain to latitudes south to 0 (to avoid problems with sin (0°))
	z = z[:, lat < 0, :]
	zclm = zclm[lat < 0, :]
	uclm = uclm[lat < 0, :]
	vclm = vclm[lat < 0, :]
	lat = lat[lat < 0]
	[nrealiz, nlats, nlons] = np.shape(z) #get dimensions
	#earth radius
	a = 6400000
	coslat = np.cos(lat * 3.14 / 180)
	sinlat = np.sin(lat * 3.14 / 180)
	#Coriolis parameter
	f = 2 * 7.29e-5 * sinlat
	#gravity
	g = 9.8
	# magnitude of basic state wind speed
	magU = np.sqrt(np.add(np.power(uclm, 2), np.power(vclm, 2)))
	#QG stream function
	psia_clim =  g / np.transpose(np.tile(f, (nlons, 1))) * zclm
	psia = g / np.transpose(np.tile(f, (nrealiz, nlons, 1)), (0, 2, 1)) * z
	psia_m = np.mean(psia, axis=0)

	#anomalies with respecto to the ensemble mean
	psiaa = psia - psia_m
	psia_m = psia_m - psia_clim
	#psia_clim = psia_clim - np.transpose(np.tile(np.mean(psia_clim, axis=1), (nlons, 1)))
	#psi derivatives
	dpsiamdlon = c_diff(psia_m, lon * 3.14 /180, 1)
	ddpsiamdlonlon = c_diff(dpsiamdlon, lon * 3.14, 1)
	dpsiamdlat = c_diff(psia_m, lat * 3.14 / 180, 0)
	ddpsiamdlatlat = c_diff(dpsiamdlat, lat * 3.14 / 180, 0)
	ddpsiamdlatlon = c_diff(dpsiamdlat, lon * 3.14 /180, 1)
	
	dpsiacdlon = c_diff(psia_clim, lon * 3.14 /180, 1)
	ddpsiacdlonlon = c_diff(dpsiacdlon, lon * 3.14, 1)
	dpsiacdlat = c_diff(psia_clim, lat * 3.14 / 180, 0)
	ddpsiacdlatlat = c_diff(dpsiacdlat, lat * 3.14 / 180, 0)
	ddpsiacdlatlon = c_diff(dpsiacdlat, lon * 3.14 /180, 1)

	dpsidlon = c_diff(psiaa, lon * 3.14 /180, 2)
	ddpsidlonlon = c_diff(dpsidlon, lon * 3.14, 2)
	dpsidlat = c_diff(psiaa, lat * 3.14 / 180, 1)
	ddpsidlatlat = c_diff(dpsidlat, lat * 3.14 / 180, 1)
	ddpsidlatlon = c_diff(dpsidlat, lon * 3.14 /180, 2)
	#linear terms
	termx1l = 2 * dpsiacdlon * dpsiamdlon - psia_clim * ddpsiamdlonlon - psia_m * ddpsiacdlonlon
	termxyl = dpsiacdlon * dpsiamdlat + dpsiamdlon * dpsiacdlat -\
		  psia_clim * ddpsiamdlatlon - psia_m * ddpsiacdlatlon
	termy1l =  2 * dpsiacdlat * dpsiamdlat - psia_clim * ddpsiamdlatlat - psia_m * ddpsiacdlatlat

	# "p" is normalized by 1000hPa
	coefcos = np.transpose(np.tile(coslat, (nlons, 1)))
	coeff = coefcos * (lev / 100000) / (2 * magU)
	#coeff = 1 / (2 * magU)
	#x-component
	pxl = coeff / (a * a * coefcos) * (uclm * termx1l / coefcos + vclm * termxyl)
	#px = coeff * (uclm * termxu + vclm * termxv)
	#y-component
	pyl = coeff / (a * a) * ( uclm / coefcos * termxyl + vclm * termy1l)
	#py = coeff * ( uclm * termxv + vclm * termyv)

	#non-linear terms
	termx1nl = dpsiamdlon ** 2 - psia_m * ddpsiamdlonlon
	termxynl = dpsiamdlat * dpsiamdlon - psia_m * ddpsiamdlatlon
	termy1nl = dpsiamdlat ** 2 - psia_m * ddpsiamdlatlat
	#x-component
	pxnl = coeff / (a * a * coefcos) * (uclm * termx1nl / coefcos + vclm * termxynl)
	#px = coeff * (uclm * termxu + vclm * termxv)
	#y-component
	pynl = coeff / (a * a) * ( uclm / coefcos * termxynl + vclm * termy1nl)
	#py = coeff * ( uclm * termxv + vclm * termyv)
	#compute divergence
	dx, dy = calc.lat_lon_grid_deltas(lon, lat)
	divl = calc.divergence(pxl, pyl, dx, dy)
	divnl = calc.divergence(pxnl, pynl, dx, dy)
	pxem = pxl + pxnl
	pyem = pyl + pynl
	divem = calc.divergence(pxem, pyem, dx, dy)

	var = {'pxl': pxl, 'pyl': pyl, 'pxnl': pxnl, 'pynl': pynl, 'divl': divl, 'divnl': divnl,
		'pxem': pxem, 'pyem': pyem, 'divem': divem,  'lat': lat}
	return var

def ComputePlumbFluxes3D(uclm, vclm, z, zclm, lat, lon, lev):
	#input: model climatology of u,v, z
	#z [nrealiz, nlevs, nlats, nlons]
	#restringe domain to latitudes south to 0 (to avoid problems with sin (0°))
	z = z[:, :, lat < 0, :]
	zclm = zclm[:, lat < 0, :]
	uclm = uclm[:, lat < 0, :]
	vclm = vclm[:, lat < 0, :]
	lat = lat[:, lat < 0]
	[nrealiz, nlevs, nlats, nlons] = np.shape(z) #get dimensions
	#earth radius
	a = 6400000
	coslat = np.cos(lat * 3.14 / 180)
	sinlat = np.sin(lat * 3.14 / 180)
	p = lev/100000
	#Coriolis parameter
	f = 2 * 7.29e-5 * sinlat
	#gravity
	g = 9.8
	#buoyanncy
	N2 = 4e-4
	f0 = 2 * 7.29e-5 * np.sin(-45 * 3.14 / 180)
	H = 6.4e3
	z = -H * np.log(p)
	# magnitude of basic state wind speed
	magU = np.tile(np.sqrt(np.add(np.power(uclm, 2), np.power(vclm, 2))), (nrealiz, 1, 1, 1))
	#QG stream function
	psia = g / np.transpose(np.tile(f, (nrealiz, nlevs, nlons, 1)), (0, 1, 3, 2)) * z
	#psi derivatives
	dpsiadlon = c_diff(psia, lon * 3.14 /180, 3)
	ddpsiadlonlon = c_diff(dpsiadlon, lon * 3.14, 3)
	dpsiadlat = c_diff(psia, lat * 3.14 / 180, 2)
	ddpsiadlatlat = c_diff(dpsiadlat, lat * 3.14 / 180, 2)
	ddpsiadlatlon = c_diff(dpsiadlat, lon * 3.14 /180, 3)
	dpsiadz = c_diff(psia, z, 1)
	ddpsiadzlat = c_diff(dpsiadz, lat * 3.14 / 180, 2)
	ddpsiadzlon = c_diff(dpsiadz, lon * 3.14 / 180, 3)
	#linear terms
	termx1 =  dpsiadlon ** 2 - psia * ddpsiadlonlon 
	termxy = dpsiadlon * dpsiadlat - psia * ddpsiadlatlon 
	termy1 = dpsiadlat **2 - psia * ddpsiamdlatlat
	termzx = dspiadlon * dpsiadz - psia * ddpsiadzlon
	termzy = dspiadlat * dpsiadz - psia * ddpsiadzlat
	# "p" is normalized by 1000hPa
	coefcos = np.transpose(np.tile(coslat, (nrealiz, nlevs, nlons, 1)), (0, 1, 3, 2))
	p = np.transpose(np.tile(p, (nrealiz, nlats, nlons, 1)), (0, 3, 1, 2))
	coeff = coefcos * p / (2 * magU)
	Uclm = np.tile(uclm, (nrealiz, 1, 1, 1))
	Vclm = np.tile(vclm, (nrealiz, 1, 1, 1))
	#coeff = 1 / (2 * magU)
	#x-component
	px = coeff / (a * a * coefcos) * (Uclm * termx1l / coefcos + Vclm * termxyl)
	#px = coeff * (uclm * termxu + vclm * termxv)
	#y-component
	py = coeff / (a * a) * ( Uclm / coefcos * termxy + Vclm * termy1)
	#py = coeff * ( uclm * termxv + vclm * termyv)
	pz = coeff * fo / N2 * (Uclm / (a * coefcos) * termzx + Vclm / a * termzy)

	var = {'px': px, 'py': py, 'pz': pz, 'lat': lat}
	return var

