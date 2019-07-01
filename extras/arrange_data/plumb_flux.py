#computation of plum fluxes. Script adapted from Kazuaki Nishii and Hisashi Nakamura
import numpy as np
import math
import xarray as xr
        
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

def ComputePlumbFluxes(u, v, z, lat, lon):
	#imput: model u, v, z with multiple realizations
	#restringe domain to latitudes south to 0 (to avoid problems with sin (0°))
	z = z[:, lat < 0, :]
	u = u[:, lat < 0, :]
	v = v[:, lat < 0, :]
	lat = lat[lat < 0]
	#compute climatologies
	uclm = np.mean(u, axis=0)
	vclm = np.mean(v, axis=0)
	zclm = np.mean(z, axis=0)
	[realiz, nlats, nlons] = np.shape(z) #get dimensions
	#earth radius
	a = 6400000
	coslat = np.cos(lat * 3.14 / 180)
	sinlat = np.sin(lat * 3.14 / 180)
	#Coriolis parameter
	f = 2 * 7.29e-5 * sinlat
	#gravity
	g = 9.8
	# unit [Pa]
	lev = 20000
	# basic state (climatology): uclm vclm zclm
	# anomalies zaa
	zaa = z - zclm
	# QG stream function
	psiaa = g / np.transpose(np.tile(f, (realiz, nlons, 1)),[0, 2, 1]) * zaa
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
	coefcos = np.transpose(np.tile(coslat, (realiz, nlons, 1)), [0, 2, 1])
	coeff = coefcos * (lev / 100000) / (2 * magU)
	#x-component
	px = coeff / (a * a * coefcos) * (uclm * termxu / coefcos + vclm * termxv)
	#y-component
	py = coeff / (a * a) * ( uclm / coefcos * termxv + vclm * termyv)
    
	return px, py, lat
