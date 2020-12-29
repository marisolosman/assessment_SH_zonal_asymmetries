import numpy as np

def doharmonic(data, ndt, kmax):
	#function [exp_var y c]=doharmonic(data,ndt,kmax)
	#data=Serie de tiempo a analizar. Datos deben ser pares.
	#ndt=cantidad de unidades en un delta t (1 día)
	#kmax=mayor número de onda que me interesa resolver (el m�ximo posible es N/2)
	data = np.transpose(data)
	#cuenta la cantidad de datos
	N = np.size(data)
	#periodo fundamental
	P = N * ndt
	#varianza de la serie original
	ovar = np.var(data) #el 1 indica que normaliza por n y no por n-1
	a = np.empty([kmax])
	b = np.empty([kmax])
	y = np.empty([N, kmax])
	c = np.empty([kmax])
	exp_var = np.empty([kmax])
	#se calculan los coeficientes Ai para los senos y Bi para los cosenos
	for i in np.arange(0, kmax):
		sum_sin = 0
		sum_cos = 0
		for t in np.arange(0, N):
			sum_sin += data[t] * np.sin((2 * np.pi * (i + 1) * (t + 1)) / P)
			sum_cos += data[t] * np.cos((2 * np.pi * (i + 1) * (t + 1)) / P)
		a[i] = (2 / N) * sum_sin
		b[i] = (2 / N) * sum_cos
	#se calculan las amplitudes de los armonicos
	for i in np.arange(0, kmax-1):
		c[i] = np.sqrt(np.power(a[i], 2) + np.power(b[i], 2))
		exp_var[i] = (np.power(c[i], 2) / (2 * ovar)) * 100

	c[kmax - 1] = b[kmax - 1] / 2
	exp_var[kmax - 1] = (np.power(c[kmax - 1], 2) / ovar) * 100
    
	#Reconstruyo los armónicos
	for k in np.arange(0, kmax):
		for j in np.arange(0, N):
			y[j, k] = a[k] * np.sin((2 * np.pi / P) * (k + 1) * (j + 1)) + b[k] *\
				  np.cos((2 * np.pi / P) * (k + 1) * (j + 1))
	return exp_var, y, c
