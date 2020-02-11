import glob
import numpy as np

lista = glob.glob("/home/users/vg140344/decile/correlaciones/*")
correlations = []
for i in lista:
	f = np.load(i)
	correlations.append(f['corr'])
	f.close()

lista = []

correlations = np.array(correlations)
#compute median and 5 and 95 percentiles
mediana = np.median(correlations, axis=0)
intervals = np.percentile(correlations, [5, 95], axis=0)


	
