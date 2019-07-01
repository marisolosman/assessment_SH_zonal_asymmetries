import numpy as np
from numpy import linalg as LA

def eofdata(x, modes):
#function [evar,E,Z,D,Zs]=eofdata(x,modes)

#==================================================                                                  #   EOF/PC analysis from the singular value decomposition of the     
#   standardized data matrix     
#                                                                                                                                                             

#    SOURCES:

#   GREEN,PE,1978. Analysing Multivariate Data. The Dryden Press. 
#                            Chapter 8

# HARTMANN, DL. Objective Analysis Class Notes (ATM 552)
#                           Chapter 4

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   X: standardized data matrix. Thought to have variables as rows (grid points),
# and samples as columns (times)
#   MODES: number of factors to retain
#   
	[vars, samps] = np.shape(x)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	corr_matrix = np.corrcoef(x)
	w, v = LA.eig(corr_matrix)
	#PCs
	u = np.matmul(np.transpose(v[:, 0:modes]), np.ma.anomalies(x, 1))
	return w[0:modes], v, u

