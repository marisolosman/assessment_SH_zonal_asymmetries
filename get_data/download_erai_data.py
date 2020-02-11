##!/usr/bin/env python
#from ecmwfapi import ECMWFDataServer
#import numpy as np
#for i in range(1981,2019):
#	output = '/storage/shared/glusterfs/acrcc/users/vg140344/sst_erai_' + str(i) + ".grib"
#	server = ECMWFDataServer()
#	server.retrieve({
#		"class": "ei",
#		"dataset": "interim",
#		"date": "".join("%s%02d01/"%(i, j) for j in range(1,13))[0:-1],
#		"expver": "1",
#		"grid": "0.75/0.75",
#		"levtype": "sfc",
#		"param": "34.128",
#		"stream": "moda",
#		"type": "an",
#		"target": output,
#		})
#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
import numpy as np
#i = 2013
#output = '/storage/shared/glusterfs/acrcc/users/vg140344/sst_erai_' + str(i) + ".grib"
#server = ECMWFDataServer()
#server.retrieve({
#	"class": "ei",
#	"dataset": "interim",
#	"date": "".join("%s%02d01/"%(i, j) for j in range(1,13))[0:-1],
#	"expver": "1",
#	"grid": "0.75/0.75",
#	"levtype": "sfc",
#	"param": "34.128",
#	"stream": "moda",
#	"type": "an",
#	"target": output,
#	})
for i in range(1981,2019):
	output = '/storage/shared/glusterfs/acrcc/users/vg140344/ERAI/HGT50_erai_' + str(i) + ".grib"
	server = ECMWFDataServer()
	server.retrieve({
		"class": "ei",
		"dataset": "interim",
		"date": "".join("%s%02d01/"%(i, j) for j in range(1,13))[0:-1],
		"expver": "1",
		"grid": "0.75/0.75",
		"levelist": "50",
		"levtype": "pl",
		"param": "129.128",
		"stream": "moda",
		"type": "an",
		"target": output,
		})

