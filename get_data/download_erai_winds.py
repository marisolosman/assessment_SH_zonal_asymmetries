#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
import numpy as np
for i in range(1981, 2019):
	output = '/storage/silver/acrcc/vg140344/ERAI/WINDS_erai_' + str(i) + ".grib"
	server = ECMWFDataServer()
	server.retrieve({
		"class": "ei",
		"dataset": "interim",
		"date": "".join("%s%02d01/"%(i, j) for j in range(1, 13))[0: -1],
		"expver": "1",
		"grid": "0.75/0.75",
		"levelist": "200/50",
		"levtype": "pl",
		"param": "131.128/132.128",
		"stream": "moda",
		"type": "an",
		"target": output,
		})

