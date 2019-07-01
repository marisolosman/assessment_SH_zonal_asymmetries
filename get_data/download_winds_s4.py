#!/usr/bin/env python
from ecmwfapi import ECMWFService
import calendar
  
server = ECMWFService("mars")

for year in range(1989, 2018):
	#Select year for data extraction
	print ('YEAR ',year)
	initmonth=11
	initbdate="%s%02d01"%(year,initmonth)
	print ("######### ERA-interim  #########")
	print ('get data from ', initbdate)
	print ("################################")		

	server.execute({
		'class'   : "od",
		'levelist': "850/200",
		'step'    : "0/to/2904/by/24",
		'number'  : "0/to/50",   
		'levtype' : "pl",
		'expver'  : "1",
		'method'  : "1",
		'origin'  : "ecmf",    
		'date'    : "%s"%(initbdate), 
		'time'    : "00:00:00",   
		'type'    : "fc",
		'param'   : "131/132", 
		'stream'  : "mmsf",
		'system'  : "4",
		'grid'    : "128"},
		 "/home/users/vg140344/datos/UV_S4Hindcasts_levels_24Hourly_%s.grib"%(initbdate))
