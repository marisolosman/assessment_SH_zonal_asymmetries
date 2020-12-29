#!/usr/bin/env python
from ecmwfapi import ECMWFService
import calendar
  
server = ECMWFService("mars")

for year in range(1981, 2018):
	#Select year for data extraction
	print ('YEAR ',year)	
	initmonth=8
	initbdate="%s%02d01"%(year,initmonth)
	print ("######### ERA-interim  #########")
	print ('get data from ', initbdate)
	print ("################################")		

	server.execute({
		'class'   : "od",
		'step'    : "0/to/5160/by/24",
		'number'  : "0/to/50",   
		'levtype' : "sfc",
		'expver'  : "1",
		'method'  : "1",
		'origin'  : "ecmf",    
		'date'    : "%s"%(initbdate), 
		'time'    : "00:00:00",   
		'type'    : "fc",
		'param'   : "31.128", 
		'stream'  : "mmsf",
		'system'  : "4",
		'grid'    : "128",},
		 "/storage/shared/glusterfs/acrcc/users/vg140344/SIF_S4Hindcasts_24Hourly_%s.grib"%(initbdate))

