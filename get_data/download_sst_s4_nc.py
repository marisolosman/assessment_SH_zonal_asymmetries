#!/usr/bin/env python
from ecmwfapi import ECMWFService
import calendar
  
server = ECMWFService("mars")

#for year in range(2016, 2018):
#	#Select year for data extraction
#	print ('YEAR ',year)
#	
#	initmonth=8
#	initbdate="%s%02d01"%(year,initmonth)
#	print ("######### ERA-interim  #########")
#	print ('get data from ', initbdate)
#	print ("################################")		
#
#	server.execute({
#		'class'   : "od",
#		'step'    : "2904/to/5160/by/24",
#		'number'  : "0/to/50",   
#		'levtype' : "sfc",
#		'expver'  : "1",
#		'method'  : "1",
#		'origin'  : "ecmf",    
#		'date'    : "%s"%(initbdate), 
#		'time'    : "00:00:00",   
#		'type'    : "fc",
#		'param'   : "34.128", 
#		'stream'  : "mmsf",
#		'system'  : "4",
#		'grid'    : "128"},
#		 "/storage/silver/acrcc/vg140344/HGT_S4Hindcasts_0108_SST_24Hourly_%s.grib"%(initbdate))
initmonth=8
year=1992
initbdate="%s%02d01"%(year,initmonth)
print ("######### ERA-interim  #########")
print ('get data from ', initbdate)
print ("################################")		
server.execute({
	'class'   : "od",
	'step'    : "2904/to/5160/by/24",
	'number'  : "0/to/50",   
	'levtype' : "sfc",
	'expver'  : "1",
	'method'  : "1",
	'origin'  : "ecmf",    
	'date'    : "%s"%(initbdate), 
	'time'    : "00:00:00",   
	'type'    : "fc",
	'param'   : "34.128", 
	'stream'  : "mmsf",
	'system'  : "4",
	'grid'    : "128"},
	 "/storage/silver/acrcc/vg140344/SST_S4Hindcasts_0108_24Hourly_%s.grib"%(initbdate))
