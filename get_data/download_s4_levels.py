#!/usr/bin/env python
from ecmwfapi import ECMWFService
import calendar
  
server = ECMWFService("mars")

for year in range(1981, 2018):
	#Select year for data extraction
	print ('YEAR ',year)
	initmonth = 8
	initbdate = "%s%02d01"%(year,initmonth)
	print ("######### ERA-interim  #########")
	print ('get data from ', initbdate)
	print ("################################")		

	server.execute({
		'class'   : "od",
		'levelist': "1/2/5/10/20/30/50/70/100/150/200/250/300/400/500/700/850/925/1000",
		'fcmonth' : "1/2/3/4/5/6/7",
		'number'  : "0/to/50",   
		'levtype' : "pl",
		'expver'  : "1",
		'method'  : "1",
		'origin'  : "ecmf",    
		'date'    : "%s"%(initbdate), 
		'time'    : "00:00:00",   
		'type'    : "fcmean",
		'param'   : "129.128", 
		'stream'  : "msmm",
		'system'  : "4",
		'grid'    : "128",},
		 "/storage/silver/acrcc/vg140344/HGT_S4Hindcasts_monthly_Mlevels_%s.grib"%(initbdate))
for year in range(2013, 2018):
	#Select year for data extraction
	print ('YEAR ',year)
	
	initmonth = 8
	initbdate = "%s%02d01"%(year,initmonth)
	print ("######### ERA-interim  #########")
	print ('get data from ', initbdate)
	print ("################################")		
	server.execute({
		'class'   : "od",
		'levelist': "1/2/5/10/20/30/50/70/100/150/200/250/300/400/500/700/850/925/1000",
		'fcmonth' : "1/2/3/4/5/6/7",
		'number'  : "0/to/50",   
		'levtype' : "pl",
		'expver'  : "1",
		'method'  : "1",
		'origin'  : "ecmf",    
		'date'    : "%s"%(initbdate), 
		'time'    : "00:00:00",   
		'type'    : "fcmean",
		'param'   : "131/132", 
		'stream'  : "msmm",
		'system'  : "4",
		'grid'    : "128",},
		 "/storage/silver/acrcc/vg140344/UV_S4Hindcasts_monthly_Mlevels_%s.grib"%(initbdate))
