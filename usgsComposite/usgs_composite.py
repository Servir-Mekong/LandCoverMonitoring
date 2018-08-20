# -*- coding: utf-8 -*-
"""
ExportLandsatSRComposite.py, SERVIR-Mekong (2017-07-30)

export landsat composites with gapfilling
________________________________________________________________________________


Usage
------

$ python  model.py {options}

{options} include:

--year (-y)	  : required
				 : year to create the Landsat composite image for
				 : in format YYYY

--season (-s)	: season to create the Landsat composite image for
				 : seasons are set for the Mekong region and will need to be \
				 : changed for other geographic areas
				 : options are 'drycool', 'dryhot', or 'rainy'


Example Usage
-------------

1) export surface reflectance composite for dryhot season of 2000 to assets:

  $ python model.py -y 2000 -s drycool 

"""

import ee
import logging
import time
import math
import argparse
from utils import *
import sun_angles
import view_angles

class environment(object):
	
	
	def __init__(self):
		"""Initialize the environment."""   
		 
		# Initialize the Earth Engine object, using the authentication credentials.
				
		ee.Initialize()
	   
		self.timeString = time.strftime("%Y%m%d_%H%M%S")

		# SEASONS:
		# '0': Dry Cool: Nov - Feb (305 - 59)
		# '1': Dry Hot: Mar - Apr (60 - 120)
		# '2': Rainy: May - Oct (121 - 304)
	
		startjulian = {'drycool':305,'dryhot':60,'rainy':121}
		endjulian = {'drycool':59,'dryhot':120,'rainy':304}
   
		# set dates
		self.startYear = int(args.year);
		self.endYear = int(args.year);
		self.startJulian = startjulian[args.season]
		self.endJulian = endjulian[args.season]
		
		if args.season == 'drycool':
			self.startYear = int(args.year)-1
		
		# filter for Nepal
		countries = ee.FeatureCollection('ft:1tdSwUL7MVpOauSgRzqVTOwdfy17KDbw-1d9omPw');
		country_names = ['Nepal'];	
		myCountry = countries.filter(ee.Filter.inList('Country', country_names));	
		aoi = myCountry.geometry().buffer(20000);
 
		self.dem = ee.Image("USGS/SRTMGL1_003")
 
		collectionName = "users/servirmekong/Nepal/" + args.season 
		
		#aoi = ee.Geometry.Polygon([[[84.9601342942816,28.412985417510814],[84.8887231614691,27.310613144962876],[85.4380395677191,27.12009683732014],[85.7841089036566,28.161451992773063],[84.9601342942816,28.412985417510814]]])
		self.region = aoi 
		#self.mekongRegion = ee.FeatureCollection('ft:1LEGeqwlBCAlN61ie5ol24NdUDqB1MgpFR_sJNWQJ')
		
		self.collection = ee.ImageCollection(collectionName)
		
		# variables for the tdom filter
		self.applyTDOM = True
		self.TDOMyears = 10
		self.shadowSumBands = ['nir','swir1'];
		self.zScoreThresh = -0.8
		self.shadowSumThresh = 0.35;
		self.dilatePixels = 2
		
		#users/servirmekong/usgs_sr_composites/drycool
		self.outputName = args.season + str(self.startYear) + "_" + str(self.endYear)

		# variable to filter cloud threshold
		self.metadataCloudCoverMax = 60
		
		# threshold for landsatCloudScore
		self.cloudThreshold = 10
		self.hazeThresh = 200
		
		# apply a filter to filter for high values
		self.filterPercentile = False
		self.filterPercentileYears = 25
		# percentiles to filter for bad data
		self.lowPercentile = 2
		self.highPercentile = 80

		# whether to use imagecolletions
		self.useL4=True
		self.useL5=True
		self.useL7=True
		self.useL7scanline = False
		self.useL8=True

		# On May 31, 2003 the Scan Line Corrector (SLC) in the ETM+ instrument failed
		self.l7Failed = ee.Date.fromYMD(2003,5,31)

		# apply cloud masks
		self.maskSR = True

		# get indicices 
		self.calcIndices = True

		# bands for tasselcap !maybe move
		self.tcInputBands =  ee.List(['blue','green','red','nir','swir1','swir2'])
		
		# bands to select
		self.bandNamesLandsat = ee.List(['blue','green','red','nir','swir1','thermal','swir2','sr_atmos_opacity','pixel_qa','radsat_qa'])
		
		# bands for export
		self.exportBands = ee.List(['blue','green','red','nir','swir1','thermal','swir2'])
		
		# bands for dividing
		self.divideBands = ee.List(['blue','green','red','nir','swir1','swir2'])
		
		#bands for stdev
		self.stdDevBands = ee.List(['blue','green','red','nir','swir1','thermal','swir2']) #,'ND_nir_red','ND_nir_swir2','ND_green_swir1']);
		self.stdDevExportsBands = ee.List(['blue_stdev','green_stdev','red_stdev','nir_stdev','swir1_stdev','thermal_stdev','swir2_stdev']) #,'ND_nir_red','ND_nir_swir2','ND_green_swir1']);
		
		# calculate stdev for indices
		self.stdIndiceDevBands = ee.List(["ND_nir_swir2","ND_green_swir1","ND_nir_red"])
		self.stdIndiceDevBandsExport = ee.List(["ND_nir_swir2_stdDev","ND_green_swir1_stdDev","ND_nir_red_stdDev"])

		# apply defringe
		self.defringe = True
		
		# pixel size
		self.pixSize = 30
		
		# user ID
		self.userID = "users/servirmekong/nepal/SR"  
		
		self.landsat4count = 0
		self.landsat5count = 0
		self.landsat7count = 0
		self.landsat8count = 0
	   
		# define the landsat bands								   
		self.sensorBandDictLandsatSR = ee.Dictionary({'L8' : ee.List([1,2,3,4,5,7,6,9,10,11]),
													  'L7' : ee.List([0,1,2,3,4,5,6,7,9,10]),
													  'L5' : ee.List([0,1,2,3,4,5,6,7,9,10]),
													  'L4' : ee.List([0,1,2,3,4,5,6,7,9,10])})

		# just placeholders for now
		self.calcMedoid = False
		self.calcMedian = True
		self.calcMean = False
		self.terrainCorrection = True
		self.brdfCorrection = True

		self.fillGaps = False
		self.fillGapYears = 2

		# threshold for defringing landsat5 and 7
		self.fringeCountThreshold = 279

		self.k = ee.Kernel.fixed(41, 41, 
								[[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
								[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]);
		
		
		
class SurfaceReflectance():
	
	def __init__(self):
		"""Initialize the Surfrace Reflectance app."""  
		
		# import the log library
		
		import logging
		
		# get the environment
		self.env = environment()
  
		# set the logfile
		logging.basicConfig(filename=str(self.env.timeString) + 'model.log', filemode='w', level=logging.DEBUG)
				
	
	def RunModel(self):
		"""Run the SR model"""  
		
		geo = self.env.region.getInfo()['coordinates']
		self.env.location = ee.Geometry.Polygon(geo) #ee.Geometry.Polygon(self.env.NgheAn)
		#self.env.location = ee.Geometry.Polygon([[104.260,22.553],[104.298,20.997],[106.556,20.997],[106.446,22.527],[104.260,22.553]])

		logging.info('starting the model the model')
				
		print self.env.startJulian, self.env.endJulian
   
		# construct date objects
		startDate = ee.Date.fromYMD(self.env.startYear,1,1).advance(self.env.startJulian,'day')
		endDate = ee.Date.fromYMD(self.env.endYear,1,1).advance(self.env.endJulian-1,'day')
		
		logging.info('startDate = ' + str(startDate.getInfo()))
		logging.info('endDatDate = ' + str(endDate.getInfo()))
		logging.info('Cloudcover filter = ' + str(self.env.metadataCloudCoverMax))		
		
		# get the images
		collection = self.GetLandsat(startDate,endDate,self.env.metadataCloudCoverMax,0,True).select(self.env.exportBands)
		

		count = ee.Image(collection.select(['blue']).count().rename('count')).int16()

		collection = collection.map(self.maskClouds)   
		
		count = collection.size();
		print('counted ' + str(count.getInfo()) +' images');	
		
		if self.env.applyTDOM:
			# years before and after
			y = int(self.env.TDOMyears)
			start = ee.Date.fromYMD(self.env.startYear-y,1,1)
			end = ee.Date.fromYMD(self.env.startYear+y,1,1)
			self.fullCollection = self.returnCollection(start,end).select(self.env.exportBands)
			collection = self.maskShadows(collection)
		
		# filter for high outliers
		if self.env.filterPercentile:
			# years before and after
			logging.info("high percentile: " + str(self.env.highPercentile))
			logging.info("low percentile: " + str(self.env.lowPercentile))
			y = int(self.env.filterPercentileYears / 2)
			start = ee.Date.fromYMD(self.env.startYear-y,1,1)
			end = ee.Date.fromYMD(self.env.startYear+y,1,1)
			self.fullCollection = self.returnCollection(start,end).select(self.env.exportBands)		 
			self.percentile = self.fullCollection.reduce(ee.Reducer.percentile([self.env.lowPercentile,self.env.highPercentile])) 
			collection = collection.map(self.MaskPercentile) 

		

		if self.env.calcIndices:
			indices = collection.map(self.addIndices)
			stdDevIndiceComposite = indices.select(self.env.stdIndiceDevBands).reduce(ee.Reducer.stdDev()).select(self.env.stdIndiceDevBandsExport)
		
		
		stdDevComposite = collection.select(self.env.stdDevBands).reduce(ee.Reducer.stdDev());

			 
		if self.env.calcMedoid:
			img = self.medoidMosaic(collection) 
			compositeStyle = "Medoid"
		
		if self.env.calcMedian:
			img = collection.median()
			compositeStyle = "Median"
		
		img = img.addBands(stdDevComposite)

		if self.env.calcIndices:
			print "add stdev indice to composite"
		   
			img = img.addBands(stdDevIndiceComposite)

		self.env.outputName = self.env.outputName + compositeStyle
		
		print "starting .. " + self.env.outputName

		if self.env.fillGaps: 
			gapfilter = ee.Image(self.env.startYear).updateMask(img.select("blue").mask())
			img = img.addBands(gapfilter.rename(['gapfill']))
			for i in range(1,self.env.fillGapYears,1):
				print i, self.env.fillGapYears
				img = self.unmaskYears(img,i)	
				img = self.unmaskFutureYears(img,i)	

		# rescale to save as int16
		img = ee.Image(self.reScaleLandsat(img))

		# export image
		self.ExportToAsset(img,self.env.outputName)
				 

	def GetLandsat(self,startDate,endDate,metadataCloudCoverMax,year,terrain=False):
		"""Get the Landsat imagery"""  
		
		logging.info('getting landsat images')
			
		# boolean to merge Landsat; when true is merges with another collection
		merge = False
		
		if self.env.useL4:   
			landsat4 =  ee.ImageCollection('LANDSAT/LT04/C01/T1_SR').filterDate(startDate,endDate).filterBounds(self.env.location)
			landsat4 = landsat4.filter(ee.Filter.calendarRange(self.env.startJulian,self.env.endJulian))	 
			landsat4 = landsat4.filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)
		
		self.env.landsat4count += int(landsat4.size().getInfo())
		if landsat4.size().getInfo() > 0:
			if self.env.defringe == True:
				 landsat4 =  landsat4.map(self.DefringeLandsat)
			if self.env.maskSR == True:
				landsat4 = landsat4.map(self.CloudMaskSR)
				landsat4 = landsat4.map(self.maskHaze)  
			if not merge:
				landsatCollection = landsat4.select(self.env.sensorBandDictLandsatSR.get('L4'),self.env.bandNamesLandsat)
				merge = True
		
		# landsat 5 image collections 
		if self.env.useL5:
			landsat5 =  ee.ImageCollection('LANDSAT/LT05/C01/T1_SR').filterDate(startDate,endDate).filterBounds(self.env.location)
			landsat5 = landsat5.filter(ee.Filter.calendarRange(self.env.startJulian,self.env.endJulian))
			landsat5 = landsat5.filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)
			self.env.landsat5count += int(landsat5.size().getInfo())
			if landsat5.size().getInfo() > 0:
				if self.env.defringe == True:
					 landsat5 =  landsat5.map(self.DefringeLandsat)
				if self.env.maskSR == True:
					landsat5 = landsat5.map(self.CloudMaskSR)
					landsat5 = landsat5.map(self.maskHaze)  
				if not merge:
					landsatCollection = landsat5.select(self.env.sensorBandDictLandsatSR.get('L5'),self.env.bandNamesLandsat)
					merge = True
				else:
					landsatCollection = landsatCollection.merge(landsat5.select(self.env.sensorBandDictLandsatSR.get('L5'),self.env.bandNamesLandsat))
		
		# landsat 7 image collections  
		l7slm = True
		if self.env.useL7:
			landsat7 =  ee.ImageCollection('LANDSAT/LE07/C01/T1_SR').filterDate(startDate,endDate).filterBounds(self.env.location)
			landsat7 = landsat7.filter(ee.Filter.calendarRange(self.env.startJulian,self.env.endJulian))
			if self.env.startYear == 2003 or self.env.endYear == 2003:
				if self.env.startJulian > 151:
						l7slm = False 
				if self.env.useL7scanline == False:
					landsat7 = landsat7.filterDate(startDate,self.env.l7Failed)
			if self.env.startYear+year > 2003 and self.env.useL7scanline == False:
				l7slm = False
			if l7slm == True:
				landsat7 = landsat7.filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)
				self.env.landsat7count += int(landsat7.size().getInfo())
				if landsat7.size().getInfo() > 0:
					if self.env.defringe == True:
							landsat7 =  landsat7.map(self.DefringeLandsat)			
					if self.env.maskSR == True:
						landsat7 = landsat7.map(self.CloudMaskSR)
						landsat7 = landsat7.map(self.maskHaze)  
					if not merge:
						landsatCollection = landsat7.select(self.env.sensorBandDictLandsatSR.get('L7'),self.env.bandNamesLandsat)
						merge = True
					else:
						landsatCollection = landsatCollection.merge(landsat7.select(self.env.sensorBandDictLandsatSR.get('L7'),self.env.bandNamesLandsat))

		# landsat8  image collections										 
		if self.env.useL8:
			landsat8 =  ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filterDate(startDate,endDate).filterBounds(self.env.location)
			landsat8 = landsat8.filter(ee.Filter.calendarRange(self.env.startJulian,self.env.endJulian))
			landsat8 = landsat8.filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)
			self.env.landsat8count += int(landsat8.size().getInfo())
			if landsat8.size().getInfo() > 0:
				if self.env.maskSR == True:
					landsat8 = landsat8.map(self.CloudMaskSRL8)	
				if not merge:
					landsatCollection = landsat8.select(self.env.sensorBandDictLandsatSR.get('L8'),self.env.bandNamesLandsat)
					merge = True
				else:
					landsatCollection = landsatCollection.merge(landsat8.select(self.env.sensorBandDictLandsatSR.get('L8'),self.env.bandNamesLandsat))			


		landsatCollection = landsatCollection.map(self.ScaleLandsat)	  
		
		
		
		if terrain:
			if self.env.brdfCorrection == True:
				print "applying brdf correction"
				landsatCollection = landsatCollection.map(self.brdf)
		
			if self.env.terrainCorrection == True:
				print "applying terrain correction"
				landsatCollection = ee.ImageCollection(landsatCollection.map(self.terrain))

		if merge:
			count = landsatCollection.size(); 
			return ee.ImageCollection(landsatCollection)
		
		else:
			return ee.ImageCollection([ee.Image(0)])
			
		
	   

	def returnCollection(self,start,end):
		"""Calculate percentiles to filter imagery"""  
		
		logging.info('calculate percentiles')
		
		#startDate = ee.Date.fromYMD(1984,1,1)
		#endDate = ee.Date.fromYMD(2099,1,1)	
		cloudCoverMax = self.env.metadataCloudCoverMax 
		
		# get the images
		collection = self.GetLandsat(start,end,cloudCoverMax,0)
		collection = collection.map(self.maskClouds)   
		
		return collection
	   
	def CloudMaskSR(self,img):
		 """apply cf-mask Landsat""" 
		 
		 QA = img.select("pixel_qa")
		 #mask = ee.Image(self.getQABits(QA,3, 5, 'internal_quality_flag')); 
		
		 return img.updateMask(QA.lt(112)).copyProperties(img)
		 #return img.addBands(mask.select('internal_quality_flag'))

	def CloudMaskSRL8(self,img):
		 """apply cf-mask Landsat""" 
		 
		 QA = img.select("pixel_qa")
		 #mask = ee.Image(self.getQABits(QA,3, 5, 'internal_quality_flag')); 
		 # clear and water
		 mask = QA.lt(325)
		 snow = QA.eq(480)
		 
		 mask = mask.add(snow)
		 
		 return img.updateMask(mask).copyProperties(img)

	def ScaleLandsat(self,img):
		"""Landast is scaled by factor 0.0001 """
		
		thermal = ee.Image(img).select(ee.List(['thermal'])).multiply(0.1)
		scaled = ee.Image(img).select(self.env.divideBands).multiply(0.0001)
		image = ee.Image(scaled.addBands(thermal))
		date = img.get("system:time_start")
	
		logging.info('return scaled image')
		return ee.Image(image.copyProperties(img)).set("system:time_start",date).set('system:footprint',img.get('system:footprint'))

	def reScaleLandsat(self,img):
		"""Landast is scaled by factor 0.0001 """
		
		thermalBand = ee.List(['thermal','thermal_stdDev'])
		#gapfillBand = ee.List(['gapfill'])
		thermal = ee.Image(img).select(thermalBand).multiply(10)
		#gapfill = ee.Image(img).select(gapfillBand)
				
		
		otherBands = ee.Image(img).bandNames().removeAll(thermalBand) #.removeAll(gapfillBand)
		scaled = ee.Image(img).select(otherBands).divide(0.0001)
		
		image = ee.Image(scaled.addBands(thermal)).int16()
		logging.info('return scaled image')
		
		return image.copyProperties(img)

	def DefringeLandsat(self,img):   
		"""remove scanlines from landsat 4 and 7 """
		
		logging.info('removing scanlines')

		m = ee.Image(img).mask().reduce(ee.Reducer.min())
		sum = m.reduceNeighborhood(ee.Reducer.sum(), self.env.k, 'kernel')
		mask = sum.gte(self.env.fringeCountThreshold)		
		#img = img.mask(img.mask().add(sum)) 
		return img.updateMask(mask)
		
	def MaskPercentile(self,img):
		"""mask light and dark pixels """
		logging.info('mask light and dark areas')
		
		# GEE adds _p_percentile to bandname
				
		upper = str(int(self.env.highPercentile))
		lower = str(int(self.env.lowPercentile))
		
		# we dont want to include the mask
		selectedBandNamesLandsat =['blue','green','red','nir','swir1','swir2']

		for b in selectedBandNamesLandsat:
			
			selectedBand = ee.List([b])

			# get the upper and lower band
			bandsUpper = ee.List([str(b)+ '_p'+ upper])#+ upper,'green_p'+ upper,'red_p'+ upper,'nir_p'+ upper,'swir1_p'+ upper,'swir2_p'+ upper])
			bandsLower = ee.List([str(b)+ '_p' + lower])#+ lower,'green_p'+ lower,'red_p'+ lower,'nir_p'+ lower,'swir1_p'+ lower,'swir2_p'+ lower])

			percentilesUp = self.percentile.select(bandsUpper,selectedBand  )
			percentilesLow = self.percentile.select(bandsLower,selectedBand )
			
			#print percentilesUp
			imgToMask = img.select(selectedBand)
							 
			darkMask = ee.Image(imgToMask.lt(percentilesLow).reduce(ee.Reducer.sum())).eq(0)
			lightMask = ee.Image(imgToMask.gt(percentilesUp).reduce(ee.Reducer.sum())).eq(0)
			
			#img = ee.Image(img.updateMask(lightMask).updateMask(darkMask).copyProperties(img))
			img = ee.Image(img.updateMask(lightMask).copyProperties(img))
			
		return img.copyProperties(img)

	def maskHaze(self,img):
		""" mask haze """
		opa = ee.Image(img.select(['sr_atmos_opacity']).multiply(0.001))
		haze = opa.gt(self.env.hazeThresh)
		return img.updateMask(haze.Not())
 

	def maskClouds(self,img):

		score = ee.Image(1.0);
		# Clouds are reasonably bright in the blue band.
		blue_rescale = img.select('blue').subtract(ee.Number(0.1)).divide(ee.Number(0.3).subtract(ee.Number(0.1)))
		score = score.min(blue_rescale);

		# Clouds are reasonably bright in all visible bands.
		visible = img.select('red').add(img.select('green')).add(img.select('blue'))
		visible_rescale = visible.subtract(ee.Number(0.2)).divide(ee.Number(0.8).subtract(ee.Number(0.2)))
		score = score.min(visible_rescale);

		# Clouds are reasonably bright in all infrared bands.
		infrared = img.select('nir').add(img.select('swir1')).add(img.select('swir2'))
		infrared_rescale = infrared.subtract(ee.Number(0.3)).divide(ee.Number(0.8).subtract(ee.Number(0.3)))
		score = score.min(infrared_rescale);

		# Clouds are reasonably cool in temperature.
		temp_rescale = img.select('thermal').subtract(ee.Number(300)).divide(ee.Number(290).subtract(ee.Number(300)))
		score = score.min(temp_rescale);

		# However, clouds are not snow.
		ndsi = img.normalizedDifference(['green', 'swir1']);
		ndsi_rescale = ndsi.subtract(ee.Number(0.8)).divide(ee.Number(0.6).subtract(ee.Number(0.8)))
		score =  score.min(ndsi_rescale).multiply(100).byte();
		mask = score.lt(self.env.cloudThreshold).rename(['cloudMask']);
		img = img.updateMask(mask);
		
		return img;

	def maskShadows(self,collection,zScoreThresh=-0.8,shadowSumThresh=0.35,dilatePixels=2):

		def TDOM(image):
			zScore = image.select(shadowSumBands).subtract(irMean).divide(irStdDev)
			irSum = image.select(shadowSumBands).reduce(ee.Reducer.sum())
			TDOMMask = zScore.lt(zScoreThresh).reduce(ee.Reducer.sum()).eq(2)\
				.And(irSum.lt(shadowSumThresh)).Not()
			TDOMMask = TDOMMask.focal_min(dilatePixels)
			
			return image.updateMask(TDOMMask)

		shadowSumBands = ['nir','swir1']

		# Get some pixel-wise stats for the time series
		irStdDev = self.fullCollection.select(shadowSumBands).reduce(ee.Reducer.stdDev())
		irMean = self.fullCollection.select(shadowSumBands).reduce(ee.Reducer.mean())

		# Mask out dark dark outliers
		collection_tdom = collection.map(TDOM)

		return collection_tdom

	def ExportToAsset(self,img,assetName):  
		"""export to asset """
		
		outputName = self.env.userID + assetName 
		logging.info('export image to asset: ' + str(outputName))   
		
		startDate = ee.Date.fromYMD(self.env.startYear,1,1)
		endDate = ee.Date.fromYMD(self.env.endYear,12,31)	

		compositingMethod = ""
		if self.env.calcMedoid == True:
			compositingMethod = "Medoid"
		if self.env.calcMedian == True:
			compositingMethod = "Median"			

		image = ee.Image(img).set({'system:time_start':startDate.millis(), \
									'startyear':self.env.startYear, \
									'endyear':self.env.endYear, \
									'startJulian':self.env.startJulian, \
									'endJulian':self.env.endJulian,
									'source':'USGS_SR',\
									'median':self.env.calcMedian,\
									'mean':self.env.calcMean,\
									'medoid':self.env.calcMedoid,\
									'count_landsat_4':str(self.env.landsat4count),\
									'count_landsat_5':str(self.env.landsat5count),\
									'count_landsat_7':str(self.env.landsat7count),\
									'count_landsat_8':str(self.env.landsat8count),\
									'cloud_threshold':self.env.metadataCloudCoverMax ,\
									'defringe':str(self.env.defringe),\
									'apply_usgs_cloud_filter': str(self.env.maskSR ),\
									'applyTDOM':self.env.applyTDOM,\
									'TDOMyears':str(self.env.TDOMyears),\
									'pixel_size':self.env.pixSize,\
									'zScoreThresh':self.env.zScoreThresh,\
									'shadowSumThresh':self.env.shadowSumThresh,\
									'dilatePixels':self.env.dilatePixels,\
									'filter_percentile':str(self.env.filterPercentile),\
									'filter_percentile_years':self.env.filterPercentileYears,\
									'upper_percentile': self.env.highPercentile,\
									'calculate_indices':str(self.env.calcIndices),\
									'gap_filling':str(self.env.fillGaps),\
									'landsat_7_scanline':str(self.env.useL7scanline),\
									'years_of_gap_filling':self.env.fillGapYears,\
									'compositingMethod':compositingMethod,\
									'version':'1.0'}).clip(self.env.region) 


		#task_ordered = ee.batch.Export.image.toAsset(image=ee.Image(img), description=str(self.env.timeString)+assetName, assetId=outputName,region=self.env.location['coordinates'], maxPixels=1e13,scale=self.env.pixSize)
		import time
		t= time.strftime("%Y%m%d_%H%M%S")

		task_ordered = ee.batch.Export.image.toAsset(image=ee.Image(image), description=assetName+t, assetId=outputName+t,region=self.env.location['coordinates'], maxPixels=1e13,scale=self.env.pixSize)
		
		# start task
		task_ordered.start() 
	
	def addIndices(self,img):
		""" Function to add common (and less common) spectral indices to an image.
			Includes the Normalized Difference Spectral Vector from (Angiuli and Trianni, 2014) """
			
		img = img.addBands(img.normalizedDifference(['green','swir1']).rename(['ND_green_swir1']));  # NDSI, MNDWI
		img = img.addBands(img.normalizedDifference(['nir','red']).rename(['ND_nir_red']));  # NDVI
		img = img.addBands(img.normalizedDifference(['nir','swir2']).rename(['ND_nir_swir2']));  # NBR, MNDVI

		return img;	   


	def unmaskYears(self,img,year):
		""" Function to unmask nodata withpixels previous year """
		
		
		if self.env.startYear-year > 1984:
			print "unmasking for year " + str(self.env.startYear-year) 
			startDate = ee.Date.fromYMD(self.env.startYear-year,1,1)
			endDate = ee.Date.fromYMD(self.env.endYear-year,12,31)	
			prev = self.GetLandsat(startDate,endDate,self.env.metadataCloudCoverMax,year)
			if int(prev.size().getInfo()) > 1:

				prev = self.maskShadows(prev.select(self.env.exportBands))
				
				if self.env.filterPercentile:
					prev = prev.map(self.MaskPercentile)
				
				stdDevComposite = prev.select(self.env.stdDevBands).reduce(ee.Reducer.stdDev());	 

				if self.env.calcMedoid:
					previmg = self.medoidMosaic(prev) 
				
				if self.env.calcMedian:
					previmg = prev.median()		 
				
				previmg = previmg.mask(previmg.select("blue").gt(0))
				print "added start minux year ", self.env.startYear-year
				gapfilter = ee.Image(self.env.startYear-year).updateMask(previmg.select("blue").mask())
				previmg = previmg.addBands(gapfilter.rename(['gapfill'])).addBands(stdDevComposite)
				
				if self.env.calcIndices:
					indices = prev.map(self.addIndices)
					stdDevIndiceComposite = indices.select(self.env.stdIndiceDevBands).reduce(ee.Reducer.stdDev()).select(self.env.stdIndiceDevBandsExport)
					previmg = previmg.addBands(stdDevIndiceComposite)
				
				img = img.unmask(previmg)
			
		return ee.Image(img)

	def unmaskFutureYears(self,img,year):
		""" Function to unmask nodata withpixels future year """
		
		
		if self.env.startYear+year < 2018:
			print "unmasking for year " + str(self.env.startYear+year) 
			startDate = ee.Date.fromYMD(self.env.startYear+year,1,1)
			endDate = ee.Date.fromYMD(self.env.endYear+year,12,31)	
			prev = self.GetLandsat(startDate,endDate,self.env.metadataCloudCoverMax,year)
			if int(prev.size().getInfo()) > 1:
				prev = self.maskShadows(prev.select(self.env.exportBands))
				stdDevComposite = prev.select(self.env.stdDevBands).reduce(ee.Reducer.stdDev()); 

				if self.env.filterPercentile:
					prev = prev.map(self.MaskPercentile) 

				if self.env.calcMedoid:
					previmg = self.medoidMosaic(prev)			 
				
				if self.env.calcMedian:
					previmg = prev.median()			 
				
				
				previmg = previmg.mask(previmg.select("blue").gt(0))
				print "added start year ", self.env.startYear+year
				gapfilter = ee.Image(self.env.startYear+year).updateMask(previmg.select("blue").mask())

				if self.env.calcIndices:
					indices = prev.map(self.addIndices)
					stdDevIndiceComposite = indices.select(self.env.stdIndiceDevBands).reduce(ee.Reducer.stdDev()).select(self.env.stdIndiceDevBandsExport)
					previmg = previmg.addBands(stdDevIndiceComposite)

				
				previmg = previmg.addBands(gapfilter.rename(['gapfill'])).addBands(stdDevComposite)
				
				if self.env.calcIndices:
					previmg.addBands(stdDevIndiceComposite)
				
				img = ee.Image(img).unmask(previmg)
		
		return ee.Image(img)

	def makeTiles(self):
		""" Function to tile up """
		# set bounds
		xmin = 91.979254;
		xmax = 114.664186;
		ymin = 5.429208;
		ymax = 28.728774;

		# number of rows and columns
		n = 1;

		for  i in range(0, n, 1):
			for j in range(0, n,1):		# x, y distance of one block
				xs = (xmax - xmin) / n
				ys = (ymax - ymin) / n
				
				xl = xmin + i * xs;
				xr = xmin + (i+1) * xs;
				
				yt = ymin + (j*ys);
				yb = ymin + (j+1)*ys;
				geom =   [[xl, yt], [xl, yb], [xr, yb], [xr, yt]];
				print geom

				col = SurfaceReflectance().RunModel(geom,i,j)


	def terrain(self,img):   
		
		
		degree2radian = 0.01745;
		thermalBand = img.select(['thermal'])
		img = img.unmask(0)
		dem = self.env.dem.unmask(0)	
		geom = ee.Geometry(img.get('system:footprint')).bounds().buffer(1000) 
		
		def topoCorr_IC(img):
				
			

			# Extract image metadata about solar position
			SZ_rad = ee.Image.constant(ee.Number(img.get('SOLAR_ZENITH_ANGLE'))).multiply(degree2radian).clip(geom); 
			SA_rad = ee.Image.constant(ee.Number(img.get('SOLAR_AZIMUTH_ANGLE'))).multiply(degree2radian).clip(geom); 
				
					
			# Creat terrain layers
			slp = ee.Terrain.slope(dem).clip(img.geometry().buffer(10000));
			slp_rad = ee.Terrain.slope(dem).multiply(degree2radian).clip(geom);
			asp_rad = ee.Terrain.aspect(dem).multiply(degree2radian).clip(geom);
	  
	  
				
			# Calculate the Illumination Condition (IC)
			# slope part of the illumination condition
			cosZ = SZ_rad.cos();
			cosS = slp_rad.cos();
			slope_illumination = cosS.expression("cosZ * cosS", \
							{'cosZ': cosZ, 'cosS': cosS.select('slope')});
				
				
			# aspect part of the illumination condition
			sinZ = SZ_rad.sin(); 
			sinS = slp_rad.sin();
			cosAziDiff = (SA_rad.subtract(asp_rad)).cos();
			aspect_illumination = sinZ.expression("sinZ * sinS * cosAziDiff", \
						{'sinZ': sinZ, \
						'sinS': sinS, \
						'cosAziDiff': cosAziDiff});
				
			# full illumination condition (IC)
			ic = slope_illumination.add(aspect_illumination);
				
				

			# Add IC to original image
			img_plus_ic = ee.Image(img.addBands(ic.rename(['IC'])).addBands(cosZ.rename(['cosZ'])).addBands(cosS.rename(['cosS'])).addBands(slp.rename(['slope'])));
				
			return ee.Image(img_plus_ic);
	 
		def topoCorr_SCSc(img):
			img_plus_ic = img;
			mask1 = img_plus_ic.select('nir').gt(-0.1);
			mask2 = img_plus_ic.select('slope').gte(5) \
					.And(img_plus_ic.select('IC').gte(0)) \
					.And(img_plus_ic.select('nir').gt(-0.1));

			img_plus_ic_mask2 = ee.Image(img_plus_ic.updateMask(mask2));

			bandList = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2']; # Specify Bands to topographically correct
		

			def applyBands(image):
				blue = apply_SCSccorr('blue').select(['blue'])
				green = apply_SCSccorr('green').select(['green'])
				red = apply_SCSccorr('red').select(['red'])
				nir = apply_SCSccorr('nir').select(['nir'])
				swir1 = apply_SCSccorr('swir1').select(['swir1'])
				swir2 = apply_SCSccorr('swir2').select(['swir2'])
				return replace_bands(image, [blue, green, red, nir, swir1, swir2])

			def apply_SCSccorr(band):
				method = 'SCSc';
				
				out = img_plus_ic_mask2.select('IC', band).reduceRegion(reducer= ee.Reducer.linearFit(), \
				geometry= ee.Geometry(geom.buffer(-5000)), \
							  scale= 100, \
							  maxPixels = 100000000); 
										  	
				out_a = ee.Number(out.get('scale'));
				out_b = ee.Number(out.get('offset'));
				out_c = ee.Number(out.get('offset')).divide(ee.Number(out.get('scale')));
						
				# apply the SCSc correction
				SCSc_output = img_plus_ic_mask2.expression("((image * (cosB * cosZ + cvalue)) / (ic + cvalue))", {
										'image': img_plus_ic_mask2.select([band]),
										'ic': img_plus_ic_mask2.select('IC'),
										'cosB': img_plus_ic_mask2.select('cosS'),
										'cosZ': img_plus_ic_mask2.select('cosZ'),
										'cvalue': out_c });
		  
				return ee.Image(SCSc_output);
																	  
																	  
		
			img_SCSccorr = applyBands(img).select(bandList).addBands(img_plus_ic.select('IC'))
			
			bandList_IC = ee.List([bandList, 'IC']).flatten();
				
			img_SCSccorr = img_SCSccorr.unmask(img_plus_ic.select(bandList_IC)).select(bandList);
				
			return img_SCSccorr.unmask(img_plus_ic.select(bandList)) 
		
		img = topoCorr_IC(img)
		img = topoCorr_SCSc(img)
			
		return img.addBands(thermalBand)		
		
	def brdf(self,img):   
		
		def _apply(image, kvol, kvol0):
			blue = _correct_band(image, 'blue', kvol, kvol0, f_iso=0.0774, f_geo=0.0079, f_vol=0.0372)
			green = _correct_band(image, 'green', kvol, kvol0, f_iso=0.1306, f_geo=0.0178, f_vol=0.0580)
			red = _correct_band(image, 'red', kvol, kvol0, f_iso=0.1690, f_geo=0.0227, f_vol=0.0574)
			nir = _correct_band(image, 'nir', kvol, kvol0, f_iso=0.3093, f_geo=0.0330, f_vol=0.1535)
			swir1 = _correct_band(image, 'swir1', kvol, kvol0, f_iso=0.3430, f_geo=0.0453, f_vol=0.1154)
			swir2 = _correct_band(image, 'swir2', kvol, kvol0, f_iso=0.2658, f_geo=0.0387, f_vol=0.0639)
			return replace_bands(image, [blue, green, red, nir, swir1, swir2])


		def _correct_band(image, band_name, kvol, kvol0, f_iso, f_geo, f_vol):
			"""fiso + fvol * kvol + fgeo * kgeo"""
			iso = ee.Image(f_iso)
			geo = ee.Image(f_geo)
			vol = ee.Image(f_vol)
			pred = vol.multiply(kvol).add(geo.multiply(kvol)).add(iso).rename(['pred'])
			pred0 = vol.multiply(kvol0).add(geo.multiply(kvol0)).add(iso).rename(['pred0'])
			cfac = pred0.divide(pred).rename(['cfac'])
			corr = image.select(band_name).multiply(cfac).rename([band_name])
			return corr


		def _kvol(sunAz, sunZen, viewAz, viewZen):
			"""Calculate kvol kernel.
			From Lucht et al. 2000
			Phase angle = cos(solar zenith) cos(view zenith) + sin(solar zenith) sin(view zenith) cos(relative azimuth)"""
			
			relative_azimuth = sunAz.subtract(viewAz).rename(['relAz'])
			pa1 = viewZen.cos() \
				.multiply(sunZen.cos())
			pa2 = viewZen.sin() \
				.multiply(sunZen.sin()) \
				.multiply(relative_azimuth.cos())
			phase_angle1 = pa1.add(pa2)
			phase_angle = phase_angle1.acos()
			p1 = ee.Image(PI().divide(2)).subtract(phase_angle)
			p2 = p1.multiply(phase_angle1)
			p3 = p2.add(phase_angle.sin())
			p4 = sunZen.cos().add(viewZen.cos())
			p5 = ee.Image(PI().divide(4))

			kvol = p3.divide(p4).subtract(p5).rename(['kvol'])

			viewZen0 = ee.Image(0)
			pa10 = viewZen0.cos() \
				.multiply(sunZen.cos())
			pa20 = viewZen0.sin() \
				.multiply(sunZen.sin()) \
				.multiply(relative_azimuth.cos())
			phase_angle10 = pa10.add(pa20)
			phase_angle0 = phase_angle10.acos()
			p10 = ee.Image(PI().divide(2)).subtract(phase_angle0)
			p20 = p10.multiply(phase_angle10)
			p30 = p20.add(phase_angle0.sin())
			p40 = sunZen.cos().add(viewZen0.cos())
			p50 = ee.Image(PI().divide(4))

			kvol0 = p30.divide(p40).subtract(p50).rename(['kvol0'])

			return (kvol, kvol0)
         
 		date = img.date()

		footprint = determine_footprint(img)

		(sunAz, sunZen) = sun_angles.create(date, footprint)
		(viewAz, viewZen) = view_angles.create(footprint)
		(kvol, kvol0) = _kvol(sunAz, sunZen, viewAz, viewZen)
		
		image = _apply(img, kvol.multiply(PI()), kvol0.multiply(PI()))

		return image
		
	def medoidMosaic(self,collection):
		""" medoid composite with equal weight among indices """
		
		# calculate the median of temp band
		thermal = ee.ImageCollection(collection.select(['thermal'])).median()
				
		collection = collection.select(self.env.divideBands)

		bandNames = self.env.divideBands;
		bandNumbers = ee.List.sequence(1,bandNames.length());
		
		# calculate medion
		median = ee.ImageCollection(collection).median()
		
		def subtractmedian(img):
			diff = ee.Image(img).subtract(median).pow(ee.Image.constant(2));
			return diff.reduce('sum').addBands(img);
		
		medoid = collection.map(subtractmedian)
  
		medoid = ee.ImageCollection(medoid).reduce(ee.Reducer.min(bandNames.length().add(1))).select(bandNumbers,bandNames);
  
		return medoid.addBands(thermal);
	 
   
if __name__ == "__main__":
  
	# set argument parsing object
	parser = argparse.ArgumentParser(description="Create Landsat image composites using Google\
												  Earth Engine.")
   
	parser.add_argument('--year','-y', type=str,required=True,
						help="Year to perform the ats correction and save to asset format in 'YYYY'")

	parser.add_argument('--season','-s', choices=['drycool','dryhot','rainy'],type=str,
						help="Season to create composite for, these align with SERVIR-Mekong's seasonal composite times")

	#parser.add_argument('--user','-u', type=str, default="servir-mekong",choices=['servir-mekong','servirmekong',"ate","biplov","quyen","atesig","adpc","KSA"],
	#					help="specify user account to run task")

	args = parser.parse_args() # get arguments  
  
	# user account to run task on
	
	year = args.year
	seasons = args.season
	
	ee.Initialize()	
	SurfaceReflectance().RunModel()
