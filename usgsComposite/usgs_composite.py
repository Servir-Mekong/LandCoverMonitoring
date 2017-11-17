# -*- coding: utf-8 -*-
"""
ExportLandsatSRComposite.py, SERVIR-Mekong (2017-07-30)

export landsat composites with gapfilling
________________________________________________________________________________


Usage
------

$ python  model.py {options}

{options} include:

--year (-y)      : required
                 : year to create the Landsat composite image for
                 : in format YYYY

--season (-s)    : season to create the Landsat composite image for
                 : seasons are set for the Mekong region and will need to be \
                 : changed for other geographic areas
                 : options are 'drycool', 'dryhot', or 'rainy'

--user (-u)      : user account used to create the composite
                 : changes the ~/.config/earthengine/credentials file
                 : dictionary is called to get credentials
                 : options are servirmekong, servir-mekong, ate, biplov .. default is servir-mekong

Example Usage
-------------

1) export surface reflectance composite for dryhot season of 2000 to assets:

  $ python model.py -y 2000 -s drycool -u Quyen

"""

import ee
import logging
import time
import math
from usercredentials import addUserCredentials
import argparse

class environment(object):
    
    
    def __init__(self):
        """Initialize the environment."""   
         
        # Initialize the Earth Engine object, using the authentication credentials.
		
        ee.Initialize()
       
	self.timeString = time.strftime("%Y%m%d_%H%M%S")

	# SEASONS:
	# '0': Dry Cool: Nov - Feb (305 - 59)
	# '1': Dry Hot: Mar - Apr (60 - 181)
	# '2': Rainy: May - Oct (182 - 304)
   
	startjulian = {'drycool':305,'dryhot':60,'rainy':182}
	endjulian = {'drycool':59,'dryhot':181,'rainy':304}
   
        # set dates
        self.startYear = int(args.year);
        self.endYear = int(args.year);
        self.startJulian = startjulian[args.season]
        self.endJulian = endjulian[args.season]
	
	if args.season == 'drycool':
	    self.startYear = int(args.year)-1
        
	self.NgheAn = [[103.876,18.552],[105.806,18.552],[105.806,19.999],[103.876,19.999],[103.876,18.552]]
	
        collectionName = "projects/servir-mekong/usgs_sr_composites/" + args.season 
	
	self.mekongRegion = ee.FeatureCollection('ft:1LEGeqwlBCAlN61ie5ol24NdUDqB1MgpFR_sJNWQJ')

	
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
        self.metadataCloudCoverMax = 70
        
        # threshold for landsatCloudScore
        self.cloudThreshold = 1
        
	# apply a filter to filter for high values
	self.filterPercentile = True
	self.filterPercentileYears = 10
        # percentiles to filter for bad data
        self.lowPercentile = 1
        self.highPercentile = 66

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
	self.getIndices = False
	self.calcIndices = False
	self.calcTasselcap = False
	self.calcTCangles = False

	# bands for tasselcap !maybe move
	self.tcInputBands =  ee.List(['blue','green','red','nir','swir1','swir2'])
        
	# bands to select
	self.bandNamesLandsat = ee.List(['blue','green','red','nir','swir1','thermal','swir2','sr_atmos_opacity','pixel_qa','radsat_qa'])
        
	# bands for export
	self.exportBands = ee.List(['blue','green','red','nir','swir1','thermal','swir2'])
	
	# bands for dividing
	self.divideBands = ee.List(['blue','green','red','nir','swir1','swir2'])
        
	# apply defringe
        self.defringe = True
        
        # pixel size
        self.pixSize = 30
        
        # user ID
        #self.userID = "users/servirmekong/assemblage/"
        self.userID = "users/servirmekong/temp/1NAmedoid_SurfacReflectance"
        #self.userID = "projects/servir-mekong/usgs_sr_composites/" + args.season + "/" 
        self.userID = "projects/servir-mekong/temp/" + args.season 

       
        # define the landsat bands			           
        self.sensorBandDictLandsatSR = ee.Dictionary({'L8' : ee.List([1,2,3,4,5,7,6,9,10,11]),
                                                      'L7' : ee.List([0,1,2,3,4,5,6,7,9,10]),
                                                      'L5' : ee.List([0,1,2,3,4,5,6,7,9,10]),
                                                      'L4' : ee.List([0,1,2,3,4,5,6,7,9,10])})

	# just placeholders for now
	self.calcMedoid = True
	self.calcMedian = False
	self.calcMean = False

	self.fillGaps = True
	self.fillGapYears = 10

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
                
    
    def RunModel(self,geo,x,y):
        """Run the SR model"""  
        
	self.env.location = ee.Geometry.Polygon(geo) #ee.Geometry.Polygon(self.env.NgheAn)
	#self.env.location = ee.Geometry.Polygon(self.env.NgheAn)

	#self.env.location = ee.Geometry.Polygon([[102.716,18.615],[105.622,18.620],[105.540,22.451],[102.650,22.466],[104.716,18.615]])
	#self.env.outputName = self.env.outputName + str(x) + str(y)
	 
        logging.info('starting the model the model')
        	
	print self.env.startJulian, self.env.endJulian
   
	# construct date objects
        startDate = ee.Date.fromYMD(self.env.startYear,1,1)
        endDate = ee.Date.fromYMD(self.env.endYear,12,31)    
        
	logging.info('startDate = ' + str(startDate.getInfo()))
	logging.info('endDatDate = ' + str(startDate.getInfo()))
	logging.info('Cloudcover filter = ' + str(self.env.metadataCloudCoverMax))	
	
        # get the images
        collection = self.GetLandsat(startDate,endDate,self.env.metadataCloudCoverMax).select(self.env.exportBands)
	collection = collection.map(self.maskClouds)   
	
	count = collection.size();
        print('counted ' + str(count.getInfo()) +' images');    
	
  		
	if self.env.applyTDOM:
	    # years before and after
	    y = int(self.env.TDOMyears / 2)
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
	
	
	
	if self.env.calcMedoid:
	    img = self.medoidMosaic(collection) 
	        
	if self.env.fillGaps: 
	    gapfilter = img.select(["blue"]).gt(0).multiply(self.env.startYear)
	    img = img.addBands(gapfilter.rename(['gapfill']))
	    for i in range(1,self.env.fillGapYears,1):
		img = self.unmaskYears(img,i)    
	        img = self.unmaskFutureYears(img,i)    

	if self.env.getIndices:
	    img = self.getAllIndices(img)
	
	# rescale to save as int16
	img = ee.Image(self.reScaleLandsat(img))
	
	# export image
	self.ExportToAsset(img,self.env.outputName)
	

	
    def getAllIndices(self,img):
	""" get all indices"""

	if self.env.calcIndices:
	    img = self.addIndices(img)
	
	if self.env.calcTasselcap:
	    img = self.getTasseledCap(img,self.env.tcInputBands )
	
	if self.env.calcTCangles:
	    img = self.addTCAngles(img)
	
	return img
             

    def GetLandsat(self,startDate,endDate,metadataCloudCoverMax):
        """Get the Landsat imagery"""  
        
        logging.info('getting landsat images')
            
        # boolean to merge Landsat; when true is merges with another collection
        merge = False
	
	self.landsat4count = 0
	self.landsat5count = 0
	self.landsat7count = 0
	self.landsat8count = 0
	
	#print startDate, endDate, self.env.startJulian,self.env.endJulian	
        # landsat4 image collections 
        if self.env.useL4:        
            landsat4 =  ee.ImageCollection('LANDSAT/LT04/C01/T1_SR').filterDate(startDate,endDate).filterBounds(self.env.location)
            landsat4 = landsat4.filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)
            landsat4 = landsat4.filter(ee.Filter.calendarRange(self.env.startJulian,self.env.endJulian))
	    self.landsat4count  = str(landsat4.size().getInfo()) 
	    if landsat4.size().getInfo() > 0:
		if self.env.defringe == True:
		     landsat4 =  landsat4.map(self.DefringeLandsat)
		if self.env.maskSR == True:
		    landsat4 = landsat4.map(self.CloudMaskSR)
		if not merge:
		    landsatCollection = landsat4.select(self.env.sensorBandDictLandsatSR.get('L4'),self.env.bandNamesLandsat)
		    merge = True
	
        # landsat 5 image collections 
        if self.env.useL5:
            landsat5 =  ee.ImageCollection('LANDSAT/LT05/C01/T1_SR').filterDate(startDate,endDate).filterBounds(self.env.location)
            landsat5 = landsat5.filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)
            landsat5 = landsat5.filter(ee.Filter.calendarRange(self.env.startJulian,self.env.endJulian))
	    self.landsat5count  = str(landsat5.size().getInfo()) 
	    if landsat5.size().getInfo() > 0:
		if self.env.defringe == True:
		     landsat5 =  landsat5.map(self.DefringeLandsat)
		if self.env.maskSR == True:
		    landsat5 = landsat5.map(self.CloudMaskSR)
		if not merge:
		    landsatCollection = landsat5.select(self.env.sensorBandDictLandsatSR.get('L5'),self.env.bandNamesLandsat)
		    merge = True
		else:
		    landsatCollection = landsatCollection.merge(landsat5.select(self.env.sensorBandDictLandsatSR.get('L5'),self.env.bandNamesLandsat))
	
        # landsat 7 image collections  
        l7slm = True
	if self.env.useL7:
            landsat7 =  ee.ImageCollection('LANDSAT/LE07/C01/T1_SR').filterDate(startDate,endDate).filterBounds(self.env.location)
	    if self.env.startYear == 2003 or self.env.endYear == 2003:
		if self.env.useL7scanline == False:
		    landsat7 = landsat7.filterDate(startDate,self.env.l7Failed)
	    if self.env.startYear > 2003 and self.env.useL7scanline == False:
		l7slm = False
	    if l7slm == True:
		landsat7 = landsat7.filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)
		landsat7 = landsat7.filter(ee.Filter.calendarRange(self.env.startJulian,self.env.endJulian))
		self.landsat7count = str(landsat7.size().getInfo()) 
		if landsat7.size().getInfo() > 0:
		    if self.env.defringe == True:
			landsat7 =  landsat7.map(self.DefringeLandsat)            
		    if self.env.maskSR == True:
			landsat7 = landsat7.map(self.CloudMaskSR)
		    if not merge:
			landsatCollection = landsat7.select(self.env.sensorBandDictLandsatSR.get('L7'),self.env.bandNamesLandsat)
			merge = True
		    else:
			landsatCollection = landsatCollection.merge(landsat7.select(self.env.sensorBandDictLandsatSR.get('L7'),self.env.bandNamesLandsat))

        # landsat8  image collections 				        
        if self.env.useL8:
            landsat8 =  ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filterDate(startDate,endDate).filterBounds(self.env.location)
            landsat8 = landsat8.filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)
            landsat8 = landsat8.filter(ee.Filter.calendarRange(self.env.startJulian,self.env.endJulian))
	    self.landsat8count =str(landsat8.size().getInfo())
	    if landsat8.size().getInfo() > 0:
		if self.env.maskSR == True:
		    landsat8 = landsat8.map(self.CloudMaskSRL8)    
		if not merge:
		    landsatCollection = landsat8.select(self.env.sensorBandDictLandsatSR.get('L8'),self.env.bandNamesLandsat)
		    merge = True
		else:
		    landsatCollection = landsatCollection.merge(landsat8.select(self.env.sensorBandDictLandsatSR.get('L8'),self.env.bandNamesLandsat))            
   
        if landsatCollection:
	    count = landsatCollection.size();
	    
	    landsatCollection = landsatCollection.map(self.ScaleLandsat)         
	    
	    # return the image collection           
	    return ee.ImageCollection(landsatCollection)
       

    def returnCollection(self,start,end):
        """Calculate percentiles to filter imagery"""  
        
        logging.info('calculate percentiles')
        
        #startDate = ee.Date.fromYMD(1984,1,1)
        #endDate = ee.Date.fromYMD(2099,1,1)    
        cloudCoverMax = self.env.metadataCloudCoverMax 
        
        # get the images
        collection = self.GetLandsat(start,end,cloudCoverMax)
        collection = collection.map(self.maskClouds)   
	
	return collection
       
    def CloudMaskSR(self,img):
         """apply cf-mask Landsat""" 
	 
	 QA = img.select("pixel_qa")
	 #mask = ee.Image(self.getQABits(QA,3, 5, 'internal_quality_flag')); 
	 #print mask
        
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
    
	logging.info('return scaled image')
        return ee.Image(image.copyProperties(img))

    def reScaleLandsat(self,img):
        """Landast is scaled by factor 0.0001 """
        
	thermalBand = ee.List(['thermal'])
	gapfillBand = ee.List(['gapfill'])
	thermal = ee.Image(img).select(thermalBand).multiply(10)
	gapfill = ee.Image(img).select('gapfill')
	
	otherBands = ee.Image(img).bandNames().removeAll(thermalBand)
	otherBands = otherBands.removeAll(gapfillBand)
        scaled = ee.Image(img).select(otherBands).divide(0.0001)
	
	image = ee.Image(scaled.addBands(thermal).addBands(gapfill)).int16()
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
	    #print bandsUpper

	    percentilesUp = self.percentile.select(bandsUpper,selectedBand  )
	    percentilesLow = self.percentile.select(bandsLower,selectedBand )
	    
	    #print percentilesUp
	    imgToMask = img.select(selectedBand)
			     
	    darkMask = ee.Image(imgToMask.lt(percentilesLow).reduce(ee.Reducer.sum())).eq(0)
	    lightMask = ee.Image(imgToMask.gt(percentilesUp).reduce(ee.Reducer.sum())).eq(0)
	    
	    #img = ee.Image(img.updateMask(lightMask).updateMask(darkMask).copyProperties(img))
	    img = ee.Image(img.updateMask(lightMask).copyProperties(img))
	    
	return img


 

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
        
        outputName = self.env.userID + assetName +self.env.timeString 
        logging.info('export image to asset: ' + str(outputName))   
	
        startDate = ee.Date.fromYMD(self.env.startYear,1,1)
        endDate = ee.Date.fromYMD(self.env.endYear,12,31)    

	image = ee.Image(img).set({'system:time_start':startDate.millis(), \
				    'startyear':self.env.startYear, \
			            'endyear':self.env.endYear, \
			            'startJulian':self.env.startJulian, \
			            'startJulian':self.env.startJulian, \
		                    'endJulian':self.env.endJulian,
				    'source':'USGS_SR',\
				    'median':self.env.calcMedian,\
				    'mean':self.env.calcMean,\
				    'medoid':self.env.calcMedoid,\
				    'count_landsat_4':self.landsat4count,\
				    'count_landsat_5':self.landsat5count,\
				    'count_landsat_7':self.landsat7count,\
				    'count_landsat_8':self.landsat8count,\
				    'cloud_threshold':self.env.metadataCloudCoverMax ,\
				    'defringe':str(self.env.defringe),\
				    'apply_usgs_cloud_filter': str(self.env.maskSR ),\
				    'applyTDOM':self.env.applyTDOM,\
				    'pixel_size':self.env.pixSize,\
				    'zScoreThresh':self.env.zScoreThresh,\
				    'shadowSumThresh':self.env.shadowSumThresh,\
				    'dilatePixels':self.env.dilatePixels,\
				    'filter_percentile':str(self.env.filterPercentile),\
				    'filter_percentile_years':self.env.filterPercentileYears,\
				    'upper_percentile': self.env.highPercentile,\
				    'calculate_indices':str(self.env.getIndices),\
				    'calculate_tasselcap':str(self.env.calcTasselcap),\
				    'calculate_tasselcap_angles':str(self.env.calcTCangles),\
				    'Gap_filling':str(self.env.fillGaps),\
				    'years_of_gap_filling':self.env.fillGapYears,\
				    'version':'1.0'}).clip(self.env.mekongRegion) 


        #task_ordered = ee.batch.Export.image.toAsset(image=ee.Image(img), description=str(self.env.timeString)+assetName, assetId=outputName,region=self.env.location['coordinates'], maxPixels=1e13,scale=self.env.pixSize)
        task_ordered = ee.batch.Export.image.toAsset(image=ee.Image(image), description=assetName, assetId=outputName,region=self.env.location['coordinates'], maxPixels=1e13,scale=self.env.pixSize)
        
        # start task
        task_ordered.start() 
    
    def addIndices(self,img):
	""" Function to add common (and less common) spectral indices to an image.
	    Includes the Normalized Difference Spectral Vector from (Angiuli and Trianni, 2014) """

	#Add Normalized Difference Spectral Vector (NDSV)
        img = img.addBands(img.normalizedDifference(['blue','green']).rename(['ND_blue_green']));
	img = img.addBands(img.normalizedDifference(['blue','red']).rename(['ND_blue_red']));
	img = img.addBands(img.normalizedDifference(['blue','nir']).rename(['ND_blue_nir']));
	img = img.addBands(img.normalizedDifference(['blue','swir1']).rename(['ND_blue_swir1']));
	img = img.addBands(img.normalizedDifference(['blue','swir2']).rename(['ND_blue_swir2']));

	img = img.addBands(img.normalizedDifference(['green','red']).rename(['ND_green_red']));
	img = img.addBands(img.normalizedDifference(['green','nir']).rename(['ND_green_nir']));  # NDWBI
	img = img.addBands(img.normalizedDifference(['green','swir1']).rename(['ND_green_swir1']));  # NDSI, MNDWI
	img = img.addBands(img.normalizedDifference(['green','swir2']).rename(['ND_green_swir2']));

	img = img.addBands(img.normalizedDifference(['red','swir1']).rename(['ND_red_swir1']));
	img = img.addBands(img.normalizedDifference(['red','swir2']).rename(['ND_red_swir2']));

	img = img.addBands(img.normalizedDifference(['nir','red']).rename(['ND_nir_red']));  # NDVI
	img = img.addBands(img.normalizedDifference(['nir','swir1']).rename(['ND_nir_swir1']));  # NDWI, LSWI, -NDBI
	img = img.addBands(img.normalizedDifference(['nir','swir2']).rename(['ND_nir_swir2']));  # NBR, MNDVI

	img = img.addBands(img.normalizedDifference(['swir1','swir2']).rename(['ND_swir1_swir2']));
  
	# Add ratios
	img = img.addBands(img.select('swir1').divide(img.select('nir')).rename(['R_swir1_nir']));  # ratio 5/4
	img = img.addBands(img.select('red').divide(img.select('swir1')).rename(['R_red_swir1']));  # ratio 3/5

	#Add Enhanced Vegetation Index (EVI)
	evi = img.expression(
	    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
	      'NIR': img.select('nir'),
	      'RED': img.select('red'),
	      'BLUE': img.select('blue')
	  }).float();
	
	img = img.addBands(evi.rename(['EVI']));
	  
	# Add Soil Adjust Vegetation Index (SAVI)
	# using L = 0.5;
	savi = img.expression(
	    '(NIR - RED) * (1 + 0.5)/(NIR + RED + 0.5)', {
	      'NIR': img.select('nir'),
	      'RED': img.select('red')
	  }).float();
	img = img.addBands(savi.rename(['SAVI']));
	  
	# Add Index-Based Built-Up Index (IBI)
	ibi_a = img.expression(
	    '2*SWIR1/(SWIR1 + NIR)', {
	      'SWIR1': img.select('swir1'),
	      'NIR': img.select('nir')
	    }).rename(['IBI_A']);
	
	ibi_b = img.expression(
	    '(NIR/(NIR + RED)) + (GREEN/(GREEN + SWIR1))', {
	      'NIR': img.select('nir'),
	      'RED': img.select('red'),
	      'GREEN': img.select('green'),
	      'SWIR1': img.select('swir1')
	    }).rename(['IBI_B']);
	
	ibi_a = ibi_a.addBands(ibi_b);
	ibi = ibi_a.normalizedDifference(['IBI_A','IBI_B']);
	img = img.addBands(ibi.rename(['IBI']));
	  
	return img;       


    def unmaskYears(self,img,year):
	""" Function to unmask nodata withpixels previous year """
	
	print "unmasking for year " + str(self.env.startYear-year) 
	startDate = ee.Date.fromYMD(self.env.startYear-year,1,1)
	endDate = ee.Date.fromYMD(self.env.endYear-year,12,31)    
	prev = self.GetLandsat(startDate,endDate,self.env.metadataCloudCoverMax).select(self.env.exportBands)
	prev = self.maskShadows(prev)
	if prev.size().getInfo() > 0:
	    prev = prev.map(self.MaskPercentile) 
	    previmg = self.medoidMosaic(prev) 
	    
	    #previmg = ee.Image(prev.median())
	    
	    previmg = previmg.mask(previmg.gt(0))
	    gapfilter = previmg.select(["blue"]).gt(0).multiply(self.env.startYear-year)
	    previmg = previmg.addBands(gapfilter.rename(['gapfill']))
	    
	    img = img.unmask(previmg)
	
	return ee.Image(img)

    def unmaskFutureYears(self,img,year):
	""" Function to unmask nodata withpixels future year """
	
	print "unmasking for year " + str(self.env.startYear+year) 
	startDate = ee.Date.fromYMD(self.env.startYear+year,1,1)
	endDate = ee.Date.fromYMD(self.env.endYear+year,12,31)    
	prev = self.GetLandsat(startDate,endDate,self.env.metadataCloudCoverMax).select(self.env.exportBands)
	prev = self.maskShadows(prev)
	prev = prev.map(self.MaskPercentile) 
	if prev.size().getInfo() > 0:
	    prev = prev.map(self.MaskPercentile) 
	    previmg = self.medoidMosaic(prev) 	    
	    previmg = previmg.mask(previmg.gt(0))
	    gapfilter = previmg.select(["blue"]).gt(0).multiply(self.env.startYear-year)
	    previmg = previmg.addBands(gapfilter.rename(['gapfill']))
	    
	    img = img.unmask(previmg)
	return ee.Image(img)

    def makeTiles(self):
	# set bounds
	xmin = 91.979254;
	xmax = 114.664186;
	ymin = 5.429208;
	ymax = 28.728774;

	# number of rows and columns
	n = 2;

	for  i in range(0, n, 1):
	    for j in range(0, n,1):	# x, y distance of one block
		xs = (xmax - xmin) / n
		ys = (ymax - ymin) / n
		
		xl = xmin + i * xs;
		xr = xmin + (i+1) * xs;
		
		yt = ymin + (j*ys);
		yb = ymin + (j+1)*ys;
		geom =   [[xl, yt], [xl, yb], [xr, yb], [xr, yt]];

		col = SurfaceReflectance().RunModel(geom,i,j)

    def getTasseledCap(self,img,bands):
	"""Function to compute the Tasseled Cap transformation and return an image"""

        logging.info('get tasselcap for computed images')
	
	coefficients = ee.Array([
	    [0.3037, 0.2793, 0.4743, 0.5585, 0.5082, 0.1863],
	    [-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800],
	    [0.1509, 0.1973, 0.3279, 0.3406, -0.7112, -0.4572],
	    [-0.8242, 0.0849, 0.4392, -0.0580, 0.2012, -0.2768],
	    [-0.3280, 0.0549, 0.1075, 0.1855, -0.4357, 0.8085],
	    [0.1084, -0.9022, 0.4120, 0.0573, -0.0251, 0.0238]
	  ]);
	
	# Make an Array Image, with a 1-D Array per pixel.
	arrayImage1D = img.select(bands).toArray()
	
	# Make an Array Image with a 2-D Array per pixel, 6x1.
	arrayImage2D = arrayImage1D.toArray(1)
	
	componentsImage = ee.Image(coefficients).matrixMultiply(arrayImage2D).arrayProject([0]).arrayFlatten([['brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth']]).float();
  
	# Get a multi-band image with TC-named bands.
  	return img.addBands(componentsImage);


    def addTCAngles(self,img):
	""" Function to add Tasseled Cap angles and distances to an image.
	    Assumes image has bands: 'brightness', 'greenness', and 'wetness'."""
	
	logging.info('add tasseled cap angles')
	
	# Select brightness, greenness, and wetness bands	
	brightness = img.select('brightness');
	greenness = img.select('greenness');
	wetness = img.select('wetness');
  
	# Calculate Tasseled Cap angles and distances
	tcAngleBG = brightness.atan2(greenness).divide(math.pi).rename(['tcAngleBG']);
	tcAngleGW = greenness.atan2(wetness).divide(math.pi).rename(['tcAngleGW']);
	tcAngleBW = brightness.atan2(wetness).divide(math.pi).rename(['tcAngleBW']);
	tcDistBG = brightness.hypot(greenness).rename(['tcDistBG']);
	tcDistGW = greenness.hypot(wetness).rename(['tcDistGW']);
	tcDistBW = brightness.hypot(wetness).rename(['tcDistBW']);
	
	img = img.addBands(tcAngleBG).addBands(tcAngleGW).addBands(tcAngleBW).addBands(tcDistBG).addBands(tcDistGW).addBands(tcDistBW);
	return img;

	
    def medoidMosaic(self,collection):
	""" medoid composite with equal weight among indices """
	
	f = ee.Image(collection.first());
	bandNames = f.bandNames();
	bandNumbers = ee.List.sequence(1,bandNames.length());
	median = ee.ImageCollection(collection).median();
	
	def subtractmedian(img):
	    diff = ee.Image(img).subtract(median).pow(ee.Image.constant(2));
	    return diff.reduce('sum').addBands(img);
	
	medoid = collection.map(subtractmedian)
  
	medoid = ee.ImageCollection(medoid).reduce(ee.Reducer.min(bandNames.length().add(1))).select(bandNumbers,bandNames);
  
	return medoid;



class Primitives():
 
    def __init__(self):
        """Initialize the Surfrace Reflectance app."""  
        
        # import the log library
        import logging
	
	# get the environment
        self.env = environment()    
        
   
if __name__ == "__main__":
  
    # set argument parsing object
    parser = argparse.ArgumentParser(description="Create Landsat image composites using Google\
                                                  Earth Engine.")
   
    parser.add_argument('--year','-y', type=str,required=True,
                        help="Year to perform the ats correction and save to asset format in 'YYYY'")

    parser.add_argument('--season','-s', choices=['drycool','dryhot','rainy'],type=str,
                        help="Season to create composite for, these align with SERVIR-Mekong's seasonal composite times")

    parser.add_argument('--user','-u', type=str, default="servir-mekong",choices=['servir-mekong','servirmekong',"ate","biplov","quyen"],
			help="specify user account to run task")

    args = parser.parse_args() # get arguments  
  
    # user account to run task on
    userName = args.user
    year = args.year
    seasons = args.season
    
    # create a new file in ~/.config/earthengine/credentials with token of user
    addUserCredentials(userName)
    geom = ''    
    #SurfaceReflectance().RunModel(geom,1,1)
    SurfaceReflectance().makeTiles()
