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

--user (-u)      : user account used to create the composite
                 : changes the ~/.config/earthengine/credentials file
                 : dictionary is called to get credentials
                 : options are servirmekong, servir-mekong, ate, biplov .. default is servir-mekong

Example Usage
-------------

1) export surface reflectance composite for dryhot season of 2000 to assets:

  $ python model.py -y 2000 -u Quyen

"""

import ee
import logging
import time
import math
from usercredentials import addUserCredentials
import argparse
import numpy as np

class environment(object):
        
	def __init__(self):
		"""Initialize the environment."""   
         
        # Initialize the Earth Engine object, using the authentication credentials.
        #ee.Initialize()
		self.timeString = time.strftime("%Y%m%d_%H%M%S")
		
		self.year = 0
		
		self.userID = "projects/servir-mekong/" 
		self.userID = "users/servirmekong/" 
		self.pixSize = 30
		
		self.nModels = 200
		
		self.year = 0

		# Load Mekong study region (Myanmar, Thailand, Laos, Vietnam, Cambodia)
		mekongBuffer = ee.FeatureCollection('ft:1LEGeqwlBCAlN61ie5ol24NdUDqB1MgpFR_sJNWQJ');
		mekongRegion = mekongBuffer.geometry();
		self.studyArea = mekongRegion;
		#self.studyArea = ee.Geometry.Polygon([[103.876,18.552],[105.806,18.552],[105.806,19.999],[103.876,19.999],[103.876,18.552]])
		
		#self.studyArea = ee.Geometry.Polygon([[105.260,21.553],[105.260,20.6],[106.556,20.6],[106.446,21.527],[105.260,21.553]])


		self.composite = "Medoid"

		self.classFieldName = 'land_class';
		
		self.modelType = 'RF'
		#self.modelType = 'SVM'
		#self.modelType = 'per'
		#self.modelType = 'nav'
		#self.modelType = 'cart'
		
		self.exportName = ''#  + self.modelType + self.composite + str(self.nModels)
		#self.exportName = 'wetlands_'  + self.modelType + self.composite + str(self.nModels)
		self.assetName = ''#  + self.modelType + self.composite + str(self.nModels)
		#self.assetName = 'wetlands_'  + self.modelType + self.composite + str(self.nModels)

		self.trainingData = "1M71zwCyPyqEEuq8IBrWewONoBVDWHl22D01siCV_"


class indices():

	def __init__(self):
		
		# list with functions to call for each index
		self.functionList = {"ND_blue_green" : self.ND_blue_green, \
							 "ND_blue_red" : self.ND_blue_red, \
							 "ND_blue_nir" : self.ND_blue_nir, \
							 "ND_blue_swir1" : self.ND_blue_swir1, \
							 "ND_blue_swir2" : self.ND_blue_swir2, \
							 "ND_green_red" : self.ND_green_red, \
							 "ND_green_nir" : self.ND_green_nir, \
							 "ND_green_swir1" : self.ND_green_swir1, \
							 "ND_green_swir2" : self.ND_green_swir2, \
							 "ND_red_swir1" : self.ND_red_swir1, \
							 "ND_red_swir2" : self.ND_red_swir2, \
							 "ND_nir_red" : self.ND_nir_red, \
							 "ND_nir_swir1" : self.ND_nir_swir1, \
							 "ND_nir_swir2" : self.ND_nir_swir2, \
							 "ND_swir1_swir2" : self.ND_swir1_swir2, \
							 "R_swir1_nir" : self.R_swir1_nir, \
							 "R_red_swir1" : self.R_red_swir1, \
							 "EVI" : self.EVI, \
							 "SAVI" : self.SAVI, \
							 "IBI" : self.IBI}


	def addAllTasselCapIndices(self,img): 
		""" Function to get all tasselCap indices """
		
		def getTasseledCap(img):
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
		
			bands=ee.List(['blue','green','red','nir','swir1','swir2'])
			
			# Make an Array Image, with a 1-D Array per pixel.
			arrayImage1D = img.select(bands).toArray()
		
			# Make an Array Image with a 2-D Array per pixel, 6x1.
			arrayImage2D = arrayImage1D.toArray(1)
		
			componentsImage = ee.Image(coefficients).matrixMultiply(arrayImage2D).arrayProject([0]).arrayFlatten([['brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth']]).float();
	  
			# Get a multi-band image with TC-named bands.
			return img.addBands(componentsImage);	
			
			
		def addTCAngles(img):

			""" Function to add Tasseled Cap angles and distances to an image. Assumes image has bands: 'brightness', 'greenness', and 'wetness'."""
		
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
	
	
		img = getTasseledCap(img)
		img = addTCAngles(img)
		return img

	def ND_blue_green(self,img):
		img = img.addBands(img.normalizedDifference(['blue','green']).rename(['ND_blue_green']));
		return img
	
	def ND_blue_red(self,img):
		img = img.addBands(img.normalizedDifference(['blue','red']).rename(['ND_blue_red']));
		return img
	
	def ND_blue_nir(self,img):
		img = img.addBands(img.normalizedDifference(['blue','nir']).rename(['ND_blue_nir']));
		return img
	
	def ND_blue_swir1(self,img):
		img = img.addBands(img.normalizedDifference(['blue','swir1']).rename(['ND_blue_swir1']));
		return img
	
	def ND_blue_swir2(self,img):
		img = img.addBands(img.normalizedDifference(['blue','swir2']).rename(['ND_blue_swir2']));
		return img

	def ND_green_red(self,img):
		img = img.addBands(img.normalizedDifference(['green','red']).rename(['ND_green_red']));
		return img
	
	def ND_green_nir(self,img):
		img = img.addBands(img.normalizedDifference(['green','nir']).rename(['ND_green_nir']));  # NDWBI
		return img
	
	def ND_green_swir1(self,img):
		img = img.addBands(img.normalizedDifference(['green','swir1']).rename(['ND_green_swir1']));  # NDSI, MNDWI
		return img
	
	def ND_green_swir2(self,img):
		img = img.addBands(img.normalizedDifference(['green','swir2']).rename(['ND_green_swir2']));
		return img
		
	def ND_red_swir1(self,img):
		img = img.addBands(img.normalizedDifference(['red','swir1']).rename(['ND_red_swir1']));
		return img
			
	def ND_red_swir2(self,img):
		img = img.addBands(img.normalizedDifference(['red','swir2']).rename(['ND_red_swir2']));
		return img

	def ND_nir_red(self,img):
		img = img.addBands(img.normalizedDifference(['nir','red']).rename(['ND_nir_red']));  # NDVI
		return img
	
	def ND_nir_swir1(self,img):
		img = img.addBands(img.normalizedDifference(['nir','swir1']).rename(['ND_nir_swir1']));  # NDWI, LSWI, -NDBI
		return img
	
	def ND_nir_swir2(self,img):
		img = img.addBands(img.normalizedDifference(['nir','swir2']).rename(['ND_nir_swir2']));  # NBR, MNDVI
		return img

	def ND_swir1_swir2(self,img):
		img = img.addBands(img.normalizedDifference(['swir1','swir2']).rename(['ND_swir1_swir2']));
		return img
  
	def R_swir1_nir(self,img):
		# Add ratios
		img = img.addBands(img.select('swir1').divide(img.select('nir')).rename(['R_swir1_nir']));  # ratio 5/4
		return img
			
	def R_red_swir1(self,img):
		img = img.addBands(img.select('red').divide(img.select('swir1')).rename(['R_red_swir1']));  # ratio 3/5
		return img

	def EVI(self,img):
		#Add Enhanced Vegetation Index (EVI)
		evi = img.expression(
			'2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
			  'NIR': img.select('nir'),
			  'RED': img.select('red'),
			  'BLUE': img.select('blue')
		  }).float();
	
		img = img.addBands(evi.rename(['EVI']));

		return img
	  
	def SAVI(self,img):
		# Add Soil Adjust Vegetation Index (SAVI)
		# using L = 0.5;
		savi = img.expression(
			'(NIR - RED) * (1 + 0.5)/(NIR + RED + 0.5)', {
			  'NIR': img.select('nir'),
			  'RED': img.select('red')
		  }).float();
		img = img.addBands(savi.rename(['SAVI']));

		return img
	  
	def IBI(self,img):
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
		
		return img


class primitives():
	
	def __init__(self):
		"""Initialize the Surfrace Reflectance app."""  
        
		# import the log library
		import logging
	
		# get the environment
		self.env = environment()
		
		# get object with indices
		self.indices = indices() 
		
	
	def importData(self):
		print "import data"

	def getIndices(self,img,covariates):	
		""" add indices to image"""
		
		# no need to add indices that are already there
		indices = self.removeDuplicates(covariates,img.bandNames().getInfo())
		
		for item in indices:
			img = self.indices.functionList[item](img)

		return img
		
	def createPrimitive(self,drycool,dryhot,rainy,y,primitive ):
		""" calculate the primitive """ 


		if primitive == "wetlands":
			self.env.exportName = 'wetlands'  + self.env.modelType + self.env.composite + str(self.env.nModels)
			self.env.assetName =  'wetlands_' + str(y) 
			self.env.userID = "projects/servir-mekong/binary_primitives/wetlands/" 
			selectedBands = ["change_abs","distCoast","drycool_blue","drycool_ND_blue_nir","drycool_ND_green_red","drycool_ND_swir1_swir2","drycool_nir","drycool_R_red_swir1",\
							 "drycool_sixth","drycool_tcAngleGW","drycool_tcDistBG","dryhot_blue","dryhot_ND_green_red","dryhot_ND_green_swir2","dryhot_ND_red_swir1","dryhot_ND_red_swir2",\
							 "dryhot_nir","dryhot_R_red_swir1","dryhot_tcAngleGW","dryhot_tcDistBG","dryhot_tcDistGW","elevation","occurrence","rainy_blue","rainy_fifth","rainy_ND_red_swir2",\
							 "rainy_nir","rainy_R_red_swir1","rainy_tcAngleGW","rainy_tcDistGW","slope","tcc","treeheight"] 
			trainingDataSet = ee.FeatureCollection("ft:1RXwtefKT7jKNV6IAWZSdRDx1q5x5RdQm9XVo-_y-")
			classNames = ee.List(['wetlands',"other"]);
		
		if primitive == "mangroves":
			self.env.exportName = 'mangroves'  + self.env.modelType + self.env.composite + str(self.env.nModels)
			self.env.assetName =  'mangroves_' + str(y)
			self.env.userID = "projects/servir-mekong/binary_primitives/mangroves/" 
			selectedBands = ["distCoast","drycool_blue","drycool_EVI","drycool_ND_blue_nir","drycool_ND_green_red","drycool_ND_swir1_swir2","drycool_sixth","drycool_tcAngleGW",\
							 "dryhot_blue","dryhot_ND_blue_nir","dryhot_ND_green_swir2","dryhot_sixth","dryhot_tcAngleGW","dryhot_tcDistBG","rainy_IBI","rainy_ND_blue_nir", \
							 "rainy_ND_swir1_swir2","rainy_tcAngleBW","rainy_tcAngleGW","tcc","treeheight"]
			trainingDataSet = ee.FeatureCollection("ft:1Ud9HdqsG-_MvatTf5PI-4h8y-SkRSWqt6_UO6TzP")
			classNames = ee.List(['mangroves',"other"]);
			
		if primitive == "barren":
			self.env.exportName = 'barren_imperv_crop_rice'  + self.env.modelType + self.env.composite + str(self.env.nModels)
			self.env.assetName =  'barren_imperv_crop_rice' + str(y) 
			self.env.userID = "projects/servir-mekong/SR_primitives/barren_imperv_crop_rice/" 
			self.env.userID = "users/servirmekong/temp/" 
			#trainingDataSet = ee.FeatureCollection("ft:14AwPqA8ADQU5kCdPN-J2HdD8qUQlREWH97fajYeZ")
			#trainingDataSet = ee.FeatureCollection("ft:1fYKlAs2Xba3HT66FTutrehZROXOMQVuGMkTvJo0c") # allyears
			trainingDataSet = ee.FeatureCollection("ft:10C65Hh-VE4K7IrzUISygmI4pX4HGi98MDKJ9lIzM") # mandaylay updated
			classNames = ee.List(['barren','imperv','crop','rice',"other"]);

		if primitive == "imperv":
			self.env.exportName = 'imperv'  + self.env.modelType + self.env.composite + str(self.env.nModels)
			self.env.assetName =  'imperv' + str(y) 
			selectedBands =  ["drycool_fourth","drycool_green","drycool_ND_blue_green","drycool_ND_blue_nir","drycool_ND_blue_swir1","drycool_ND_swir1_swir2","drycool_nir","drycool_tcAngleGW",\
							  "drycool_tcDistBG","drycool_tcDistGW","dryhot_fourth","dryhot_ND_blue_nir","dryhot_ND_blue_swir1","dryhot_nir","dryhot_tcDistBG","dryhot_tcDistGW","dryhot_wetness",\
							  "rainy_fourth","rainy_ND_blue_green","rainy_ND_green_red","rainy_ND_green_swir2","rainy_nir","rainy_R_swir1_nir","rainy_tcAngleBW","rainy_tcAngleGW","rainy_tcDistBG","rainy_tcDistGW"]
			self.env.userID = "projects/servir-mekong/binary_primitives/imperv/" 
			trainingDataSet = ee.FeatureCollection("ft:16FHVOH2RVLyM3jzBmpy0g7GHdLQ3t9rHsg3oKgD7") 
			classNames = ee.List(['imperv',"other"]);



			
		if primitive == "grass":
			self.env.exportName = 'grass'  + self.env.modelType + self.env.composite + str(self.env.nModels)
			self.env.assetName =  'grass_' + str(y) 
			self.env.userID = "projects/servir-mekong/SR_primitives/grass/" 
			selectedBands = ["distCoast","drycool_EVI","drycool_ND_blue_green","drycool_ND_blue_swir1","drycool_ND_blue_swir2","drycool_ND_red_swir1","drycool_R_red_swir1","drycool_tcAngleGW",\
							 "drycool_thermal","dryhot_blue","dryhot_EVI","dryhot_ND_green_red","dryhot_ND_red_swir2","dryhot_R_red_swir1","dryhot_thermal","dryhot_wetness","elevation", \
							 "rainy_R_red_swir1","slope","tcc","treeheight"]
			trainingDataSet = ee.FeatureCollection("ft:1ajeo9UVM9o7fGMWoARfMNCxXVzrULSvdSFGsv6Ov")
			classNames = ee.List(['grass',"other"]);			
			
		if primitive == "cropplantation":
			self.env.exportName = 'cropplantation'  + self.env.modelType + self.env.composite + str(self.env.nModels)
			self.env.assetName =  'cropplantation_' + str(y) 
			selectedBands = ["distCoast","drycool_EVI","drycool_fourth","drycool_IBI","drycool_ND_blue_nir","drycool_ND_green_red","drycool_ND_red_swir2","drycool_nir","drycool_R_red_swir1",\
							 "drycool_sixth","drycool_tcAngleBW","drycool_tcAngleGW","drycool_tcDistBG","drycool_thermal","dryhot_ND_blue_swir2","dryhot_ND_red_swir2","dryhot_R_red_swir1",\
							 "dryhot_tcAngleGW","dryhot_tcDistBG","elevation","rainy_blue","rainy_EVI","rainy_ND_blue_nir","rainy_ND_blue_swir2","rainy_ND_red_swir2","rainy_nir","rainy_R_red_swir1",\
							 "rainy_tcAngleBW","rainy_tcAngleGW","rainy_tcDistBG","rainy_tcDistGW","slope","tcc","treeheight"]
			self.env.userID = "projects/servir-mekong/SR_primitives/cropplantation/" 
			trainingDataSet = ee.FeatureCollection("ft:1uD_uqVQc_dTS0ZcBHfsyoeLrQtmoBcQOZOsJFFt0")
			classNames = ee.List(['cropplantation',"other"]);				
		
		
		if primitive == "irrigated":
			self.env.exportName = 'irrigated'  + self.env.modelType + self.env.composite + str(self.env.nModels)
			self.env.assetName =  'irrigated_' + str(y) 
			self.env.userID = "users/servirmekong/temp/" 
			
			selectedBands = ["auto","distCoast","drycool_blue","drycool_fourth","drycool_greenness","drycool_ND_blue_green","drycool_ND_blue_swir1","drycool_ND_green_red","drycool_ND_red_swir2",\
							 "drycool_ND_swir1_swir2","drycool_nir","drycool_R_swir1_nir","drycool_sixth","drycool_tcAngleGW","drycool_tcDistGW","drycool_thermal","dryhot_ND_blue_swir1",\
							 "dryhot_ND_green_swir2","dryhot_ND_red_swir2","dryhot_nir","dryhot_tcAngleGW","dryhot_tcDistBG","dryhot_tcDistGW","dryhot_thermal","elevation","R2_cycle1","R2_cycle2",\
							 "R2_cycle3","rainy_fifth","rainy_ND_green_swir2","rainy_ND_red_swir2","rainy_ND_swir1_swir2","rainy_R_swir1_nir","rainy_sixth","rainy_swir2","rainy_tcAngleGW"]
			trainingDataSet = ee.FeatureCollection("ft:1-lkLTrW_ea6eGtnqe5pyQV_6vmPBO864-OWDA-LO")
			classNames = ee.List(['irrigated',"other"]);	

		if primitive == "rice":
			self.env.exportName = 'rice'  + self.env.modelType + self.env.composite + str(self.env.nModels)
			self.env.assetName =  'rice_' + str(y) 
			self.env.userID = "projects/servir-mekong/binary_primitives/rice/" 
			
			selectedBands = ["distCoast","drycool_blue","drycool_ND_swir1_swir2","drycool_sixth","dryhot_EVI","dryhot_ND_blue_nir","dryhot_ND_blue_swir1","dryhot_ND_green_swir2","dryhot_tcAngleGW",\
							 "dryhot_tcDistBG","dryhot_thermal","elevation","rainy_blue","rainy_ND_blue_green","rainy_ND_blue_nir","rainy_ND_green_red","rainy_ND_green_swir2","rainy_tcAngleBW",\
							 "rainy_tcAngleGW","rainy_thermal","treeheight"]
			trainingDataSet = ee.FeatureCollection("ft:1gQxjd15imuF5lK_9efRCrh4E7pQ7vhd9fMjoGSwi")
			classNames = ee.List(['rice',"other"]);	

		if primitive == "barren":
			self.env.exportName = 'barren'  + self.env.modelType + self.env.composite + str(self.env.nModels)
			self.env.assetName =  'barren_' + str(y) 
			self.env.userID = "projects/servir-mekong/binary_primitives/barren/" 
			
			selectedBands = ["distCoast","drycool_blue","drycool_EVI","drycool_fourth","drycool_ND_blue_nir","drycool_ND_blue_red","drycool_ND_green_red","drycool_ND_swir1_swir2","drycool_R_red_swir1",\
							 "drycool_R_swir1_nir","drycool_sixth","drycool_tcAngleBW","drycool_tcAngleGW","drycool_tcDistBG","drycool_tcDistGW","dryhot_ND_blue_nir","dryhot_ND_blue_swir1","dryhot_ND_green_swir2",\
							 "dryhot_ND_red_swir1","dryhot_ND_red_swir2","dryhot_nir","dryhot_R_red_swir1","dryhot_tcDistBG","dryhot_tcDistGW","dryhot_thermal_stdDev","elevation","rainy_fourth","rainy_ND_blue_green",\
							 "rainy_ND_blue_nir","rainy_ND_blue_red","rainy_R_red_swir1","rainy_tcAngleBW","rainy_tcAngleGW","slope","tcc","treeheight"]
			trainingDataSet = ee.FeatureCollection("ft:12a03BLHsJ9wsDC8DvySPlbNGuA-jyBcpNNeqE5S4")
			classNames = ee.List(['barren',"other"]);	


		self.env.trainingData = trainingDataSet
		self.env.startYear = y		

		covariates = ["ND_blue_green","ND_blue_red","ND_blue_nir","ND_blue_swir1","ND_blue_swir2","ND_green_red","ND_green_nir","ND_green_swir1","ND_green_swir2","ND_red_swir1","ND_red_swir2","ND_nir_red","ND_nir_swir1","ND_nir_swir2","ND_swir1_swir2","R_swir1_nir","R_red_swir1","EVI","SAVI","IBI"]

		dryhot = self.ScaleBands(dryhot)
		dryhot = self.indices.addAllTasselCapIndices(dryhot)
		dryhot = self.getIndices(dryhot,covariates)
	
		drycool = self.ScaleBands(drycool)
		drycool = self.indices.addAllTasselCapIndices(drycool)
		drycool = self.getIndices(drycool,covariates)
	
		rainy = self.ScaleBands(rainy)
		rainy = self.indices.addAllTasselCapIndices(rainy)
		rainy = self.getIndices(rainy,covariates)
		
		# Get the JRC water 
		water = ee.Image('JRC/GSW1_0/GlobalSurfaceWater').mask(ee.Image(1));
		
		# rename the bands

		dryhot = self.renameImageBands(dryhot,"dryhot") 
		drycool = self.renameImageBands(drycool,"drycool")
		rainy = self.renameImageBands(rainy,"rainy") 
		
		# construct the composites
		composite = drycool.addBands(dryhot).addBands(rainy).addBands(water);
		composite = self.addTopography(composite);
		composite = self.addJRC(composite)
		#composite = self.addNightLights(composite,y)
		
		#auto = ee.Image("users/servirmekong/autocor/autocor_year_2015").select("autocorrelation")
		#composite = composite.addBands(auto.rename(["auto"]))
		
		#cycle = ee.Image("users/servirmekong/seasons/seasons_year_2015")
		#composite = composite.addBands(cycle)
		
		distCoast = ee.Image('projects/servir-mekong/Primitives/DistancetoCoast_1k').float().rename(['distCoast']);
		composite = composite.addBands(distCoast)
		
		''' selectedBands = ['drycool_blue', 'drycool_green', 'drycool_red', 'drycool_nir', 'drycool_swir1', 'drycool_swir2', 'drycool_blue_stdDev', 'drycool_green_stdDev', 'drycool_red_stdDev', 'drycool_nir_stdDev', \
						 'drycool_swir1_stdDev', 'drycool_swir2_stdDev', 'drycool_ND_nir_swir2_stdDev', 'drycool_ND_green_swir1_stdDev', 'drycool_ND_nir_red_stdDev','drycool_thermal', \
						 'drycool_thermal_stdDev', 'drycool_brightness', 'drycool_greenness', 'drycool_wetness', 'drycool_fourth', 'drycool_fifth', 'drycool_sixth', 'drycool_tcAngleBG', \
						 'drycool_tcAngleGW', 'drycool_tcAngleBW', 'drycool_tcDistBG', 'drycool_tcDistGW', 'drycool_tcDistBW', 'drycool_ND_blue_green', 'drycool_ND_blue_red', 'drycool_ND_blue_nir', \
						 'drycool_ND_blue_swir1', 'drycool_ND_blue_swir2', 'drycool_ND_green_red', 'drycool_ND_green_nir', 'drycool_ND_green_swir1', 'drycool_ND_green_swir2', 'drycool_ND_red_swir1', \
						 'drycool_ND_red_swir2', 'drycool_ND_nir_red', 'drycool_ND_nir_swir1', 'drycool_ND_nir_swir2', 'drycool_ND_swir1_swir2', 'drycool_R_swir1_nir', 'drycool_R_red_swir1', 'drycool_EVI', \
						 'drycool_SAVI', 'drycool_IBI', 'dryhot_blue', 'dryhot_green', 'dryhot_red', 'dryhot_nir', 'dryhot_swir1', 'dryhot_swir2', 'dryhot_blue_stdDev', 'dryhot_green_stdDev', 'dryhot_red_stdDev', \
						 'dryhot_nir_stdDev', 'dryhot_swir1_stdDev', 'dryhot_swir2_stdDev', 'dryhot_ND_nir_swir2_stdDev', 'dryhot_ND_green_swir1_stdDev', 'dryhot_ND_nir_red_stdDev', 'dryhot_thermal', 'dryhot_thermal_stdDev', \
						 'dryhot_brightness', 'dryhot_greenness', 'dryhot_wetness', 'dryhot_fourth', 'dryhot_fifth', 'dryhot_sixth', \
						 'dryhot_tcAngleBG', 'dryhot_tcAngleGW', 'dryhot_tcAngleBW', 'dryhot_tcDistBG', 'dryhot_tcDistGW', 'dryhot_tcDistBW', 'dryhot_ND_blue_green', 'dryhot_ND_blue_red', 'dryhot_ND_blue_nir', \
						 'dryhot_ND_blue_swir1', 'dryhot_ND_blue_swir2', 'dryhot_ND_green_red', 'dryhot_ND_green_nir', 'dryhot_ND_green_swir1', 'dryhot_ND_green_swir2', 'dryhot_ND_red_swir1', \
						 'dryhot_ND_red_swir2', 'dryhot_ND_nir_red', 'dryhot_ND_nir_swir1', 'dryhot_ND_nir_swir2', 'dryhot_ND_swir1_swir2', 'dryhot_R_swir1_nir', 'dryhot_R_red_swir1', 'dryhot_EVI', \
						 'dryhot_SAVI', 'dryhot_IBI', 'rainy_blue', 'rainy_green', 'rainy_red', 'rainy_nir', 'rainy_swir1', 'rainy_swir2', 'rainy_blue_stdDev', 'rainy_green_stdDev', 'rainy_red_stdDev', \
						 'rainy_nir_stdDev', 'rainy_swir1_stdDev', 'rainy_swir2_stdDev', 'rainy_ND_nir_swir2_stdDev', 'rainy_ND_green_swir1_stdDev', 'rainy_ND_nir_red_stdDev', 'rainy_thermal', \
						 'rainy_thermal_stdDev','rainy_brightness', 'rainy_greenness', 'rainy_wetness', 'rainy_fourth', 'rainy_fifth', 'rainy_sixth', 'rainy_tcAngleBG', 'rainy_tcAngleGW', \
						 'rainy_tcAngleBW', 'rainy_tcDistBG', 'rainy_tcDistGW', 'rainy_tcDistBW', 'rainy_ND_blue_green', 'rainy_ND_blue_red', 'rainy_ND_blue_nir', 'rainy_ND_blue_swir1', 'rainy_ND_blue_swir2', \
						 'rainy_ND_green_red', 'rainy_ND_green_nir', 'rainy_ND_green_swir1', 'rainy_ND_green_swir2', 'rainy_ND_red_swir1', 'rainy_ND_red_swir2', 'rainy_ND_nir_red', 'rainy_ND_nir_swir1', \
						 'rainy_ND_nir_swir2', 'rainy_ND_swir1_swir2', 'rainy_R_swir1_nir', 'rainy_R_red_swir1', 'rainy_EVI', 'rainy_SAVI', 'rainy_IBI', 'occurrence', 'change_abs', 'change_norm', 'seasonality', \
						 'transition', 'max_extent', 'elevation', 'slope', 'aspect', 'eastness', 'northness','distCoast'] '''
		start = ee.Date.fromYMD(y, 1, 1)
		end  = ee.Date.fromYMD(y, 12,31)
		tcc = ee.Image(ee.ImageCollection("projects/servir-mekong/Primitives/P_canopy").filterDate(start,end).first()).rename(["tcc"])
		composite = composite.addBands(tcc)
		treeheight = ee.Image(ee.ImageCollection("projects/servir-mekong/Primitives/P_tree_height").filterDate(start,end).first()).rename(["treeheight"])
		composite = composite.addBands(treeheight)
		
		composite = composite.select(selectedBands)
				
		#print composite.bandNames().getInfo()
		
		if primitive == "mangroves":
			mask  = composite.select('distCoast').lt(50)
			composite = composite.updateMask(mask)
			
			
			
		#classNames = ee.List(['wetlands',"other"]);

		# run the model		
		"""classification = self.runRandomForestModel(composite, \
						     trainingDataSet, \
						     composite.bandNames(), \
						     self.env.nModels, \
						     self.env.classFieldName, \
						     classNames, \
						     self.env.modelType);"""

		classification = self.runRandomForestModel(composite, \
						     trainingDataSet, \
						     composite.bandNames(), \
						     self.env.nModels, \
						     self.env.classFieldName, \
						     classNames, \
						     self.env.modelType);

		# export the classification
		
		classification = classification

		self.ExportToAsset(self.env.exportName,classification)
		
	def addTopography(self,img):
		"""  Function to add 30m SRTM elevation and derived slope, aspect, eastness, and 
		northness to an image. Elevation is in meters, slope is between 0 and 90 deg,
		aspect is between 0 and 359 deg. Eastness and northness are unitless and are
		between -1 and 1. """

		# Import SRTM elevation data
		elevation = ee.Image("USGS/SRTMGL1_003");
		
		# Calculate slope, aspect, and hillshade
		topo = ee.Algorithms.Terrain(elevation);
		
		# From aspect (a), calculate eastness (sin a), northness (cos a)
		deg2rad = ee.Number(math.pi).divide(180);
		aspect = topo.select(['aspect']);
		aspect_rad = aspect.multiply(deg2rad);
		eastness = aspect_rad.sin().rename(['eastness']).float();
		northness = aspect_rad.cos().rename(['northness']).float();
		
		# Add topography bands to image
		topo = topo.select(['elevation','slope','aspect']).addBands(eastness).addBands(northness);
		img = img.addBands(topo);
		return img;

	def addJRC(self,img):
		""" Function to add JRC Water layers: 'occurrence', 'change_abs', 
			'change_norm', 'seasonality','transition', 'max_extent' """
		
		jrcImage = ee.Image("JRC/GSW1_0/GlobalSurfaceWater")
		
		img = img.addBands(jrcImage.select(['occurrence']).rename(['occurrence']))
		img = img.addBands(jrcImage.select(['change_abs']).rename(['change_abs']))
		img = img.addBands(jrcImage.select(['change_norm']).rename(['change_norm']))
		img = img.addBands(jrcImage.select(['seasonality']).rename(['seasonality']))
		img = img.addBands(jrcImage.select(['transition']).rename(['transition']))
		img = img.addBands(jrcImage.select(['max_extent']).rename(['max_extent']))
		
		return img
		
	def addNightLights(self,img,y):
		""" Function to add nighlights to the composite' """
		
		startDate = ee.Date.fromYMD(y-2, 1, 1)
		endDate = ee.Date.fromYMD(y-2, 12, 31)
		
		if y < 2012:
		
			nightLights = ee.Image(ee.ImageCollection("NOAA/DMSP-OLS/NIGHTTIME_LIGHTS").filterDate(startDate,endDate).mean())	
			img = img.addBands(nightLights.select(["stable_lights"]).rename(["stable_lights"]))
		
		if y >= 2012:
			nightLights = ee.Image(ee.ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG").filterDate(startDate,endDate).mean())	
			img = img.addBands(nightLights.select(["avg_rad"]).rename(["stable_lights"]))
		
		return img


	def runSVMModel(self,image,data,bands,nModels,classFieldName,classNames,modelType):

		#////////////////////////////////////////////////////////////////////////////////
		#// Function to rescale an image such that each band is has a mean of 0 and a 
		#// standard deviation of 1. The function uses a dictionary with the means and 
		#// standard deviations for each band.
		def standardizeImage(image,dict):
			# Get band names
			bands = image.bandNames();
		
			def getMean(band):
				band_mean = ee.String(band).cat('_mean');
				return dict.get(band_mean);	
		
			def getstdDev(band):
				# Convert mean and stdDev values into rasters
				band_stdDev = ee.String(band).cat('_stdDev');
				return dict.get(band_stdDev);   
			
			image_mean = bands.map(getMean)
			image_stdDev = bands.map(getstdDev)
		
			image_mean = ee.Image.constant(image_mean).rename(bands);
	  
			image_stdDev = ee.Image.constant(image_stdDev).rename(bands);
	  
			scaled_image = image.subtract(image_mean).divide(image_stdDev);
	  
			return scaled_image;

		#////////////////////////////////////////////////////////////////////////////////
		#// Function to standardize a feature using a dictionary of means and standard
		#// deviations for each property in the feature.
		def standardizeFeatures(feature):
			feature = ee.Feature(feature);
	  
			# Get properties and save class value for later
			properties = feature.propertyNames().removeAll(['system:index',self.env.classFieldName]);
			classVal = feature.get(self.env.classFieldName);
	  
			#Get list of standardized values
			def getStandardizedValues(prop):
				val = ee.Number(feature.get(prop));
				property_mean = ee.String(prop).cat('_mean');
				property_stdDev = ee.String(prop).cat('_stdDev');
				mean_val = dictionary.get(property_mean);
				stdDev_val = dictionary.get(property_stdDev);
				newval = (val.subtract(mean_val)).divide(stdDev_val);
				return newval;
	  
			newvals = properties.map(getStandardizedValues)
	  
			# Update feature
			updateDict = ee.Dictionary.fromLists(properties,newvals);
			feature = feature.setMulti(updateDict);
			feature = feature.set(self.env.classFieldName,classVal);
			return feature;

		dictionary = image.reduceRegion(reducer = ee.Reducer.mean().combine(ee.Reducer.stdDev(),None,True), \
										geometry=self.env.studyArea, \
										scale= 100, \
										maxPixels= 1e13);
			
		#image = standardizeImage(image,dictionary);
								  
		training = data.select(training_bands.add(self.env.classFieldName));
		
		#print training.bandNames().getInfo()	
	
		training = training.map(standardizeFeatures)
		
			
		classifier = ee.Classifier.svm(kernelType='RBF').train(training,self.env.classFieldName,training_bands);
		classification = image.classify(classifier,'Mode');

		return ee.Image(classification)


	def runRandomForestModel(self,image,data,bands,nModels,classFieldName,classNames,modelType):
		""" Function to perform bagged classification on an image, using either Random
			Forest or Support Vector Machines (SVM). The data is bagged such that each 
			class has the same number of observations, equal to the minimum number of 
			observations of any class. """

		print image.bandNames().getInfo()
		training_bands = bands # ee.List(["distCoast", "rainy_ND_nir_swir2", "drycool_wetness", \
					#			  "elevation", "dryhot_ND_nir_swir1", "rainy_ND_swir1_swir2", "rainy_tcAngleGW", \
					#			  "dryhot_ND_swir1_swir2", "drycool_tcAngleGW", "rainy_R_red_swir1", \
					#			  "dryhot_tcAngleGW", "drycool_ND_swir1_swir2", "rainy_swir1_stdDev", \
					#			  "dryhot_greenness", "drycool_nir", "drycool_R_red_swir1", "dryhot_thermal", \
					#			  "dryhot_ND_red_swir1", "dryhot_blue", "dryhot_R_red_swir1", \
					#			  "rainy_ND_green_swir1_stdDev", "rainy_nir"]);
		
		data = data.select(training_bands.add(self.env.classFieldName))
		# Find the average number of observations (avgMin = 5000/nClasses)
		nTrees = ee.Number(nModels);
		
		classifier = ee.Classifier.randomForest( nTrees,0).train(data,classFieldName,training_bands);
		
		image = image.select(training_bands)
		classification = image.classify(classifier,'Mode');
  
		return classification


	def newCollectionToImage(self,collection):
		""" Helper function to convert image collection into stack of image bands"""
		
		def createStack(img,prev):
			return ee.Image(prev).addBands(img)
		
		stack = ee.Image(collection.iterate(createStack, ee.Image(1)));

		stack = stack.select(ee.List.sequence(1, stack.bandNames().size().subtract(1)));
		return stack;

       
	def removeDuplicates(self,covariateList,bands):
		""" function to remove duplicates, i.e. existing bands do not need to be calculated """
		
		return [elem for elem in covariateList if elem not in bands]

	def ScaleBands(self,img):
		"""Landsat is scaled by factor 0.0001 """
		
		thermalBand = ee.List(['thermal','thermal_stdDev'])
		gapfillBand = ee.List(['gapfill'])
		
		thermal = ee.Image(img).select(thermalBand).divide(10)
		gapfill = ee.Image(img).select(gapfillBand)
		
		otherBands = ee.Image(img).bandNames().removeAll(thermalBand)
		otherBands = otherBands.removeAll(gapfillBand)
		scaled = ee.Image(img).select(otherBands).multiply(0.0001)
		
		image = ee.Image(scaled.addBands(thermal).addBands(gapfill))
		
		return ee.Image(image.copyProperties(img))

 
	def renameImageBands(self,image,prefix):
		""" Function to add a prefix to all bands in an image """
		
		bandNames = list(image.bandNames().getInfo());
				
		newNames = []
		for band in bandNames:
			bandName = prefix + "_" + band
			newNames.append(bandName)
		
		image = image.rename(ee.List(newNames));
		return image

	def renameBands(self,bandNames,prefix):
		""" Function to add a prefix to all bands in an image """
					
		newNames = []
		for band in bandNames:
			bandName = prefix + "_" + band
			newNames.append(bandName)
		
		
		return newNames

				
	def ExportToAsset(self,name,img):  
		"""export to asset """
		
		outputName = self.env.userID + self.env.assetName  #+ self.env.timeString
		logging.info('export image to asset: ' + str(outputName)) 
		
		region = ee.Geometry.Polygon(self.env.studyArea.bounds().getInfo()['coordinates'])
		
		startDate = ee.Date.fromYMD(self.env.startYear,1,1)
		#endDate = ee.Date.fromYMD(self.env.endYear,12,31)    

		#image = ee.Image(img).set({'system:time_start':startDate.millis(), \
		
		img = img.multiply(100).int16().clip(self.env.studyArea)
		
		img = img.set({'system:time_start':startDate.millis(), \
					   'modeltype':self.env.modelType, \
					   'composite':self.env.composite, \
					   'nmodels':self.env.nModels, \
					   'trainingTable':self.env.trainingData})
				
		task_ordered = ee.batch.Export.image.toAsset(image=ee.Image(img), description=self.env.assetName, assetId=outputName,region=region['coordinates'], maxPixels=1e13,scale=self.env.pixSize )
        
        # start task
		task_ordered.start() 
		


if __name__ == "__main__":
  
	
    # set argument parsing object
	parser = argparse.ArgumentParser(description="Create primitive composite using Google Earth Engine.")
   
	parser.add_argument('--year','-y', type=str,required=True, \
                        help="Year to perform the ats correction and save to asset format in 'YYYY'")

	parser.add_argument('--user','-u', type=str, default="servir-mekong",choices=['servir-mekong','servirmekong',"ate","biplov","quyen","atesig","adpc","KSA"], \
						help="specify user account to run task")

	parser.add_argument('--primitive','-p', required=True,type=str, default="None",choices=["wetlands","mangroves","imperv","grass","cropplantation","rice","irrigated","barren"], \
						help="specify user account to run task")

	args = parser.parse_args() # get arguments  
  
	# user account to run task on
	userName = args.user
	y = int(args.year)
	#self.env.year = year

	# create a new file in ~/.config/earthengine/credentials with token of user
	addUserCredentials(userName)

	ee.Initialize()
	
	env = environment()
   
	# import the images
	dryhot = ee.Image("projects/servir-mekong/usgs_sr_composites/drycool/SC_drycool" + str(int(y)-1) + "_"  + str(y) + env.composite).clip(env.studyArea)
	drycool = ee.Image("projects/servir-mekong/usgs_sr_composites/dryhot/SC_dryhot" + str(y) + "_"  + str(y) + env.composite) .clip(env.studyArea)
	rainy = ee.Image("projects/servir-mekong/usgs_sr_composites/rainy/SC_rainy" + str(y) + "_"  + str(y) + env.composite).clip(env.studyArea)
	
	#if env.composite == "Median":
		#trainingData = ee.FeatureCollection("14AwPqA8ADQU5kCdPN-J2HdD8qUQlREWH97fajYeZ") # all barren, ..
		#trainingData = ee.FeatureCollection("ft:1M71zwCyPyqEEuq8IBrWewONoBVDWHl22D01siCV_") # wetlands
	    #trainingData = ee.FeatureCollection("ft:1pLatAgQLPdzU__wFLZ74591o-t-RwnHNSAayeket") # mangroves
#	selectedBandsMedian = ["drycool_green","drycool_ND_blue_nir","drycool_ND_green_red","drycool_ND_swir1_swir2","drycool_nir","drycool_R_red_swir1","drycool_tcAngleBW","drycool_tcAngleGW","drycool_tcDistBG","dryhot_blue","dryhot_ND_blue_swir1","dryhot_ND_green_red","dryhot_ND_green_swir2","dryhot_nir","dryhot_R_red_swir1","dryhot_tcAngleGW","dryhot_tcDistBG","dryhot_tcDistGW","elevation","rainy_blue","rainy_ND_blue_green","rainy_ND_blue_nir","rainy_ND_blue_red","rainy_ND_blue_swir2","rainy_ND_green_red","rainy_ND_green_swir2","rainy_nir","rainy_R_red_swir1","rainy_R_swir1_nir","rainy_tcAngleBW","rainy_tcAngleGW","rainy_tcDistBG","slope"]
	
	#if env.composite == "Medoid":
	#    trainingData = ee.FeatureCollection("ft:1QyUHohRqWW7HEqag0Yio6z99_LLFJNfmyF7XXEYq")
#	selectedBandsMedoid = ["drycool_blue","drycool_ND_blue_nir","drycool_ND_blue_swir2","drycool_ND_green_red","drycool_ND_swir1_swir2","drycool_nir","drycool_R_red_swir1","drycool_R_swir1_nir","drycool_tcAngleBW","drycool_tcAngleGW","drycool_tcDistBG","dryhot_green","dryhot_ND_blue_swir1","dryhot_ND_blue_swir2","dryhot_ND_green_red","dryhot_ND_green_swir2","dryhot_nir","dryhot_R_red_swir1","dryhot_tcAngleGW","dryhot_tcDistBG","dryhot_tcDistGW","elevation","rainy_blue","rainy_ND_blue_green","rainy_ND_blue_nir","rainy_ND_blue_red","rainy_ND_blue_swir2","rainy_ND_green_red","rainy_ND_green_swir2","rainy_ND_red_swir2","rainy_nir","rainy_R_red_swir1","rainy_tcAngleBW","rainy_tcAngleGW","rainy_tcDistBG","slope"]
	
	
	#trainingData = ee.FeatureCollection("ft:" + env.trainingData)

	primitives().createPrimitive(drycool,dryhot,rainy,y,args.primitive)#,selectedBandsMedoid )

