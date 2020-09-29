import ee, math

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


	def getIndices(self,img,covariates):	
		""" add indices to image"""
		
		# no need to add indices that are already there
		indices = self.removeDuplicates(covariates,img.bandNames().getInfo())
		
		for item in indices:
			img = self.functionList[item](img)

		return img


	def removeDuplicates(self,covariateList,bands):
		""" function to remove duplicates, i.e. existing bands do not need to be calculated """
		
		return [elem for elem in covariateList if elem not in bands]

	def renameBands(self,image,prefix):
		'rename bands with prefix'
		
		bandnames = image.bandNames();

		def mapBands(band):
			band = ee.String(prefix).cat('_').cat(band);
			return band;
				
		bandnames = bandnames.map(mapBands)
		
		image = image.rename(bandnames);

		return image;

	def addModis(self,img,year):
		
		start = ee.Date.fromYMD(year,1,1)
		end = ee.Date.fromYMD(year+1,1,1)
	
		auto = ee.Image(ee.ImageCollection("users/servirmekong/autocor").filterDate(start,end).first()).rename(["auto"])
		img = img.addBands(auto)
		cycle = ee.Image(ee.ImageCollection("users/servirmekong/seasons").filterDate(start,end).first())
		img = img.addBands(cycle)
		return img

	def addDistCoast(self,img):
		distCoast = ee.Image('projects/servir-mekong/Primitives/DistancetoCoast_1k').float().rename(['distCoast']);
		img = img.addBands(distCoast)
		return img


	def addForest(self,img,year):
		start = ee.Date.fromYMD(year, 1, 1)
		end  = ee.Date.fromYMD(year, 12,31)
		
		tcc = ee.Image(ee.ImageCollection("projects/servir-mekong/Primitives/P_canopy").filterDate(start,end).first()).rename(["tcc"])
		img =img.addBands(tcc)
		treeheight = ee.Image(ee.ImageCollection("projects/servir-mekong/Primitives/P_tree_height").filterDate(start,end).first()).rename(["treeheight"])
		img = img.addBands(treeheight)
		return img
		
	def addWater(self,img,y): 
		
		geom = img.geometry()
		jrc = ee.ImageCollection("JRC/GSW1_0/MonthlyHistory")
		start = ee.Date.fromYMD(y, 1, 1)
		end  = ee.Date.fromYMD(y, 12,31)

		if y > 2015:
			start = ee.Date.fromYMD(2014,1,1)
			end = ee.Date.fromYMD(2016,1,1)
				
		jrc = jrc.filterBounds(geom).filterDate(start, end)
		
		
		def getObs(img):
			obs = img.gt(0)
			return img.addBands(obs.rename(['obs']).set('system:time_start', img.get('system:time_start')));  
		
		def getWater(img):
			water = img.eq(2);
			return img.addBands(water.rename(['onlywater']).set('system:time_start', img.get('system:time_start')));  

		totalObs = jrc.map(getObs).select(["obs"]).sum().toFloat()
		totalWater = jrc.map(getWater).select(["onlywater"]).sum().toFloat()
		
		returnTime = totalWater.divide(totalObs).multiply(100).unmask(0)
		
		return img.addBands(ee.Image(returnTime).rename(["water"]))

	def addOther(self,img):
		protected = ee.Image("projects/servir-mekong/staticMaps/protectedArea").rename(["protected"])
		distRoad = ee.Image("projects/servir-mekong/staticMaps/distRoads").rename(["distRoad"])
		distBuildings = ee.Image("projects/servir-mekong/staticMaps/distBuildings").rename(["distBuildings"])
		distStream = ee.Image("projects/servir-mekong/staticMaps/distStream").rename(["distStream"])
		eco = ee.Image("projects/servir-mekong/staticMaps/ecoRegions").rename(["eco"])
		ecoforest = ee.Image("projects/servir-mekong/staticMaps/ecoRegionsForest").rename(["ecoForest"])
		hand = ee.Image("users/gena/GlobalHAND/90m-global/fa").rename(["hand"])
		img = img.addBands(protected).addBands(distRoad).addBands(distBuildings).addBands(distStream).addBands(eco).addBands(ecoforest).addBands(hand)
		return img
		
def returnCovariates(img,year):
	
	
	# hard coded for now
	bands = ['blue','green','red','nir','swir1', 'swir2']	
	bandLow = ['p20_blue','p20_green','p20_red','p20_nir','p20_swir1', 'p20_swir2']
	bandHigh = ['p80_blue','p80_green','p80_red','p80_nir','p80_swir1', 'p80_swir2']

	"""Calculate the urban, builtup cropland rice and barren primitives """
	covariates = ["ND_blue_green","ND_blue_red","ND_blue_nir","ND_blue_swir1","ND_blue_swir2", \
				  "ND_green_red","ND_green_nir","ND_green_swir1","ND_green_swir2","ND_red_swir1",\
				  "ND_red_swir2","ND_nir_red","ND_nir_swir1","ND_nir_swir2","ND_swir1_swir2",\
				  "R_swir1_nir","R_red_swir1","EVI","SAVI","IBI"]
		
	index = indices()
		
	def addIndices(img,prefix):
		#image = scaleBands(composite)
		img = index.addAllTasselCapIndices(img)
		img = index.getIndices(img,covariates)	
		if len(prefix) > 0:	
			img = index.renameBands(img,prefix)
	
		else:
			year = 2017
			img = index.addJRC(img).unmask(0)
			img = index.addTopography(img).unmask(0)
			img = index.addModis(img,year).unmask(0)
			img = index.addDistCoast(img).unmask(0)
			img = index.addForest(img,year).unmask(0)
			img = index.addWater(img,year).unmask(0)
			img = index.addOther(img).unmask(0)
		return img
			
		
	down = addIndices(img.select(bandLow,bands),"p20")
	middle = addIndices(img.select(bands),"")
	up = addIndices(img.select(bandHigh,bands),"p80")
		
	img = down.addBands(middle).addBands(up)
	return img




	
if __name__ == "__main__":
	ee.Initialize()
	
	primi = "shrub"

	mekongBuffer = ee.FeatureCollection('ft:1LEGeqwlBCAlN61ie5ol24NdUDqB1MgpFR_sJNWQJ');
	mekongRegion = mekongBuffer.geometry();
	aoi = mekongRegion; 

	if primi == "urban":
		data = "ft:1DdRpHL-l-le1K1lQvN3xqm1LP8wPqYF46OqzsoyU"
		bandNames = ee.List(["ND_blue_nir","ND_blue_swir1","ND_swir1_swir2","distRoad","fifth","fourth","nir","p20_ND_blue_red","p20_ND_blue_swir1","p20_ND_green_swir1","p20_ND_swir1_swir2","p20_brightness","p20_fifth","p20_fourth","p20_nir","p20_swir1","p20_tcDistBG","p20_tcDistBW","p20_tcDistGW","p80_ND_blue_green","p80_ND_blue_nir","p80_ND_blue_swir1","p80_ND_nir_swir1","p80_ND_nir_swir2","p80_ND_swir1_swir2","p80_blue","p80_fifth","p80_fourth","p80_nir","p80_swir2","p80_tcAngleGW","p80_tcDistGW","swir1","tcDistBG","tcDistGW","wetness"])	
	if primi == "cropland":
		data = "ft:1f86pTtSqPz1ovbicN-4MvUocAc3TL0SMSdw5xy3g"
		bandNames = ee.List(["ND_green_swir2","ND_red_swir2","R2_cycle3","auto","brightness","distBuildings","distCoast","elevation","fifth","p20_ND_blue_swir2","p20_ND_green_swir1","p20_ND_green_swir2","p20_ND_red_swir1","p20_ND_red_swir2","p20_R_red_swir1","p20_brightness","p20_green","p20_swir1","p20_swir2","p20_tcAngleBW","p20_tcDistBG","p20_tcDistBW","p20_tcDistGW","p20_wetness","p80_ND_red_swir1","p80_red","slope","swir1","swir2","tcAngleBW","tcDistBW","tcc","treeheight"])
	if primi == "rice":
		data = "ft:1H9MPI9ICM-sBwgp2w2wFtAAGPMkt9zrO7_Rci6ik"
		bandNames = ee.List(["R2_cycle1","R2_cycle3","auto","distBuildings","distStream","elevation","green","p20_red","p80_ND_blue_swir2","p80_green","treeheight"])
	if primi == "water":
		data = "ft:1FzxTC6L9AS_qw9lzTBIJJfNCa4s3zRAG2QSsL52A"
		bandNames = ee.List(["ND_blue_swir1","ND_blue_swir2","ND_green_nir","ND_green_swir1","ND_green_swir2","ND_red_swir1","ND_red_swir2","R_red_swir1","auto","distStream","elevation","fifth","hand","nir","p20_ND_blue_swir1","p20_ND_blue_swir2","p20_ND_red_swir1","p20_ND_red_swir2","p20_R_red_swir1","p20_brightness","p20_nir","p20_tcAngleBG","p80_ND_blue_swir1","p80_ND_blue_swir2","p80_ND_green_swir1","p80_ND_green_swir2","p80_ND_red_swir1","p80_ND_red_swir2","p80_brightness","p80_swir1","p80_tcAngleBW","p80_tcAngleGW","p80_tcDistBG","swir1","tcAngleGW","tcDistBG","water"])
	if primi == "shrub":
		data = "ft:1qcSADJ6pmUiQ00KPh1WKXOrN2WTNj3P1VWDdSwlC"
		bandNames = ee.List(["ND_blue_nir","ND_blue_swir1","ND_green_nir","ND_green_swir1","ND_green_swir2","ND_red_swir1","ND_red_swir2","R2_cycle1","R_red_swir1","auto","blue","distBuildings","distRoad","eco","ecoForest","elevation","green","p20_ND_blue_nir","p20_ND_green_nir","p20_ND_red_swir2","p20_R_red_swir1","p20_green","p20_tcAngleBW","p20_tcDistGW","p20_wetness","p80_ND_blue_swir1","p80_ND_green_swir1","p80_ND_green_swir2","p80_ND_red_swir1","p80_ND_red_swir2","p80_R_red_swir1","p80_swir1","p80_tcAngleBW","p80_tcAngleGW","p80_wetness","slope","tcDistGW","tcc","treeheight","wetness"])
	if primi == "mangrove":
		data = "ft:1-3qUQWv55p-bRgPVm2eQXc6e70kQ31l8c6Ohi8nP"
		bandNames = ee.List(["EVI","ND_green_nir","ND_nir_red","ND_nir_swir1","ND_nir_swir2","ND_swir1_swir2","R2_cycle2","R2_cycle3","R_swir1_nir","distCoast","eco","elevation","p20_EVI","p20_IBI","p20_ND_blue_nir","p20_ND_nir_red","p20_ND_nir_swir1","p20_ND_nir_swir2","p20_ND_swir1_swir2","p20_R_swir1_nir","p20_SAVI","p20_blue","p20_greenness","p20_tcAngleBG","p20_tcAngleBW","p20_wetness","p80_ND_nir_red","p80_green","p80_tcAngleBG","swir2","tcAngleBG","treeheight"])
	if primi == "barren":
		data = "ft:18P2eT3DEBRJf9amNLvFL0SZWpOQ72j0iM-zEr2T8"
		bandNames = ee.List(["SAVI","blue","brightness","green","p20_blue","p20_brightness","p20_green","p20_red","p20_swir1","p20_swir2","p20_tcDistBG","p20_tcDistBW","p80_brightness","p80_green","p80_red","p80_swir1","p80_swir2","p80_tcAngleGW","p80_tcDistBG","p80_tcDistBW","p80_wetness","red","swir1","swir2","tcAngleGW","tcDistBG","tcDistBW"])
	if primi == "aquaculture":
		data = "ft:1dEFTpCnUDaD7O4VjjjrBKCutPPjsPDX1QlXI3Wdx"
		bandNames = ee.List(["ND_green_swir1","ND_red_swir1","R2_cycle3","R_red_swir1","aspect","brightness","change_abs","change_norm","distCoast","distRoad","eco","ecoForest","elevation","nir","occurrence","p20_ND_blue_swir1","p20_ND_green_swir1","p20_ND_red_swir1","p20_ND_red_swir2","p20_R_red_swir1","p20_nir","p20_swir1","p20_tcDistBG","p20_tcDistBW","p20_tcDistGW","p80_ND_blue_red","p80_brightness","p80_nir","p80_tcDistBG","p80_tcDistBW","tcDistBG","tcDistBW","tcDistGW","transition"])
	if primi == "wetlands":
		data = "ft:1231-TOPjLj3yFp0cw8vl4SI69mySc_Mb-4QJGEhz"
		bandNames = ee.List(["ND_nir_swir2","auto","change_norm","distBuildings","distCoast","distRoad","eco","elevation","fifth","max_extent","occurrence","p20_ND_swir1_swir2","p20_brightness","p20_fifth","p20_red","p20_swir2","p20_tcDistBG","p20_tcDistBW","p80_ND_green_nir","p80_fifth","p80_green","p80_swir2","sixth","slope","swir2","tcAngleBG","tcc","transition","treeheight"])
	if primi == "plantations":
		data = "ft:1GJysfZ8PLtFDhx9xd-vblrQnhLQxH44AsPje8sqi"
		bandNames = ee.List(["ND_blue_nir","ND_blue_swir1","ND_green_nir","ND_green_swir1","ND_red_swir1","ND_red_swir2","R2_cycle2","R2_cycle3","R_red_swir1","blue","elevation","greenness","p20_EVI","p20_ND_blue_nir","p20_ND_blue_swir1","p20_ND_green_nir","p20_ND_nir_red","p20_ND_red_swir1","p20_ND_red_swir2","p20_R_red_swir1","p20_SAVI","p20_greenness","p20_tcAngleBG","p20_tcDistGW","p80_ND_blue_nir","p80_ND_blue_swir1","p80_ND_green_nir","p80_ND_green_swir1","p80_ND_red_swir1","p80_ND_red_swir2","p80_R_red_swir1","p80_SAVI","p80_greenness","slope","tcc","treeheight"])
	if primi == "grass":
		data = "ft:131tdMB6OU3c9z2ZiCD3XITHkskEatzzm-yUZ5Fhk"
		bandNames = ee.List(["ND_blue_swir1","ND_blue_swir2","ND_green_swir1","ND_green_swir2","ND_nir_red","ND_red_swir1","ND_red_swir2","R2_cycle1","R2_cycle2","R2_cycle3","R_red_swir1","auto","distCoast","distRoad","elevation","green","p20_ND_blue_swir1","p20_ND_blue_swir2","p20_ND_red_swir1","p20_ND_red_swir2","p20_R_red_swir1","p80_ND_blue_swir1","p80_ND_blue_swir2","p80_ND_red_swir2","p80_ND_swir1_swir2","p80_R_red_swir1","p80_SAVI","p80_green","p80_red","protected","red","slope","tcc","treeheight"])
	if primi == "floodedForest":
		data = "ft:1w0qgGJCZhTtNcJzmJabBmLSeLT8g0EtTFkbmwvbs"
		bandNames = ee.List(["EVI","ND_swir1_swir2","SAVI","auto","distBuildings","distCoast","distRoad","eco","elevation","fifth","green","p20_ND_blue_red","p20_ND_swir1_swir2","p20_brightness","p20_fifth","p20_swir1","p20_swir2","p20_tcDistBW","p80_EVI","p80_SAVI","p80_red","p80_tcAngleBG","red","slope","swir2","treeheight"])
	if primi == "tidal":
		data = "ft:1zYwQeeuIvWUC2YvPMUSyRTIeVgEupeY8PGp6MhG6"
		bandNames = ee.List(["ND_green_nir","ND_nir_red","ND_red_swir1","R_red_swir1","SAVI","distCoast","elevation","max_extent","nir","occurrence","p20_ND_blue_nir","p20_ND_blue_swir1","p20_ND_green_nir","p20_ND_red_swir1","p20_ND_red_swir2","p20_R_red_swir1","p20_nir","p80_ND_red_swir1","p80_ND_red_swir2","p80_R_red_swir1","p80_fifth","p80_nir","p80_tcAngleGW","seasonality","slope","tcAngleBG","transition","water"])
	if primi == "evergreen":
		data = "ft:16YokjCqm9FMLDN7Qg6YLfpPoaBb3XEKyHkwbVeHx"
		bandNames = ee.List(["ND_green_nir","ND_nir_red","ND_swir1_swir2","R2_cycle1","R2_cycle2","R2_cycle3","SAVI","distRoad","eco","ecoForest","elevation","green","p20_ND_green_nir","p20_ND_nir_red","p20_R_swir1_nir","p20_SAVI","p20_green","p20_swir1","p20_tcAngleBG","p20_tcAngleGW","p80_ND_green_nir","p80_ND_red_swir1","p80_SAVI","p80_red","p80_tcAngleBG","red","sixth","slope","swir2","tcAngleBG","tcc","treeheight"])		
	if primi == "deciduous":
		data = "ft:1e2ovhYnbh4KoFvcKz_VkRCqSs7i6nzmAH_d-sakk"
		bandNames = ee.List(["EVI","ND_green_swir1","ND_green_swir2","ND_red_swir1","ND_swir1_swir2","R2_cycle1","R2_cycle2","R_red_swir1","SAVI","auto","distBuildings","distRoad","eco","ecoForest","elevation","fifth","p20_EVI","p20_SAVI","p20_fifth","p20_swir1","p20_swir2","p20_tcAngleBG","p20_wetness","p80_EVI","p80_ND_red_swir1","p80_ND_swir1_swir2","p80_R_red_swir1","slope","tcAngleBW","tcAngleGW","tcc","treeheight"])		
	if primi == "mixedForest":
		data = "ft:14HDvA4uMLpNo5nCLnQNWIkBFjCSL6RPozlbTTYb_"
		bandNames = ee.List(["EVI","ND_green_nir","ND_green_swir1","ND_green_swir2","ND_swir1_swir2","R_red_swir1","auto","brightness","elevation","fifth","green","p20_EVI","p20_ND_green_nir","p20_brightness","p20_fifth","p20_green","p20_red","p20_tcDistBW","p80_ND_red_swir1","p80_ND_swir1_swir2","p80_R_red_swir1","p80_SAVI","p80_green","p80_wetness","red","slope","swir2","tcc","treeheight","wetness"])
	if primi == "snow":
		data = "ft:1egnJK-fxIhKmERy4jO7dyzUPpiAVTKSQ2_A3e4_V"
		bandNames = ee.List(["ND_blue_swir2","ND_green_swir1","ND_nir_swir1","R_red_swir1","R_swir1_nir","distBuildings","distRoad","ecoForest","elevation","fourth","p20_ND_blue_swir1","p20_ND_blue_swir2","p20_ND_red_swir1","p20_ND_red_swir2","p20_R_red_swir1","p20_blue","p20_fourth","p20_greenness","p20_red","p20_tcAngleBW","p20_tcDistGW","p20_wetness","p80_ND_blue_swir2","p80_ND_green_nir","p80_ND_green_swir1","p80_ND_green_swir2","p80_ND_red_swir2","p80_green","p80_nir","p80_sixth","p80_tcDistGW","p80_wetness","sixth","slope","wetness"])
		
	"""	
	for year in range(2018,2019,1):
		print primi, year
		img = returnCovariates(ee.Image("projects/servir-mekong/yearlyComposites/composite" + str(year)),year)
		print(img.bandNames().getInfo())
		#training_bands = img.bandNames()
		classifier = ee.Classifier.randomForest(100,0).setOutputMode('PROBABILITY').train(ee.FeatureCollection(data),"land_class",bandNames);
		
		classification = img.classify(classifier,'Mode');
    
		startDate = ee.Date.fromYMD(year,1,1)

		classification = ee.Image(classification.multiply(100).int16().set({"system:time_start": startDate.millis(),"method":"random Forest","trees" : 100,"data":data})).clip(aoi)
		assetid = "projects/servir-mekong/yearly_primitives/" + primi + "/" + primi +"_"+ str(year)
		task_ordered = ee.batch.Export.image.toAsset(image=classification, description=primi+str(year), assetId=assetid,region=img.geometry().getInfo()['coordinates'], maxPixels=1e13,scale=30 )
        
		# start task
		task_ordered.start() 

	"""
	
