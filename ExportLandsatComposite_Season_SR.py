"""
ee_ats_correct.py, SERVIR-Mekong (2017-07-30)

Perfoms atmospheric corrections on Landsat collections
______________________________________________________________________


Usage
------

$ python ee_ats_correct.py {options}

{options} include:

--path (-p)      : required
                 : year to create the Landsat composite image for
                 : in format YYYY

--year (-y)      : required
                 : year to create the Landsat composite image for
                 : in format YYYY-MM-DD

--season (-s)    : required
                 : season to create the Landsat composite image for
                 : seasons are set for the Mekong region and will need to be changed for other geographic areas
                 : options are 'drycool', 'dryhot', or 'rainy'


Example Usage
-------------

1) export surface reflectance composite for dryhot season of 2000 to assets:

  $ python ee_ats_correct.py -p ~/eeAtsCorLut/ -y 2000 -s dryhot

2) export surface reflectance composite for dryhot season of 2005 to assets:

  $ python ee_ats_correct.py -p ~/eeAtsCorLut/ -y 2005 -s drycool

3) export surface reflectance composite for rainy season of 2010 to assets:

  $ python ee_ats_correct.py -p ~/eeAtsCorLut/ -y 2010 -s rainy

"""
from __future__ import division, print_function
import ee
import math
import pickle
import argparse
import datetime

ee.Initialize()

class eeAtsCorrection(object):
    def __init__(self,iniDate,endDate,region,bands,luts,startJulian=1,endJulian=365):
        self.iniDate = iniDate
        self.endDate = endDate
        self.startJulian = startJulian
        self.endJulian = endJulian
        self.region = region
        self.bands = bands
        self.luts = luts
        self.feature = 0

    def getCollection(self,metadataCloudCoverMax=100):

        sensorBandDictLandsatTOA = {'L8': [1,2,3,4,5,9,6],
                                    'L7': [0,1,2,3,4,5,7],
                                    'L5': [0,1,2,3,4,5,6],
                                    'L4': [0,1,2,3,4,5,6]}
        self.bandNamesLandsatTOA = ['blue','green','red','nir','swir1','temp','swir2']

        lt4 = ee.ImageCollection('LANDSAT/LT4_L1T_TOA')\
            .filterBounds(self.region).filterDate(self.iniDate,self.endDate)\
            .filter(ee.Filter.calendarRange(self.startJulian,self.endJulian))\
            .filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)\
            .select(sensorBandDictLandsatTOA['L4'],self.bandNamesLandsatTOA)
        lt5 = ee.ImageCollection('LANDSAT/LT5_L1T_TOA')\
            .filterBounds(self.region).filterDate(self.iniDate,self.endDate)\
            .filter(ee.Filter.calendarRange(self.startJulian,self.endJulian))\
            .filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)\
            .select(sensorBandDictLandsatTOA['L5'],self.bandNamesLandsatTOA)
        le7 = ee.ImageCollection('LANDSAT/LE7_L1T_TOA')\
            .filterBounds(self.region).filterDate(self.iniDate,self.endDate)\
            .filter(ee.Filter.calendarRange(self.startJulian,self.endJulian))\
            .filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)\
            .select(sensorBandDictLandsatTOA['L7'],self.bandNamesLandsatTOA)
        lc8 = ee.ImageCollection('LANDSAT/LC8_L1T_TOA')\
            .filterBounds(self.region).filterDate(self.iniDate,self.endDate)\
            .filter(ee.Filter.calendarRange(self.startJulian,self.endJulian))\
            .filterMetadata('CLOUD_COVER','less_than',metadataCloudCoverMax)\
            .select(sensorBandDictLandsatTOA['L8'],self.bandNamesLandsatTOA)

        collection = ee.ImageCollection(lt4.merge(lt5).merge(le7).merge(lc8))
        self.collectionMeta = collection.getInfo()['features']
        cmasked_collection = collection.map(self.maskClouds)
        masked_collection = self.maskShadows(cmasked_collection)
        return masked_collection


    def maskClouds(self,img,cloudThresh=20):
        def rescale(img,exp,thresholds):
            outimg = img.expression(exp, {img: img})\
                    .subtract(thresholds[0]).divide(thresholds[1] - thresholds[0])
            return outimg

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
        temp_rescale = img.select('temp').subtract(ee.Number(300)).divide(ee.Number(290).subtract(ee.Number(300)))
        score = score.min(temp_rescale);

        # However, clouds are not snow.
        ndsi = img.normalizedDifference(['green', 'swir1']);
        ndsi_rescale = ndsi.subtract(ee.Number(0.8)).divide(ee.Number(0.6).subtract(ee.Number(0.8)))
        score =  score.min(ndsi_rescale).multiply(100).byte();
        score = score.lt(cloudThresh).rename(['cloudMask']);
        img = img.updateMask(score);
        return img.addBands(score);

    def maskShadows(self,collection,zScoreThresh=-0.8,shadowSumThresh=0.35,dilatePixels=2):

        def TDOM(image):
            zScore = image.select(shadowSumBands).subtract(irMean).divide(irStdDev)
            irSum = image.select(shadowSumBands).reduce(ee.Reducer.sum())
            TDOMMask = zScore.lt(zScoreThresh).reduce(ee.Reducer.sum()).eq(2)\
                .And(irSum.lt(shadowSumThresh)).Not()
            TDOMMask = TDOMMask.focal_min(dilatePixels)
            return image.addBands(TDOMMask.rename(['TDOMMask']))

        def mask(image):
            outimg = image.updateMask(image.select(['TDOMMask']))
            return outimg

        shadowSumBands = ['nir','swir1']

        # Get some pixel-wise stats for the time series
        irStdDev = collection.select(shadowSumBands).reduce(ee.Reducer.stdDev())
        irMean = collection.select(shadowSumBands).reduce(ee.Reducer.mean())

        # Mask out dark dark outliers
        collection_tdom = collection.map(TDOM)

        return collection_tdom.map(mask)

    def addIndices(img):
        # Add Normalized Difference Spectral Vector (NDSV)
        img = img.addBands(img.normalizedDifference(['blue','green']).rename('ND_blue_green'))
        img = img.addBands(img.normalizedDifference(['blue','red']).rename('ND_blue_red'))
        img = img.addBands(img.normalizedDifference(['blue','nir']).rename('ND_blue_nir'))
        img = img.addBands(img.normalizedDifference(['blue','swir1']).rename('ND_blue_swir1'))
        img = img.addBands(img.normalizedDifference(['blue','swir2']).rename('ND_blue_swir2'))

        img = img.addBands(img.normalizedDifference(['green','red']).rename('ND_green_red'))
        img = img.addBands(img.normalizedDifference(['green','nir']).rename('ND_green_nir')) # NDWBI
        img = img.addBands(img.normalizedDifference(['green','swir1']).rename('ND_green_swir1')) # NDSI, MNDWI
        img = img.addBands(img.normalizedDifference(['green','swir2']).rename('ND_green_swir2'))

        img = img.addBands(img.normalizedDifference(['red','swir1']).rename('ND_red_swir1'))
        img = img.addBands(img.normalizedDifference(['red','swir2']).rename('ND_red_swir2'))

        img = img.addBands(img.normalizedDifference(['nir','red']).rename('ND_nir_red')) # NDVI
        img = img.addBands(img.normalizedDifference(['nir','swir1']).rename('ND_nir_swir1')) # NDWI, LSWI, -NDBI
        img = img.addBands(img.normalizedDifference(['nir','swir2']).rename('ND_nir_swir2')) # NBR, MNDVI

        img = img.addBands(img.normalizedDifference(['swir1','swir2']).rename('ND_swir1_swir2'))

        # Add ratios
        img = img.addBands(img.select('swir1').divide(img.select('nir')).rename('R_swir1_nir')) # ratio 5/4
        img = img.addBands(img.select('red').divide(img.select('swir1')).rename('R_red_swir1')) # ratio 3/5

        # Add Enhanced Vegetation Index (EVI)
        evi = img.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
          'NIR': img.select('nir'),
          'RED': img.select('red'),
          'BLUE': img.select('blue')
        }).float()
        img = img.addBands(evi.rename('EVI'))

        # Add Soil Adjust Vegetation Index (SAVI)
        # using L = 0.5;
        savi = img.expression(
        '(NIR - RED) * (1 + 0.5)/(NIR + RED + 0.5)', {
          'NIR': img.select('nir'),
          'RED': img.select('red')
        }).float()
        img = img.addBands(savi.rename('SAVI'))

        # Add Index-Based Built-Up Index (IBI)
        ibi_a = img.expression(
        '2*SWIR1/(SWIR1 + NIR)', {
          'SWIR1': img.select('swir1'),
          'NIR': img.select('nir')
        }).rename('IBI_A')
        ibi_b = img.expression(
        '(NIR/(NIR + RED)) + (GREEN/(GREEN + SWIR1))', {
          'NIR': img.select('nir'),
          'RED': img.select('red'),
          'GREEN': img.select('green'),
          'SWIR1': img.select('swir1')
        }).rename('IBI_B')
        ibi_a = ibi_a.addBands(ibi_b)
        ibi = ibi_a.normalizedDifference(['IBI_A','IBI_B'])
        img = img.addBands(ibi.rename('IBI'))

        return img

    def getTasseledCap(image,bands) :
        # Kauth-Thomas coefficients for Thematic Mapper data
        coefficients = ee.Array([
                [0.3037, 0.2793, 0.4743, 0.5585, 0.5082, 0.1863],
                [-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800],
                [0.1509, 0.1973, 0.3279, 0.3406, -0.7112, -0.4572],
                [-0.8242, 0.0849, 0.4392, -0.0580, 0.2012, -0.2768],
                [-0.3280, 0.0549, 0.1075, 0.1855, -0.4357, 0.8085],
                [0.1084, -0.9022, 0.4120, 0.0573, -0.0251, 0.0238]
                ]);
        # Make an Array Image, with a 1-D Array per pixel.
        arrayImage1D = image.select(bands).toArray();

        # Make an Array Image with a 2-D Array per pixel, 6x1.
        arrayImage2D = arrayImage1D.toArray(1);

        componentsImage = ee.Image(coefficients)\
        .matrixMultiply(arrayImage2D)\
        .arrayProject([0])\
        .arrayFlatten( # Get a multi-band image with TC-named bands.
          [['brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth']])\
        .float();

        return image.addBands(componentsImage)

    def correctAts(self,image):
        esun = {'TM': [1958.00,1827.00,1551.00,1036.00,214.90,80.65],
                'ETM':[1969.00,1840.00,1551.00,1044.00,225.70,82.07],
                'OLI':[2004.57,1820.75,1549.49, 951.76,247.55,85.46]}

        exportBands = ['blue','green','red','nir','swir1','swir2']

        metaData = self.collectionMeta[self.feature]['properties']

        sensorId = str(metaData['SENSOR_ID'])
        if len(sensorId) > 3:
            sensorId = sensorId.split('_')[0]

        sunAngle = 90-float(metaData['SUN_ELEVATION'])

        llLat = float(metaData['CORNER_LL_LAT_PRODUCT'])
        llLon = float(metaData['CORNER_LL_LON_PRODUCT'])
        ulLat = float(metaData['CORNER_UL_LAT_PRODUCT'])
        ulLon = float(metaData['CORNER_UL_LON_PRODUCT'])
        lrLat = float(metaData['CORNER_LR_LAT_PRODUCT'])
        lrLon = float(metaData['CORNER_LR_LON_PRODUCT'])
        urLat = float(metaData['CORNER_UR_LAT_PRODUCT'])
        urLon = float(metaData['CORNER_UR_LON_PRODUCT'])

        imgBounds = ee.Geometry.Polygon([[llLon,llLat],[ulLon,ulLat],[urLon,urLat],[lrLon,lrLat],[llLon,llLat]])

        srtm = ee.Image('USGS/SRTMGL1_003')
        alt = srtm.reduceRegion(reducer=ee.Reducer.mean(),geometry=imgBounds,bestEffort=True).get('elevation').getInfo()/1000.

        acqDate = metaData['DATE_ACQUIRED'].split('-')
        dateob = datetime.date(int(acqDate[0]),int(acqDate[1]),int(acqDate[2]))
        eeDate = ee.Date.fromYMD(int(acqDate[0]),int(acqDate[1]),int(acqDate[2]))
        doy = dateob.timetuple().tm_yday
        orbit_correction = 0.03275104*math.cos(doy/59.66638337) + 0.96804905

        dsun = (1 - 0.01672 * math.cos(0.9856 * (doy - 4))) ** 2

        aot_ic =  ee.ImageCollection('MODIS/MOD08_M3_051')\
                        .filterDate(ee.Date.fromYMD(int(acqDate[0]),int(acqDate[1]),1),
                                    ee.Date.fromYMD(int(acqDate[0]),int(acqDate[1]),1)
                                        .advance(1,'month'))
        aot_img = ee.Image(aot_ic.first())
        try:
            aot = aot_img.reduceRegion(reducer=ee.Reducer.mean(), geometry=imgBounds).get('AOT_550').getInfo()
        except:
            aot = 0.1

        water_ic = ee.ImageCollection('NCEP_RE/surface_wv')\
            .filterDate(eeDate,eeDate.advance(1,'month'))
        water_img = ee.Image(water_ic.first())
        water = water_img.reduceRegion(reducer=ee.Reducer.mean(), geometry=imgBounds).get('pr_wtr')
        try:
            h2o = ee.Number(water).divide(10).getInfo()
        except:
            h2o = 0.25

        ozone_ic = ee.ImageCollection('users/samsammurphy/public/ozone_fill').toList(366)
        ozone_img = ee.Image(ozone_ic.get(doy))
        o3DU = ozone_img.reduceRegion(reducer=ee.Reducer.mean(), geometry=imgBounds).get('ozone')
        o3 = ee.Number(o3DU).divide(1000).getInfo()

        outsr = ee.Image(0)

        for i in range(len(self.bands)):
            a,b = self.luts[sensorId][self.bands[i]](sunAngle,h2o,o3,aot,alt)
            band = image.select(self.bandNamesLandsatTOA[i])
            denom = esun[sensorId][i]*math.cos(math.radians(sunAngle))
            numer = math.pi * dsun
            radiance = band.multiply(ee.Number(denom)).divide(ee.Number(numer))
            a *= orbit_correction
            b *= orbit_correction
            outsr = outsr.addBands((radiance.subtract(ee.Number(a)).divide(ee.Number(b)).rename([self.bands[i]])))

        self.feature += 1

        return outsr

    def correct(self):
        collection = self.getCollection() # make image collection with masks
        atsCorrected = collection.map(self.correctAts) # atmospheric correction
        return atsCorrected

def main():
    # set argument parsing object
    parser = argparse.ArgumentParser(description="Perform atmospheric corrections on\
                                                  Landsat image collections using Google\
                                                  Earth Engine and the 6S Emulator code.")

    # specify arguments for script and help informations
    parser.add_argument('--path','-p', type=str,required=True,
                        help="Path to lookup tables to correct the atmosphere")
    parser.add_argument('--year','-y', type=str,required=True,
                        help="Year to perform the ats correction and save to asset format in 'YYYY'")
    parser.add_argument('--season','-s', choices=['drycool','dryhot','rainy'],type=str,
                        help="Season to create composite for, these align with SERVIR-Mekong's seasonal composite times")
    parser.add_argument('--month','-m', choices=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'],type=str,
                        help="Month to create composite for")

    args = parser.parse_args() # get arguments

    # input path for the atmospheric correction lookup tables
    inpath = args.path

    # set sensor and band variables
    sensors = ['TM','ETM','OLI']
    bands = ['B1','B2','B3','B4','B5','B7']
    exportBands = ['blue','green','red','nir','swir1','temp','swir2']

    # blank dictionary to store lookup tables
    luts = {}

    # loop over all of the sensors and bands to get lookup tables
    for i in range(len(sensors)):
        temp = {} # blank temp dict
        for j in range(len(bands)):
            with open(inpath+'LANDSAT_{0}_{1}.ilut'.format(sensors[i],bands[j]), 'rb') as f:
                lut = pickle.load(f)
            temp['{}'.format(bands[j])] = lut # set band lookup table

        luts['{}'.format(sensors[i])] = temp # set sensor lookup tables for all bands

    # get the Mekong coutries + buffer region
    mekongBuffer = ee.FeatureCollection('ft:1LEGeqwlBCAlN61ie5ol24NdUDqB1MgpFR_sJNWQJ')
    mekongRegion = mekongBuffer.geometry()

    # handle seasons and the date ranges
    seasons = {'drycool':[305,59],'dryhot':[60,120],'rainy':[121,304]}


    months = {'jan':[1,31],'feb':[32,59]}

    if args.season == None:
        timing = args.month
        month = months[args.month]
        startJulian = month[0]
        endJulian = month[1]
    else:
        timing = args.season
        season = seasons[args.season]
        startJulian = season[0]
        endJulian = season[1]

    # read in year argument
    yr = int(args.year)

    # make data strings
    iniDate = '{}-01-01'.format(yr-1)
    endDate = '{}-12-31'.format(yr+1)

    # define parameters to create composite
    metadataCloudCoverMax = 100;
    cloudThresh = 20;
    dilatePixels = 2;
    cloudHeights = ee.List.sequence(200,5000,500);
    zScoreThresh = -0.8;
    shadowSumThresh = 0.35;

    print(startJulian, endJulian)

    # initialize atmospheric correction class and start process
    processor = eeAtsCorrection(iniDate,endDate,mekongRegion,bands,luts,startJulian,endJulian)
    collection_sr = processor.correct()

    # get the median SR value for compposite
    season_sr = collection_sr.select(exportBands).median().clip(mekongRegion)

    # scale image to interger values
    compositeSR = season_sr.multiply(10000).int16()

    # set up export paths
    exportName = 'Landsat_' + timing + '_' + iniDate.split('-')[0] + '_' + endDate.split('-')[0]
    #exportPath = 'projects/servir-mekong/Composites/Landsat_'  + args.season + '/' + exportName;
    exportPath = 'users/kelmarkert/' + exportName

    # define metadata for asset
    snippetName = 'ExportComposite_Seasonal_SR'
    projectName = 'SERVIR-Mekong'
    versionNumber = 1

    # get serializable geometry for export
    exportRegion = mekongRegion.bounds().getInfo()['coordinates']

    # set metadata information
    compositeSR = compositeSR.set({
      'system:time_start': ee.Date.fromYMD(yr,6,1).millis(),
      'date': ee.Date.fromYMD(yr,6,1),
      'source':'SR',
      'snippet': snippetName,
      'project': projectName,
      'collection': exportName,
      'version': versionNumber})

    exportTask = "Landsat_SR_{0}_{1}_{2}".format(args.season,yr-1,yr+1)

    # create export task
    export = ee.batch.Export.image.toAsset(compositeSR,
      description=exportTask,
      assetId= exportPath,
      scale= 30,
      region=exportRegion,
      maxPixels=1e13,
      crs='EPSG:4326'
    )
    export.start() # begin task

    print('\nExport has begun for image {0} to collection {1}\n'.format(exportName,exportPath.split('/')[-1]))

    return

if __name__ == "__main__":
    main()
