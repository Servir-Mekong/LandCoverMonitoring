"""
ExportLandsatSRComposite.py, SERVIR-Mekong (2017-07-30)

Perfoms atmospheric corrections on Landsat collections
________________________________________________________________________________


Usage
------

$ python ee_ats_correct.py {options}

{options} include:

--path (-p)      : required
                 : path to the folder containing the atmospheric corrections LUTs
                 : must end with "/"

--year (-y)      : required
                 : year to create the Landsat composite image for
                 : in format YYYY

--season (-s)    : season to create the Landsat composite image for
                 : seasons are set for the Mekong region and will need to be \
                 : changed for other geographic areas
                 : options are 'drycool', 'dryhot', or 'rainy'

--month (-m)     : month to create the Landsat composite image for
                 : month input strings should be in the 3-letter abbreviations \
                 : i.e. January = jan, February = feb, and so on
                 : options are 'jan', 'feb', 'mar', 'apr', 'may', 'jun',
                 :             'jul', 'aug', 'sep', 'oct', 'nov', 'dec'


Example Usage
-------------

1) export surface reflectance composite for dryhot season of 2000 to assets:

  $ python ExportLandsatSRComposite.py -p ~/eeAtsCorLut/ -y 2000 -s dryhot

2) export surface reflectance composite for the month of july 2005 to assets:

  $ python ExportLandsatSRComposite.py -p ./eeAtsCorLut/ -y 2005 -m jul

3) export surface reflectance composite for rainy season of 2010 to assets:

  $ python ExportLandsatSRComposite.py -p ~/eeAtsCorLut/ -y 2010 -s rainy

"""
from __future__ import division, print_function
import ee
import sys
import math
import pickle
import argparse
import datetime
import numpy as np
from atmospheric import Atmospheric

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
        self.cloudThresh = 10

    def getCollection(self,metadataCloudCoverMax=90):

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
            .select(sensorBandDictLandsatTOA['L5'],self.bandNamesLandsatTOA)\
            .map(self.defringeLandsat)
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
        cmasked_collection = collection.map(self.maskClouds)
        masked_collection = self.maskShadows(cmasked_collection)
        self.collectionMeta = masked_collection.getInfo()['features']
        return masked_collection


    def maskClouds(self,img,cloudThresh=10):
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
        mask = score.lt(self.cloudThresh).rename(['cloudMask']);
        img = img.updateMask(mask);
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

    def correctAts(self,image):

        keepBands = ['temp']

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

        elv = ee.Image('CGIAR/SRTM90_V4')
        srtm = ee.Image(0).where(elv.gt(0),elv)
        #minAlt = np.abs(srtm.reduceRegion(reducer=ee.Reducer.min(),scale=90,geometry=imgBounds,maxPixels=9e9).get('elevation').getInfo())
        maxAlt = np.abs(srtm.reduceRegion(reducer=ee.Reducer.max(),scale=90,geometry=imgBounds,maxPixels=9e9).get('constant').getInfo())

        acqDate = metaData['DATE_ACQUIRED'].split('-')
        dateob = datetime.date(int(acqDate[0]),int(acqDate[1]),int(acqDate[2]))
        eeDate = ee.Date.fromYMD(int(acqDate[0]),int(acqDate[1]),int(acqDate[2]))
        doy = dateob.timetuple().tm_yday
        orbit_correction = 0.03275104*math.cos(doy/59.66638337) + 0.96804905

        dsun = (1 - 0.01672 * math.cos(0.9856 * (doy - 4))) ** 2

        atmoData = Atmospheric(imgBounds,eeDate)

        aot = atmoData.aerosol().getInfo()
        h2o = atmoData.water().getInfo()
        o3 = atmoData.ozone().getInfo()

        imgList = []

        altRange = np.arange(0,maxAlt+60,15)

        for i in range(altRange.size):
            elvMask = srtm.gte(altRange[i]-7.5).And(srtm.lt(altRange[i]+7.5))
            temp = image.updateMask(elvMask)
            outsr = image.select(keepBands)

            for j in range(len(self.bands)):
                a,b = self.luts[sensorId][self.bands[j]](sunAngle,h2o,o3,aot,(altRange[i]/1000.))
                band = temp.select(self.bands[j])
                denom = esun[sensorId][j]*math.cos(math.radians(sunAngle))
                numer = math.pi * dsun
                radiance = band.multiply(ee.Number(denom)).divide(ee.Number(numer))
                a *= orbit_correction
                b *= orbit_correction
                outsr = outsr.addBands((radiance.subtract(ee.Number(a)).divide(ee.Number(b))\
                             .rename([self.bands[j]])))
            imgList.append(outsr.where(outsr.lt(0),0))

        outImg = ee.ImageCollection(imgList).median().copyProperties(image,['system:time_start'])

        self.feature += 1

        return outImg

    def defringeLandsat(self,img):
        k = ee.Kernel.fixed(41, 41,
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
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
            ,-1,-1,True);
        # Find any pixel without data
        m = img.reduce(ee.Reducer.min());

        # Find pixels without data in all pixels in the kernel
        sum41 = m.reduceNeighborhood(ee.Reducer.sum(), k, 'kernel');
        sum41 = sum41.gte(1);
        img = img.mask(sum41);
        return img;


    def correct(self):
        collection = self.getCollection() # make image collection with masks
        atsCorrected = collection.map(self.correctAts) # atmospheric correction
        return atsCorrected

def metoidMosiac(inCollection):
    def calcDistance(img):
        diff = ee.Image(img).subtract(median).pow(ee.Image.constant(2));
        return diff.reduce('sum').addBands(img);

    f = ee.Image(inCollection.first())
    bandNames = f.bandNames()
    bandNumbers = ee.List.sequence(1,bandNames.length())

    median = ee.ImageCollection(inCollection).median()

    medoid = inCollection.map(calcDistance)

    medoid = ee.ImageCollection(medoid).reduce(ee.Reducer.min(bandNames.length().add(1))).select(bandNumbers,bandNames);

    return medoid

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
    parser.add_argument('--box','-b',nargs='*',required=True,default=[103,12,104,13],
                        help="Bounding box to create export for in format 'xMin yMin xMax yMax'")
    args = parser.parse_args() # get arguments

    # input path for the atmospheric correction lookup tables
    inpath = args.path

    # check correct input for bounding box
    if len(args.box) != 4:
        print('ERROR: Bounding box argument requires 4 inputs. Exiting...')
        sys.exit(1)
    else:
        for i in range(len(args.box)):
            args.box[i] = int(args.box[i])

    # set sensor and band variables
    sensors = ['TM','ETM','OLI']
    #bands = ['B1','B2','B3','B4','B5','B7']
    bands = ['blue','green','red','nir','swir1','swir2']

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
    #mekongBuffer = ee.FeatureCollection('projects/servir-hkh/LCMS/Mizoram_Pilot/mizoram_area')
    mekongRegion = ee.Geometry.Rectangle([args.box[0],args.box[1],args.box[2],args.box[3]])

    # handle seasons and the date ranges
    seasons = {'drycool':[1,1,305,59],'dryhot':[4,1,60,120],'rainy':[8,1,121,304]}

    # handle months and the date ranges
    months = {'jan':[1,15,1,31],'feb':[2,15,32,59],'mar':[3,15,60,90],'apr':[4,15,91,120],
              'may':[5,15,121,151],'jun':[6,15,152,181],'jul':[7,15,182,212],'aug':[8,15,213,243],
              'sep':[9,15,244,273],'oct':[10,15,274,304],'nov':[11,15,305,334],'dec':[12,15,335,365]}

    # read in year argument
    yr = int(args.year)
    sYr = yr-1
    eYr = yr+1

    if args.month:
        timing = args.month
        month = months[args.month]
        startJulian = month[2]
        endJulian = month[3]
    elif args.season:
        timing = args.season
        month = seasons[args.season]
        startJulian = month[2]
        endJulian = month[3]
    else:
        print('WARNING: No seasonal or monthly information provided\n'\
              'creating composite for entire year')
        timing = 'yearly'
        month = [6,15,1,1]
        startJulian = 1
        endJulian = 365
        sYr = yr
        eYr = yr

    # make data strings
    iniDate = '{}-01-01'.format(sYr)
    endDate = '{}-12-31'.format(eYr)

    # define parameters to create composite
    metadataCloudCoverMax = 100;
    cloudThresh = 20;
    dilatePixels = 2;
    cloudHeights = ee.List.sequence(200,5000,500);
    zScoreThresh = -0.8;
    shadowSumThresh = 0.35;

    # define list of bands to composite
    # medianBands = ['blue','green','red','nir','swir1','temp','swir2',
    #     'ND_blue_green','ND_blue_red','ND_blue_nir','ND_blue_swir1','ND_blue_swir2',
    #     'ND_green_red','ND_green_nir','ND_green_swir1','ND_green_swir2','ND_red_swir1',
    #     'ND_red_swir2','ND_nir_red','ND_nir_swir1','ND_nir_swir2','ND_swir1_swir2',
    #     'R_swir1_nir','R_red_swir1','EVI','SAVI','IBI','brightness','greenness',
    #     'wetness','fourth', 'fifth', 'sixth','tcAngleBG','tcDistBG','tcAngleGW',
    #     'tcDistGW','tcAngleBW','tcDistBW']
    medianBands = ['blue','green','red','nir','swir1','temp','swir2']

    stdDevBands = ['blue','green','red','nir','swir1','temp','swir2']

    # initialize atmospheric correction class and start process
    processor = eeAtsCorrection(iniDate,endDate,mekongRegion,bands,luts,startJulian,endJulian)
    collection_sr = processor.correct()

    # get number of observations for the composite
    countPost = collection_sr.select(['blue']).count().rename(['count']);

    # get the median SR value for compposite
    #srImg = metoidMosiac(collection_sr.select(medianBands)).clip(mekongRegion)
    srImg = metoidMosiac(collection_sr.select(medianBands)).clip(mekongRegion)

    #get the standard deviation SR values
    stdDevComposite = collection_sr.select(stdDevBands).reduce(ee.Reducer.stdDev())

    # add the std dev bands and count band to the export image
    srImg = srImg.addBands(stdDevComposite).addBands(countPost)

    # scale the images based on bands
    compositeBands = sum([medianBands, ['blue_stdDev','green_stdDev','red_stdDev','nir_stdDev','swir1_stdDev','swir2_stdDev','count']], []) #srImg.bandNames().getInfo();
    nonDivideBands = ['temp','temp_stdDev','count','SR_mask','R_swir1_nir','R_red_swir1'];
    composite10000 = srImg.select([b for b in compositeBands if b not in nonDivideBands]).multiply(10000);
    composite10 = srImg.select(['temp'],['thermal']).multiply(10);
    composite1 = srImg.select(['count']);
    compositeSR= composite10000.addBands(composite10).addBands(composite1).int16();

    # set up export paths
    exportName = 'Landsat_{0}_{1}_{2}_{3}_{4}_{5}_{6}'.format(timing,iniDate.split('-')[0],endDate.split('-')[0],args.box[0],args.box[1],args.box[2],args.box[3])
    exportPath = 'projects/servir-mekong/Composites/Landsat_'  + timing + '/' + exportName

    # define metadata for asset
    snippetName = 'ExportComposite_SR'
    projectName = 'SERVIR-MK'
    versionNumber = 1

    # get serializable geometry for export
    exportRegion = mekongRegion.bounds().getInfo()['coordinates']

    # set metadata information
    compositeSR = compositeSR.set({
      'system:time_start': ee.Date.fromYMD(yr,month[0],month[1]).millis(),
      'date': ee.Date.fromYMD(yr,month[0],month[1]),
      'source':'SR',
      'snippet': snippetName,
      'project': projectName,
      'collection': exportName,
      'version': versionNumber})

    exportTask = "Landsat_SR_{0}_{1}_{2}".format(timing,sYr,eYr)

    # create export task
    export = ee.batch.Export.image.toAsset(compositeSR,
      description=exportTask,
      assetId= exportPath,
      scale = 30,
      region=exportRegion,
      maxPixels=1e13,
      crs='EPSG:4326'
    )
    export.start() # begin task

    print('\nExport has begun for image {0} to collection {1}\n'.format(exportName,exportPath.split('/')[-2]))

    return

if __name__ == "__main__":
    main()
