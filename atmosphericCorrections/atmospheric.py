"""
atmospheric.py
Created by Sam Murphy (2016-10-26)
Updated by Kel Markert (2017-9-27)

Atmospheric water vapour, ozone and AOT from GEE

Usage:

atmoData = Atmospheric(geom,eeDate)

aot = atmoData.aerosol().getInfo()
h2o = atmoData.water().getInfo()
o3 = atmoData.ozone().getInfo()

"""


import ee

class Atmospheric(object):

  def __init__(self,geom,date):
      self.geom = geom
      self.date = date

  def round_date(self,date,xhour):
    """
    rounds a date of to the closest 'x' hours
    """
    y = date.get('year')
    m = date.get('month')
    d = date.get('day')
    H = date.get('hour')
    HH = H.divide(xhour).round().multiply(xhour)
    return date.fromYMD(y,m,d).advance(HH,'hour')

  def round_month(self,date):
    """
    round date to closest month
    """
    # start of THIS month
    m1 = date.fromYMD(date.get('year'),date.get('month'),ee.Number(1))

    # start of NEXT month
    m2 = m1.advance(1,'month')

    # difference from date
    d1 = ee.Number(date.difference(m1,'day')).abs()
    d2 = ee.Number(date.difference(m2,'day')).abs()

    # return closest start of month
    return ee.Date(ee.Algorithms.If(d2.gt(d1),m1,m2))



  def water(self):
    """
    Water vapour column above target at time of image aquisition.

    (Kalnay et al., 1996, The NCEP/NCAR 40-Year Reanalysis Project. Bull.
    Amer. Meteor. Soc., 77, 437-471)
    """

    # H2O datetime is in 6 hour intervals
    H2O_date = self.round_date(self.date,6)

    # filtered water collection
    water_ic = ee.ImageCollection('NCEP_RE/surface_wv').filterDate(H2O_date, H2O_date.advance(1,'month'))

    # water image
    water_img = ee.Image(water_ic.first())

    # water_vapour at target
    water = water_img.reduceRegion(reducer=ee.Reducer.mean(), geometry=self.geom).get('pr_wtr')

    # convert to Py6S units (Google = kg/m^2, Py6S = g/cm^2)
    water_Py6S_units = ee.Number(water).divide(10)

    return water_Py6S_units


  def ozone(self):
    """
    returns ozone measurement from merged TOMS/OMI dataset

    OR

    uses our fill value (which is mean value for that latlon and day-of-year)

    """

    def ozone_measurement(geom,O3_date):

      # filtered ozone collection
      ozone_ic = ee.ImageCollection('TOMS/MERGED').filterDate(O3_date, O3_date.advance(1,'month'))

      # ozone image
      ozone_img = ee.Image(ozone_ic.first())

      # ozone value IF TOMS/OMI image exists ELSE use fill value
      ozone = ee.Algorithms.If(ozone_img,\
      ozone_img.reduceRegion(reducer=ee.Reducer.mean(), geometry=geom).get('ozone'),\
      ozone_fill(self.geom,O3_date))

      return ozone

    def ozone_fill(geom,O3_date):
      """
      Gets our ozone fill value (i.e. mean value for that doy and latlon)

      you can see it
      1) compared to LEDAPS: https://code.earthengine.google.com/8e62a5a66e4920e701813e43c0ecb83e
      2) as a video: https://www.youtube.com/watch?v=rgqwvMRVguI&feature=youtu.be

      """

      # ozone fills (i.e. one band per doy)
      ozone_fills = ee.ImageCollection('users/samsammurphy/public/ozone_fill').toList(366)

      # day of year index
      jan01 = ee.Date.fromYMD(O3_date.get('year'),1,1)
      doy_index = self.date.difference(jan01,'day').toInt()# (NB. index is one less than doy, so no need to +1)

      # day of year image
      fill_image = ee.Image(ozone_fills.get(doy_index))

      # return scalar fill value
      return fill_image.reduceRegion(reducer=ee.Reducer.mean(), geometry=geom).get('ozone')

    # O3 datetime in 24 hour intervals
    O3_date = self.round_date(self.date,24)

    # TOMS temporal gap
    TOMS_gap = ee.DateRange('1994-11-01','1996-08-01')

    # avoid TOMS gap entirely
    ozone = ee.Algorithms.If(TOMS_gap.contains(O3_date),ozone_fill(self.geom,O3_date),ozone_measurement(self.geom,O3_date))

    # fix other data gaps (e.g. spatial, missing images, etc..)
    ozone = ee.Algorithms.If(ozone,ozone,ozone_fill(self.geom,O3_date))

    #convert to Py6S units
    ozone_Py6S_units = ee.Number(ozone).divide(1000)# (i.e. Dobson units are milli-atm-cm )

    return ozone_Py6S_units


  def aerosol(self):
    """
    Aerosol Optical Thickness.

    try:
      MODIS Aerosol Product (monthly)
    except:
      fill value
    """

    def aerosol_fill(date):
      """
      MODIS AOT fill value for this month (i.e. no data gaps)
      """
      jan01 = ee.Date.fromYMD(date.get('year'),1,1);
      doy_index = date.difference(jan01,'day');

      newaodFill = ee.ImageCollection('users/kelmarkert/public/MOD04L2_GRIDDED_CLIMATOLOGY')
      aodFill = newaodFill.toList(366)

      fillImage = ee.Image(aodFill.get(doy_index.toUint16()))\
            .select(['aod550']).rename(['AOT_550']).divide(1000)

      return fillImage


    def aerosol_this_day(date):
      """
      MODIS AOT original data product for this month (i.e. some data gaps)
      """
      # image for this month
      img =  ee.Image(\
                      ee.ImageCollection('users/kelmarkert/public/MOD04L2_GRIDDED_DAILY')\
                        .filterDate(date)\
                        .first()\
                     )

      # fill missing month (?)
      img = ee.Algorithms.If(img,\
                               # all good
                               img\
                               .select(['b1'])\
                               .divide( 1000)\
                               .rename(['AOT_550']),\
                              # missing month
                                aerosol_fill(date))

      return img


    def get_AOT(AOT_band,geom):
      """
      AOT scalar value for target
      """
      return ee.Image(AOT_band).reduceRegion(reducer=ee.Reducer.mean(),\
                                 geometry=geom)\
                                .get('AOT_550')


    after_modis_start = self.date.difference(ee.Date('2000-03-01'),'month').gt(0)

    AOT_band = ee.Algorithms.If(after_modis_start, aerosol_this_day(self.date), aerosol_fill(self.date))

    AOT = get_AOT(AOT_band,self.geom)

    AOT = ee.Algorithms.If(AOT,AOT,get_AOT(aerosol_fill(self.date),self.geom))
    # i.e. check reduce region worked (else force fill value)

    return AOT
