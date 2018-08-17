# LandCoverMonitoring

The USGSComposite folder contains scripts to:
1. Create the composite (usgs_composite.py)
2. Create a training data set (trainingData.py)
3. Create the primivites (calculatePrimitives.py)

Usage:
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
