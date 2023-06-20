#!/usr/bin/env python
# -*- coding: utf8 -*-

# Import standard modules
import os
import os.path
import math
import numpy as np
from scipy import ndimage
import sys

# Modules for plotting
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Import astronomic modules
import aplpy
import pyregion
import astropy.io.fits as pyfits
from astropy import wcs
from astropy.stats import sigma_clip, sigma_clipped_stats
from photutils import CircularAperture
from photutils import aperture_photometry



# -----------------------------------------------------------------

# Disable astropy logging except for warnings and errors
from astropy import log
log.setLevel("WARNING")

# -----------------------------------------------------------------

# Do not show warnings, to block Canopy's UserWarnings from spoiling the console log
import warnings
warnings.filterwarnings("ignore")

# -----------------------------------------------------------------

#*** Colour fonts ***
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

def zoom_fits(fitsfile, scalefactor, preserve_bad_pixels=True, set_wcs=True, **kwargs):
    """
    Zoom in on a FITS image by interpolating using `~scipy.ndimage.interpolation.zoom`

    Parameters
    ----------
    fitsfile : str
        FITS file name
    scalefactor : float
        Zoom factor along all axes
    preserve_bad_pixels : bool
        Try to set NAN pixels to NAN in the zoomed image.  Otherwise, bad
        pixels will be set to zero
    """

    arr = pyfits.getdata(fitsfile)
    h = pyfits.getheader(fitsfile)
    if set_wcs==True:
      try:
	h['CRPIX1'] = (h['CRPIX1']-1)*scalefactor + scalefactor/2. + 0.5
	h['CRPIX2'] = (h['CRPIX2']-1)*scalefactor + scalefactor/2. + 0.5
	if 'CD1_1' in h:
	    for ii in (1, 2):
		for jj in (1, 2):
		    k = "CD%i_%i" % (ii, jj)
		    if k in h:  # allow for CD1_1 but not CD1_2
			h[k] = h[k]/scalefactor
	elif 'CDELT1' in h:
	    h['CDELT1'] = h['CDELT1']/scalefactor
	    h['CDELT2'] = h['CDELT2']/scalefactor
      except:
	print 'ERROR: WCS was not added!'
    bad_pixels = np.isnan(arr) + np.isinf(arr)

    arr[bad_pixels] = 0

    upscaled = ndimage.zoom(arr, scalefactor, **kwargs)

    if preserve_bad_pixels:
        bp_up = ndimage.zoom(bad_pixels, scalefactor,
                                   mode='constant', cval=np.nan, order=0)
        upscaled[bp_up] = np.nan

    up_hdu = pyfits.PrimaryHDU(data=upscaled, header=h)

    return up_hdu


## This function downsamples the frame by a specified zooming factor
def downsample(input_image, factor, output_image, set_wcs=True):
	print bcolors.OKBLUE+ 'Rebining the input image by a factor of:' + bcolors.ENDC, str(factor) 
        # Open the input file and read it
	hdulist_in = pyfits.open(input_image)
	data_in = hdulist_in[0].data
	header_in = hdulist_in[0].header
	
	# Use zoom_fits function from FITS_tools:  http://fits-tools.readthedocs.org/en/latest/_modules/FITS_tools/hcongrid.html
        hdu_out = zoom_fits(input_image, 1.0/factor, set_wcs=set_wcs)
        
        # Save the HDU
	hdu_out.writeto(output_image,clobber=True)
	print 'Done!'
	#return output_image
	



if __name__ == '__main__':
	# Factor, e.g., 2 equals zooming 1/2
	input_image = str(sys.argv[1])
	factor = float(sys.argv[2])
	if len(sys.argv)>3:
	  output_image = str(sys.argv[3])
	else:
	  output_image = input_image.split('.fits')[0] + '_rebin.fits'
	downsample(input_image, factor, output_image)
