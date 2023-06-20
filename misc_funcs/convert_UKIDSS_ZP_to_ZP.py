#! /usr/bin/env python

import sys
from os.path import exists
from os import remove
from scipy.ndimage.filters import convolve
from numpy import max, sum
from astropy.io import fits as pyfits
from math import *

def main(input_image):
    hdulist = pyfits.open(input_image)
    image = hdulist[1].data
    header0 = hdulist[0].header
    header1 = hdulist[1].header
    header_zero_point = float(header1['MAGZPT'])
    header_exptime = float(header0['EXP_TIME'])
    header_extinction = float(header1['EXTINCT'])
    header_airmass = (float(header0['AMSTART']) + float(header0['AMEND']))/2.
    
    gain = float(header1['GAIN'])
    RON = float(header1['READNOIS'])
    FWHM = float(header1['SEEING'])
    NCOMBINE = int(header0['NEXP'])
    
    true_zero_point = header_zero_point - (2.5*log10(1./header_exptime)) - (header_extinction*(header_airmass-1.0))# - vega_to_AB
    print('Zero point is %.4f' % (true_zero_point))
    return true_zero_point,gain,RON,FWHM,NCOMBINE,header_exptime

if __name__ == "__main__":
    input_image = sys.argv[1]

    main(input_image)

    
