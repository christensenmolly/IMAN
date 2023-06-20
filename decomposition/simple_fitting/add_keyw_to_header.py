#!/usr/bin/python
# DESCRIPTION:
# Add keywords to header.
# MINIMAL USAGE: python add_keyw_to_header.py [input_image] [--KEYWORD] [VALUE]


# Import standard modules
import numpy as np
import argparse
from astropy.io import fits as pyfits

def write_keyw(header, keyword, value):
    value = float(value)
    if np.isnan(value)==False:
        if keyword in header.keys():
          header[keyword] = value
        else:
          header.append((keyword,value),end=True)
    return header

def add_to_header(input_image,EXPTIME,GAIN,NCOMBINE,RON,m0,pix2sec,fwhm,sky_level=0.,sky_subtr=1,xc=None,yc=None):
        """ The function to update the fits header adding some keywords to it.
        Theses are required for Galfit if you do not have a weight image!"""
        hdulist = pyfits.open(input_image, mode='update')
        prihdr = hdulist[0].header

        prihdr = write_keyw(prihdr,'EXPTIME', EXPTIME)
        prihdr = write_keyw(prihdr,'GAIN', GAIN)
        prihdr = write_keyw(prihdr,'NCOMBINE', NCOMBINE)
        prihdr = write_keyw(prihdr,'RDNOISE', RON)
        #prihdr = write_keyw(prihdr,'M0', m0)
        prihdr = write_keyw(prihdr,'SCALE', pix2sec)
        prihdr = write_keyw(prihdr,'FWHM', fwhm)
        prihdr = write_keyw(prihdr,'SKY_LEVEL', sky_level)
        prihdr = write_keyw(prihdr,'SKY_SUBTR', sky_subtr)
        if xc!=None and yc!=None:
          prihdr = write_keyw(prihdr,'XC', xc)
          prihdr = write_keyw(prihdr,'YC', yc)
        hdulist.flush()

'''
if __name__ == '__main__':
    # TODO:
    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Add keywords to header")
    #parser.add_argument("input_image", help="Input fits image")
    
    args = parser.parse_args()
    print(args)
'''
        



'''
input_image = 'galaxy_rot_crop.fits'
EXPTIME = 300.
GAIN = 0.84    ####
NCOMBINE = 12.0
RON = 11.910    ####
m0 = 29.45153 ####
pix2sec = 0.83
fwhm = 5.0
sky_level=0.
sky_subtr=1
xc=620
yc=396
add_to_header(input_image,EXPTIME,GAIN,NCOMBINE,RON,m0,pix2sec,fwhm,sky_level=sky_level,sky_subtr=sky_subtr,xc=xc,yc=yc)
'''