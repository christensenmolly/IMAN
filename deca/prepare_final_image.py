#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **       DECA -- DEComposition Analysis of galaxy images       **
# **         © Astronomical Observatory, Ghent University        **
# *****************************************************************


# Import standard modules
import pylab
import sys
import os
import shutil
import math
import numpy as np
import scipy as sp
from numpy import *
from pylab import *
import subprocess
import argparse
import pyfits

import add_keyw_to_header
import cosmo

# -----------------------------------------------------------------
# FUNCTION TO ADD NEW KEYWORDS WHICH ARE REQUIRED BY DECA
# AND CORRECT FOR GALACTIC EXTINCTION AND K-CORRECTION
def main(input_image,EXPTIME,GAIN,NCOMBINE,
	 RON,m0,pix2sec,fwhm,Aext,filter_name,colour_name,colour_value,redshift,output_image=None):
  #print 'here'
  #exit()
  if output_image==None:
    output_image = input_image
  
  hdulist = pyfits.open(input_image)
  image_data = hdulist[0].data
  primhdr = hdulist[0].header

  if np.isnan(Aext)==False:
    # Correct for Galactic extinction
    image_data = image_data * 10**(0.4*Aext) 
    print 'Image has been corrected for Galactic extinction!'  

  if np.isnan(colour_value)==False and np.isnan(redshift)==False:
    # K-correction 
    image_data = image_data * cosmo.A_kcor(filter_name, redshift, colour_name, colour_value)
    print 'Image has been K-corrected!' 

  hdu = pyfits.PrimaryHDU(image_data,primhdr)
  hdu.writeto(output_image,clobber=True)

  # Add these keywords to the image header
  add_keyw_to_header.add_to_header(output_image,EXPTIME,GAIN,NCOMBINE,RON,m0,pix2sec,fwhm,0.,1)
  print 'New keywords have been added to the header!'
# -----------------------------------------------------------------



if __name__ == '__main__':
    print 'FUNCTION TO ADD NEW KEYWORDS WHICH ARE REQUIRED BY DECA AND CORRECT FOR GALACTIC EXTINCTION AND K-CORRECTION'
    
    input_image = str(raw_input('Please enter the fits image'))
    EXPTIME = float(raw_input('Please enter the exposition time'))
    GAIN = float(raw_input('Please enter the gain'))
    NCOMBINE = float(raw_input('Please enter the number of combined images'))
    RON = float(raw_input('Please enter read-out-noise value'))
    m0 = float(raw_input('Please enter the Zero point'))
    pix2sec = float(raw_input('Please enter the scale'))
    fwhm = float(raw_input('Please enter the FWHM'))
    Aext = float(raw_input('Please enter the extinction (%s)'% (str('nan')))) or str('nan') 
    filter_name = str(raw_input('Please enter the photometric band (%s)'% (str('g')))) or str('g') 
    colour_name = str(raw_input('Please enter the colour (%s)'% (str('g - r')))) or str('g - r') 
    colour_value = float(raw_input('Please enter the colour_value (%s)'% (str('nan')))) or str('nan') 
    redshift = float(raw_input('Please enter the redshift (%s)'% (str('nan')))) or str('nan')    
    output_image = str(raw_input('Please enter the name for the output fits image (%s)'% (str(input_image)))) or str(input_image) 
    
    main(input_image,EXPTIME,GAIN,NCOMBINE,
	 RON,m0,pix2sec,fwhm,Aext,filter_name,colour_name,colour_value,redshift,output_image)
  