#!/usr/bin/python
# DESCRIPTION:
# Script to prepare galaxy images
# MINIMAL USAGE: python IMP.py [input_image]

from astropy.io import fits as pyfits
import numpy as np
import math
import itertools
import matplotlib.pyplot as plt
from scipy import ndimage
import sys
from matplotlib.path import Path
from itertools import product
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import os
import subprocess
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats
#from photutils import data_properties, properties_table
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from astropy import wcs
from astroquery import ned



def main(input_image, SBiso=24.0, output_mask=None, input_mask=None, m0=20., pix2sec=1.):
    hdulist_data = pyfits.open(input_image)
    data = hdulist_data[0].data
    header = hdulist_data[0].header
    ySize, xSize = data.shape
    
    if input_mask is not None:
        hdulist_mask = pyfits.open(input_mask)
        mask = hdulist_mask[0].data
    else:
        mask = np.zeros(np.shape(data))
    
    
    for k in range(ySize):
        for i in range(xSize):
            if data[k,i]>10**(0.4*(m0-SBiso+5.*math.log10(pix2sec))):
                mask[k,i] = 1.
    if output_mask is None:
        output_mask = 'mask_iso.fits'
        
    outHDU = pyfits.PrimaryHDU(mask, header=header)
    outHDU.writeto(output_mask, clobber=True)     
    print('Done!')

main('galaxy.fits', SBiso=24.0, output_mask='mask_iso.fits', input_mask='mask.fits', m0=28.5839, pix2sec=0.833)
#28.5839	0.0233	6.2319	0.833