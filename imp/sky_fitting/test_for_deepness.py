#!/usr/bin/python
# DESCRIPTION:
# Script to estimate the limiting SB for a given aperture. Also, it calculates std and mean (median) values of the background.
# MINIMAL USAGE: python test_for_deepness.py [input_image] [mask_image] [ZP] [scale]


# Import the necessary modules
import astropy.io.fits as pyfits
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
from astropy import coordinates
import astropy.units as u
from astropy import wcs
import collections
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clip


from astropy.stats import sigma_clipped_stats
from random import randint


def sky_in_boxes(input_image, m0, pix2sec, mask_image=None, box_size_arcsec=10., Nboxes=1000, n_sigma=3, units='mag', upper=False, verbosity=True):
    hdulist_data = pyfits.open(input_image)
    data = hdulist_data[0].data       
    ny,nx = np.shape(data)  
    
    #subprocess.call("ds9 %s" % (mask_image), shell=True)
    if mask_image is not None:
        hdulist_mask = pyfits.open(mask_image)
        #if len(hdulist_mask)>1:
        #    mask = hdulist_mask[1].data  #### WARNING: Why 1????
        #else:
        mask = hdulist_mask[0].data     
        
    else:
       mask = np.zeros(shape=(ny,nx)) 
    

    N = 0
    STD = []; MEAN = []; MEDIAN = []

    box_size = int(box_size_arcsec/pix2sec)

    while N<Nboxes:
       xc = randint(int(box_size/2.),int(nx-box_size-1))
       yc = randint(int(box_size/2.),int(ny-box_size-1))
       x = range(xc-int(math.ceil(box_size/2.)),xc+int(math.ceil(box_size/2.)),1)
       y = range(yc-int(math.ceil(box_size/2.)),yc+int(math.ceil(box_size/2.)),1)
       I = []
       for k in range(len(y)):
         for i in range(len(x)):
           X = x[i]; Y = y[k]
           if mask[Y,X] == 0.:
             I.append(data[Y,X])
       try:
         mean, median, std = sigma_clipped_stats(I, sigma=3.0)
         STD.append(std)
         MEAN.append(mean)
         MEDIAN.append(median)
         N = N + 1
       except:
         z=1
    
    mean_STD, median_STD, std_STD = sigma_clipped_stats(STD, sigma=3.0)
    mean_MEAN, median_MEAN, std_MEAN = sigma_clipped_stats(MEAN, sigma=3.0)
    mean_MEDIAN, median_MEDIAN, std_MEDIAN = sigma_clipped_stats(MEDIAN, sigma=3.0)
    
    if units=='mag':
        limit = m0-2.5*math.log10(n_sigma*median_STD/(pix2sec**2 * math.sqrt(box_size*box_size)))
        if verbosity: print('SB limit within %ix%i arcsec^2 at %i*sigma: %.2f mag arcsec-2' % (int(box_size_arcsec),int(box_size_arcsec),n_sigma,limit))
    else:
        limit = n_sigma*median_STD/( math.sqrt(box_size*box_size))
        if verbosity: print('SB limit at %ix%i arcsec^2 %i*sigma: %.2f ADU' % (int(box_size_arcsec),int(box_size_arcsec),n_sigma,limit))
    
    #if verbosity: print(limit, median_MEAN, median_MEDIAN, std_MEDIAN, median_STD) 
    
    if upper==False:
        return limit, median_MEAN, median_MEDIAN, std_MEDIAN, median_STD 
    else:
        Imax = 0.
        for k in range(ny):
            for i in range(nx):
                 if mask[k,i]==0. and data[k,i]>Imax:
                        Imax = data[k,i]

        if units=='mag':
            Imax =  m0-2.5*math.log10(Imax/pix2sec**2)
        
        return limit, Imax, median_MEAN, median_MEDIAN, std_MEDIAN, median_STD 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Test for deepness of the image")
    parser.add_argument("input_image", help="Input image")    
    parser.add_argument("m0", help="Zero-point in mag/arcsec^2", type=float)
    parser.add_argument("pix2sec", help="Pixel scale [arcsec/pix]", type=float)

    parser.add_argument("--mask", help="Input mask image", type=str, default=None)
    parser.add_argument("--box_size", help="Optional: Box size in arcsec", type=int, default=30) 
    parser.add_argument("--nboxes", help="Optional: Number of boxes", type=int, default=1000) 
    parser.add_argument("--nsigma", help="Optional: Number of sigma", type=int, default=3) 

    args = parser.parse_args()

    input_image = args.input_image
    mask_image = args.mask
    m0 = args.m0
    pix2sec = args.pix2sec
    box_size_arcsec = args.box_size
    nboxes = args.nboxes
    nsigma = args.nsigma
    
    sky_in_boxes(input_image, m0, pix2sec, mask_image=None, box_size_arcsec=box_size_arcsec, Nboxes=nboxes, n_sigma=nsigma)
  
  
  
