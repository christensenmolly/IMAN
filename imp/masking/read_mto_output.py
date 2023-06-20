# Import the necessary modules
from astropy.io import fits as pyfits
#import pyfits
import numpy as np
import math
import itertools
from scipy import ndimage
import sys
from itertools import product
from matplotlib.path import Path
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import os
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from astropy import wcs


def main(mto_csv, input_image=None):
    # ID,X,Y,A,B,theta,total_flux,mu_max,mu_median,mu_mean,R_fwhm,R_e,R10,R90
    try:
        ID,X,Y,A,B,theta,total_flux,mu_max,mu_median,mu_mean,R_fwhm,R_e,R10,R90 = np.loadtxt(mto_csv, usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13],dtype=float, unpack=True,skiprows=1,delimiter=',')
        ID = np.array(ID, dtype = int)
        Z = len(X)
    except:
        ID,X,Y,A,B,theta,total_flux,mu_max,mu_median,mu_mean,R_fwhm,R_e,R10,R90 = np.loadtxt(mto_csv, usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13],dtype=float, unpack=True,skiprows=1,delimiter=',',ndmin=2)
        ID = np.array(ID, dtype = int)        
    
    
    RA = []; DEC = []
    if input_image is not None:
        hdulist = pyfits.open(input_image) # open FITS file  
        Data = hdulist[0].data
        header = hdulist[0].header
        
        if 'COMMENT' in header:
                        del header['COMMENT']
        if 'HISTORY' in header:
                        del header['HISTORY']
        if '' in header:
                        del header['']  

        w = wcs.WCS(header)
        
        for k in range(len(X)):
            pixcrd = np.array([[X[k],Y[k]]], dtype = float)
            world = w.wcs_pix2world(pixcrd, 1) 
            RA.append(world[0][0])
            DEC.append(world[0][1])
            print(X[k],Y[k],RA[k],DEC[k])
    else:
        for k in range(len(X)):
            RA.append(float('nan'))
            DEC.append(float('nan'))
    
    return ID,X,Y,A,B,theta,total_flux,mu_max,mu_median,mu_mean,R_fwhm,R_e,R10,R90,RA,DEC

#main('parameters.csv', input_image='frame-r-007713-1-0269.fits')
