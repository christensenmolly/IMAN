#!/usr/bin/python
# DESCRIPTION:
# Script to fit galaxy centre

import math
import numpy as np
from astropy.io import fits as pyfits
import argparse
import numpy.ma
import os
import scipy.optimize
import subprocess
import fit_gauss    
import sys

LOCAL_DIR = "/imp/cropping"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'decomposition/simple_fitting'))
import sersic_fitting

    
def ds9_to_np(x):
    '''
    Function to convert ds9 pixel coordinates to numpy format, i.e
    0.5 pix is 0, 1.5 pix is 1
    '''
    return int(round(x)) - 1    

def np_to_ds9(x):
    '''
    Function to convert numpy coordinates to ds9 format, i.e
    0 is between 0.5 and 1.5 (not incl), 1 is between 1.5 (incl) and
    2.5 (not incl)
    '''
    return x+1.    
    
def fit_gaussian(data, xc, yc, mask_data=None, r=5):
        # xc,yc - numpy coordinates
        # r - radius of the region where to find center
        # returned new_xc, new_yc - numpy coordinates
        r = int(r)
        
        if mask_data is not None:
            fit_mask_data = mask_data[yc-r:yc+r, xc-r:xc+r]
        
        fit_data = data[yc-r:yc+r, xc-r:xc+r]
        
        
        res = fit_gauss.FitGauss2D(fit_data, ip=None)
        new_xc = res[0][1]
        new_yc = res[0][2]

        # Convert into np:
        #new_xc = ds9_to_np(new_xc)
        #new_yc = ds9_to_np(new_yc)
        
        # Go back to data:
        new_xc = ds9_to_np(np_to_ds9(xc) - r + new_xc)
        new_yc = ds9_to_np(np_to_ds9(yc) - r + new_yc)

        #print(new_xc,new_yc)
        #exit()

        return new_xc,new_yc

    
    


def main(galaxy_image, mask_image=None, sigma_image=None, psf_image=None, xc=None, yc=None, hdu=0, fit_method='sersic'):
    if True:
        hdulist = pyfits.open(galaxy_image)
        data = hdulist[hdu].data
        if mask_image is not None:
            hdulist_mask = pyfits.open(mask_image)
            mask_data = hdulist_mask[0].data        
        else:
            mask_data = None
        
        nx, ny = data.shape[1], data.shape[0]   
        if xc is None:
            xc = nx/2.
        if yc is None:
            yc = ny/2.

        # Convert coordinates to numpy
        xc = ds9_to_np(xc)
        yc = ds9_to_np(yc)  
        

        if fit_method == 'gauss':
            # Fit galaxy center within a square with a gaussian function
            xc,yc = fit_gaussian(data, xc, yc, mask_data, r=5)
        elif fit_method == 'center':
            I_max = np.max(data[yc-5:yc+5,xc-5:xc+5])
            [yc1,xc1] = list(np.where(data==I_max))
            xc = xc1[0]
            yc = yc1[0]
        elif fit_method=='sersic':
            xc,yc,mtot,re,n,q,PA,C0 = sersic_fitting.main(galaxy_image, sigma_image, psf_image, mask_image, 0.0, 1.0, xc=np_to_ds9(xc), yc=np_to_ds9(yc), sampling=1, output_dir = './fitting')
            print(xc,yc,mtot,re,n,q,PA,C0)
            exit()
            xc = ds9_to_np(xc)
            yc = ds9_to_np(yc)   
    else:
        # Convert coordinates to numpy
        xc = ds9_to_np(xc)
        yc = ds9_to_np(yc)          

