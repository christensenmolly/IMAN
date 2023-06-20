#!/usr/bin/python

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
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from astropy import wcs
from astroquery import ned
import pickle
import collections
import glob

LOCAL_DIR = "/imp/Pipelines"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'detect_objects'))
sys.path.append(os.path.join(IMAN_DIR, 'decomposition/simple_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/masking'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/sky_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))
sys.path.append(os.path.join(IMAN_DIR, 'Ellipse_photometry'))
sys.path.append(os.path.join(IMAN_DIR, 'plotting/1dprofile'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/psf'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/add_astrometry'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/phot_calibration'))
sys.path.append(os.path.join(IMAN_DIR, 'imp'))
sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))
sys.path.append(os.path.join(IMAN_DIR, 'sextractor'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/ds9_regions'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/misc'))

import replace_nans
import mask_nans
import merge_masks
import sersic_fitting
import get_galaxy_center_ned
from IMP import IMP_main
import crop_nan_values
import crop_image
import crop_galaxy_image

def get_galaxy_center(input_image, galaxy_name=None, RA=None, DEC=None):
        # Determine galaxy center:
        xc,yc = get_galaxy_center_ned.main(input_image, name=galaxy_name, RA=RA, DEC=DEC)
        return xc,yc





def main(input_image, sigma_image, galaxy_name, ZP, pixel_scale=0.6, sky_degree=1, user_interact = True, verbosity = False):
    '''
    # Initial cropping
    xmin,ymin,xmax,ymax = crop_nan_values.main(input_image, input_image.split('.fits')[0]+'_cropped.fits')
    input_image = input_image.split('.fits')[0]+'_cropped.fits'

    crop_image.main(sigma_image, xmin, ymin, xmax, ymax, output_image=sigma_image.split('.fits')[0]+'_cropped.fits', hdu=0)
    replace_nans.main(sigma_image.split('.fits')[0]+'_cropped.fits',sigma_image.split('.fits')[0]+'_cropped_no_nan.fits', 0.)
    sigma_image = sigma_image.split('.fits')[0]+'_cropped_no_nan.fits'  


    # General mask
    mto_gain = IMP_main(input_image, 'image_masking', region_file = 'general_mask.reg', mask_image='general_mask.fits', user_interact=user_interact, verbosity=verbosity, masking_tool='sextractor')

    # Convert nan to 0 and add its mask
    mask_nans.main(input_image, 'mask_nonans.fits', float('nan'))
    merge_masks.main(['mask_nonans.fits','general_mask.fits'], 'final_nan_mask.fits')
    
    # Sky fitting
    output_image, mean, median, std, sky_degree = IMP_main(input_image, 'sky_correction', mask_image='final_nan_mask.fits', output_image='galaxy.fits', sky_degree=sky_degree, user_interact=user_interact, verbosity=verbosity) # NOTE: sky_degree!!!

    # Final mask   
    open('galaxy_ellipse.reg', 'a').close()
    if user_interact:
            subprocess.call("ds9 %s -scale histequ -regions %s" % ('galaxy.fits', 'galaxy_ellipse.reg'), shell=True)
        
    # Do final masking: inner and outer
    replace_nans.main('galaxy.fits', 'galaxy_no_nan.fits', 0.)

    IMP_main('galaxy_no_nan.fits', 'final_masking', region_file='galaxy_ellipse.reg', psf_image='psf_060.fits', output_image='galaxy_clean.fits', mask_image='mask.fits',sky_std=std, user_interact=user_interact, verbosity=verbosity) # 'galaxy_mask.reg' will be created        
  

    # PSF for 0.6

    # Convert nan to 0 and add its mask
    #mask_nans.main('galaxy.fits', 'mask_nonans.fits', float('nan'))
    merge_masks.main(['mask_nonans.fits','mask.fits'], 'final_mask.fits')

    #replace_nans.main('galaxy.fits', 'galaxy_no_nan.fits', 0.)
    
    xc,yc = get_galaxy_center('galaxy.fits', galaxy_name)

    xc1,yc1,mtot,re,n,q,PA,C0 = sersic_fitting.main('galaxy_no_nan.fits', sigma_image, 'psf_060.fits', 'final_mask.fits', ZP, pixel_scale, xc=xc, yc=yc, sampling=1, output_dir = './fitting', C0=None, verbose=False, initial_pars=None)
    '''      
    xc2, yc2, sma, smb, PA = crop_galaxy_image.read_region('galaxy_ellipse.reg')
   
    xc,yc = get_galaxy_center('galaxy.fits', galaxy_name)
    
    IMP_main('galaxy_no_nan.fits', 'ellipse_fitting', xc=xc, yc=yc, sma=10., smb=10.*0.9, PA=8, mask_image='mask.fits', m0=ZP, scale=pixel_scale, step=0.03, linear='no', user_interact=user_interact, verbosity=verbosity)
    
    print('Done!')


#main('NGC0925_Spitzer_3.6.fits', 'NGC0925_Spitzer_3.6_Error.fits', 'NGC925', 8.90, pixel_scale=0.6, sky_degree=1, user_interact = True, verbosity = False)
#main('NGC6946_Spitzer_3.6.fits', 'NGC6946_Spitzer_3.6_Error.fits', 'NGC6946', 8.90, pixel_scale=0.6, sky_degree=1, user_interact = True, verbosity = False)
main('IC0342_Spitzer_3.6.fits', 'IC0342_Spitzer_3.6_Error.fits', 'IC0342', 8.90, pixel_scale=0.6, sky_degree=0, user_interact = True, verbosity = False)

#main('NGC2403_Spitzer_3.6.fits', 'NGC2403_Spitzer_3.6_Error.fits', 'NGC2403', 8.90, pixel_scale=0.6, sky_degree=1, user_interact = True, verbosity = False)

#IMP_main('galaxy_no_nan.fits', 'final_masking', region_file='galaxy_ellipse.reg', psf_image='psf_060.fits', output_image='galaxy_clean.fits', mask_image='mask.fits',sky_std=4.48955e-06, user_interact=True, verbosity=False) # 'galaxy_mask.reg' will be created   4.48955e-06



#IMP_main('galaxy_no_nan.fits', 'final_masking', region_file='galaxy_ellipse.reg', psf_image='psf_060.fits', output_image='galaxy_clean.fits', mask_image='mask.fits',sky_std=2.04922e-06, user_interact=True, verbosity=False, masking_tool='combined') # 'galaxy_mask.reg' will be created   4.48955e-06