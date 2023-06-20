#!/usr/bin/python
# DESCRIPTION:
# A script to fit sky background in an annulus around galaxy. An outer galaxy ellipse --galaxy_ellipse (where outermost isophotes end) is required, either as a DS9 region file (ellipse region in images coordinates), or in the format xc,yc,sma,smb,PA. The width of the annulus is controlled by the key --annulus_width (in pix). To check the annulus in DS9 before starting the fitting the key --manual should be given. For now, only constant background level (--degree 0) within the annulus is computed.
# MINIMAL USAGE: python sky_around_galaxy.py [input_image] --galaxy_ellipse [region_file OR xc,yc,sma,smb,PA]

from astropy.modeling import models, fitting
import warnings
warnings.filterwarnings("ignore")
from astropy.stats import sigma_clipped_stats
import sys
import subprocess
from astropy.io import fits as pyfits
from pylab import *
import itertools
import os
from os.path import exists
from os import remove
from scipy.spatial import cKDTree
from scipy.optimize import fmin_tnc, fmin
import argparse
from astropy.modeling import models, fitting
from datetime import datetime
import shutil
from astropy import log
import warnings
import numpy as np
from photutils import EllipticalAnnulus
import matplotlib.pyplot as plt
from types import SimpleNamespace
import random
import scipy as sp

LOCAL_DIR = "/imp/sky_fitting"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/masking'))

import convert_reg_to_mask


def main(input_image, region_file, sky_value, output_image='fixed.fits', mask_value=1):
    # Convert region file to mask image
    convert_reg_to_mask.mask(input_image, region_file, output_image=None, output_mask='inner_mask.fits', mask_value=mask_value, show_running=True, mask_DN=None, verbosity=True)


    hdulist = pyfits.open(input_image, do_not_scale_image_data=True, mode='update')
    data = hdulist[0].data 
    header = hdulist[0].header
    ny,nx = np.shape(data)
    
    fix_data = np.copy(data)
    
    # Read in created mask image
    hdulist_mask = pyfits.open('inner_mask.fits', do_not_scale_image_data=True, mode='update')
    mask = hdulist_mask[0].data    
    
    fix_data[mask==mask_value] = data[mask==mask_value] + sky_value

    outHDU = pyfits.PrimaryHDU(fix_data, header=header)
    outHDU.writeto(output_image, overwrite=True)       
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky background estimation")
    parser.add_argument("inputImage", help="Input fits image with the object") 
    parser.add_argument("--o", nargs='?', const=1, help="Optional: Output image with subtracted sky", type=str, default=None) 
    parser.add_argument("--sky", nargs='?', const=1, help="Sky value to be subtracted", type=float, default=0.)
    parser.add_argument("--r", nargs='?', const=1, help="Region file where the sky should be sutracted from", type=str, default=None)
    
    args = parser.parse_args()

    input_image = args.inputImage
    output_image = args.o
    sky_value = args.sky
    region_file = args.r
    
    main(input_image, region_file, sky_value, output_image=output_image)