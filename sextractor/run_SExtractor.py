#!/usr/bin/python
# DESCRIPTION:
# Script to run SExtractor.
# MINIMAL USAGE: python run_SExtractor.py [input_image]

import sys
import subprocess
from pylab import *
import itertools
import os
from os.path import exists
from os import remove
from scipy.spatial import cKDTree
from scipy.optimize import fmin_tnc, fmin
#import arithm_operations
import argparse
from astropy.modeling import models, fitting
from datetime import datetime
import shutil
# Disable astropy logging except for warnings and errors
from astropy import log
from astropy.io import fits as pyfits
log.setLevel("WARNING")

import warnings
warnings.filterwarnings("ignore")

def fast_background(input_image):
    from astropy.stats import sigma_clipped_stats
    
    hdulist = pyfits.open('segm.fits') # open a FITS file
    mask_infr = hdulist[0].data
    mask_astropy = (mask_infr>0.)

    hdulist = pyfits.open(input_image) # open a FITS file
    scidata = hdulist[0].data	
    ny,nx = np.shape(scidata)

    mean, median, std = sigma_clipped_stats(scidata, sigma=3.0, iters=5, mask=mask_astropy)
    return median, std   


def call_SE(input_image, snr=None, min_pix=None, sextr_dir=None, sextr_setup='cold.sex', sextr_param='default.param', output_cat='field.cat', checkimage_type='SEGMENTATION',checkimage_name='segm.fits', sextr_add_string=None,determine_sky=False, m0=28.0, GAIN=10., pix2sec=1., fwhm=4., verbosity=True):
    '''
    Function to run SExtractor
    hot.sex, cold.sex, hot_objects.sex, cold_objects.sex, default.sex
    '''
    if verbosity: print('Running SExtractor...')
    
    # detect the name of SExtractor
    if subprocess.call("which sex >/dev/null", shell=True) == 0:
        callString = "sex %s " % input_image
        CallString = "sex %s " % input_image
    elif subprocess.call("which sextractor >/dev/null", shell=True) == 0:
        callString = "sextractor %s " % input_image
        CallString = "sextractor %s " % input_image
    elif subprocess.call("which source-extractor >/dev/null", shell=True) == 0:
        callString = "source-extractor %s " % input_image
        CallString = "source-extractor %s " % input_image
    else:
        if verbosity: print("SExtractor was not found. Exiting...")
        exit(1)

    if sextr_dir is None:
        sextr_dir = os.path.abspath(os.path.dirname(__file__) + '/sex_lib')
        if not os.path.exists(sextr_dir):
            sextr_dir = os.path.abspath(os.path.abspath(__file__).split('/run_SExtractor.py')[0] + '/sex_lib')
    

    callString += "-c %s/%s -CATALOG_NAME %s " % (sextr_dir, sextr_setup, output_cat)
    if snr is not None:
        callString += " -DETECT_THRESH %f -ANALYSIS_THRESH %f " % (snr,snr)
    
    if min_pix is not None:
        callString += "-DETECT_MINAREA %f " % (min_pix)
    
    callString += " -PARAMETERS_NAME %s/%s -FILTER_NAME %s/default.conv -STARNNW_NAME %s/default.nnw -CHECKIMAGE_TYPE %s -CHECKIMAGE_NAME %s " % (sextr_dir,sextr_param,sextr_dir,sextr_dir,checkimage_type,checkimage_name)
    
    callString += " -MAG_ZEROPOINT " + str(m0) + " -GAIN " + str(GAIN) + " -PIXEL_SCALE " + str(pix2sec) + " -SEEING_FWHM " + str(fwhm)   
    
    callString += " -PSF_NAME %s/%s " % (sextr_dir,'default.psf')
    
    callString += " -VERBOSE_TYPE=QUIET"

    if sextr_add_string is not None:
        callString += sextr_add_string
    #print(callString)
    #exit()
    subprocess.call(callString, shell=True)
    if verbosity: print('Done!')
    if determine_sky:
        median, std = fast_background(input_image)
        if verbosity: print('Sky level and rms: %f %f' % (median, std))
        return median, std
    else:
        return float('nan'),float('nan')
    
    
 
    '''
	arithm_operations.main('segm_cold.fits','segm_hot.fits','add','segm.fits')
	filenames = ['field_hot.cat', 'field_cold.cat']
	with open('field.cat', 'w') as outfile:
	    for fname in filenames:
		with open(fname) as infile:
		    for line in infile:
			outfile.write(line)    
	remove('segm_hot.fits')
	remove('segm_cold.fits')
	remove('field_hot.cat')
	remove('field_cold.cat')
    '''

def get_SE_results():
    objects = []
    for line in open("field.cat"):
        if line.startswith("#"):
            continue
        params = line.split()
        n = int(params[0])
        xCen = float(params[1])
        yCen = float(params[2])
        kron = float(params[8])
        ellA = kron * float(params[4]) * Cov_coeff
        ellB = kron * float(params[5]) * Cov_coeff
        PA = float(params[6])
        objects.append(ObjParams(xCen, yCen, ellA, ellB, PA))
    return objects




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky background estimation")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("--snr", nargs='?', const=1, help="Optional: Signal-to-noise ratio of detected pixels",type=float,default=2.) 
    parser.add_argument("--min_pix", nargs='?', const=1, help="Optional: Number of joint pixels",type=int,default=5)
    parser.add_argument("--sextr_dir", nargs='?', const=1, help="Optional:Sextractor directory",type=str,default=None)
    parser.add_argument("--sextr_setup", nargs='?', const=1, help="Optional: Sextractor setup file",type=str,default='cold.sex')
    parser.add_argument("--sextr_param", nargs='?', const=1, help="Optional: Sextractor parameter file",type=str,default='default.param')
    parser.add_argument("--sky", help="Determine sky level", action="store_true", default=False)
    args = parser.parse_args()

    input_image = args.inputImage
    snr = args.snr
    min_pix = args.min_pix
    sextr_dir = args.sextr_dir
    sextr_setup = args.sextr_setup
    sextr_param = args.sextr_param
    determine_sky = args.sky
    
    call_SE(input_image, snr=snr, min_pix=min_pix, sextr_dir=sextr_dir, sextr_setup=sextr_setup, sextr_param=sextr_param, output_cat='field.cat', checkimage_type='SEGMENTATION',checkimage_name='segm.fits', sextr_add_string=None,determine_sky=determine_sky)
