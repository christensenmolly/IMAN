
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

LOCAL_DIR = "/misc_funcs"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]
sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))
import rebin_image

def swarp_name(verbosity):
        # Check what name has SWarp package on this system
        rCode = subprocess.call("which swarp >/dev/null", shell=True)
        if rCode == 0:
            swarpName = "swarp"
        else:
            rCode = subprocess.call("which SWarp >/dev/null", shell=True)
            if rCode == 0:
                swarpName = "SWarp"
            else:
                if verbosity: print("\033[31m Error: SWarp was not found on your system.\033[0m")
                if verbosity: print("\033[31m The command has to be either 'swarp' or 'SWarp'\033[0m")
                if verbosity: print("\033[31m Intall SWarp package or try to run this script without -s option.\033[0m")
                exit(1)
        return swarpName

def main(input_images, output_image='averaged.fits', verbosity=True):
                        # SWarping the image to reduce all image to a similar wcs:
                        if verbosity: print('SWarping to fix the image...')
                        swarpName = swarp_name(verbosity)
                        callSt = "%s -BACK_TYPE MANUAL -COMBINE_TYPE AVERAGE -VERBOSE_TYPE QUIET -BACK_DEFAULT 0.0 " % (swarpName)
                        callSt += " ".join(["%s[0]" % (s) for s in input_images])
                        subprocess.call(callSt, shell="True")
                        if verbosity: print('Done!')
                        
                        # Rename and remove SWarp tmp files 
                        shutil.move('coadd.fits', output_image)
                        os.remove('coadd.weight.fits')
                        os.remove('swarp.xml')
                        
                        imageHDU = pyfits.open(output_image, mode='update')
                        header = imageHDU[0].header
                        if 'GAIN' in header:
                            if header['GAIN'] == 0.:
                                header['GAIN'] = 10000.
                                
                        if 'EXPTIME' in header:
                            header['EXPTIME'] = float(header['EXPTIME'])/float(len(input_images))
                        
                        imageHDU.flush()
'''                        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="WCS correction with SWarp")
    parser.add_argument("input_images", help="Input images separated by comma")    
    parser.add_argument("output_image", help="Output image")
    
    args = parser.parse_args()

    input_images = args.input_images.split(',')
    output_image = args.output_image

    
    main(input_images, output_image, verbosity=True)
'''

for band in ['g','r','i']:
    image = 'NGC4469_0_%s.fits' % (band)
    input_images = [image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,'NGC4469_%s_trim.fits' % (band)]
    main(input_images, output_image='averaged.fits', verbosity=True)
    
    image = 'NGC4469_2_%s.fits' % (band)
    input_images = [image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,'averaged.fits']
    main(input_images, output_image='averaged1.fits', verbosity=True)

    image = 'NGC4469_3_%s.fits' % (band)
    input_images = [image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,image,'averaged1.fits']
    main(input_images, output_image='averaged2.fits', verbosity=True)


    rebin_image.rebin('NGC4469_%s_trim.fits' % (band), 'averaged2.fits', output_image='NGC4469_%s_trim_final.fits' % (band), hdu_ref=0, hdu_inp=0, preserve_bad_pixels=True, no_interp=False)