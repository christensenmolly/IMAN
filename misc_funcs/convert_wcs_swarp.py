
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

def main(input_image, verbosity=True):
                        # SWarping the image to reduce all image to a similar wcs:
                        if verbosity: print('SWarping to fix the image...')
                        swarpName = swarp_name(verbosity)
                        callSt = "%s -BACK_TYPE MANUAL -COMBINE_TYPE AVERAGE -VERBOSE_TYPE QUIET -BACK_DEFAULT 0.0 -INTERPOLATE N -RESAMPLE N " % (swarpName)
                        callSt += " ".join(["%s[0]" % (s) for s in [input_image,input_image]])
                        subprocess.call(callSt, shell="True")
                        if verbosity: print('Done!')
                        
                        # Rename and remove SWarp tmp files 
                        shutil.move('coadd.fits', input_image)
                        os.remove('coadd.weight.fits')
                        os.remove('swarp.xml')
                        
                        imageHDU = pyfits.open(input_image, mode='update')
                        header = imageHDU[0].header
                        if 'GAIN' in header:
                            if header['GAIN'] == 0.:
                                header['GAIN'] = 10000.
                                
                        if 'EXPTIME' in header:
                            header['EXPTIME'] = float(header['EXPTIME'])/2.
                        
                        imageHDU.flush()
                        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="WCS correction with SWarp")
    parser.add_argument("input_image", help="Input image")    

    args = parser.parse_args()

    input_image = args.input_image

    
    main(input_image, verbosity=True)