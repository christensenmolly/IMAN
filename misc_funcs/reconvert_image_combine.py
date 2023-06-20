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
    

def main(initial_image, combined_image, output_image, verbosity=True):
        hdulist0 = pyfits.open(initial_image)
        frame0 = hdulist0[0].data
        header0 = hdulist0[0].header

        hdulist1 = pyfits.open(combined_image)
        frame1 = hdulist1[0].data
        ny,nx = np.shape(frame1)
        
        header0['NAXIS1'] = nx 
        header0['NAXIS1'] = ny
        
        hdu = pyfits.PrimaryHDU(frame1, header0)
        hdu.writeto('tmp.fits', overwrite=True) 
        os.rename('tmp.fits', output_image)
                    
                        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Re-wcs")
    parser.add_argument("initial_image", help="Initial image before stacking")    
    parser.add_argument("combined_image", help="Combined image (returned from imcombine)")
    parser.add_argument("output_image", help="Output image") 
    
    args = parser.parse_args()

    initial_image = args.initial_image
    combined_image = args.combined_image
    output_image = args.output_image
    
    
    main(initial_image, combined_image, output_image, verbosity=True)