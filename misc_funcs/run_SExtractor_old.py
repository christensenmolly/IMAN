#!/usr/bin/python
# DESCRIPTION:
# Script to run SExtractor.
# MINIMAL USAGE: python sextr.py [input_image]

import sys
import math
import numpy as np
from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import NullFormatter
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
from astropy.io import fits as pyfits
import argparse

tmp_out = sys.stdout

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))

# Colors to highlight the output text
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

def read_input_files(se_file):
    with open(se_file,'r') as f:
        lines = f.readlines()
    f.close()
    files = []
    for line in lines:
      if 'PARAMETERS_NAME' in line:
        files.append(line.split()[1])
      if 'FILTER_NAME' in line:
        files.append(line.split()[1])
      if 'STARNNW_NAME' in line:
        files.append(line.split()[1])
    return files

def run_sextr(file_in,m0,GAIN,pix2sec,fwhm,se_file='default.sex'):
    # Check how to call Sextractor:
    if subprocess.call("which sex >/dev/null", shell=True) == 0:
        callSE = "sex %s " % file_in
    elif subprocess.call("which sextractor >/dev/null", shell=True) == 0:
        callSE = "sextractor %s " % file_in
    else:
        print(bcolors.FAIL+ "SExtractor was not found. Exiting..." + bcolors.ENDC)
        exit(1)
    
    # Copy input default files for Sextractor
    files = read_input_files(PATH_TO_SCRIPT+'/'+se_file)
    files.append('default.psf')
    for file in files:
      shutil.copy(PATH_TO_SCRIPT+'/'+file,file)

    print(bcolors.OKBLUE+ "Running Sextractor" + bcolors.ENDC)


    # Run Sextractor:
    callSE += "-MAG_ZEROPOINT " + str(m0) + " -GAIN " + str(GAIN) + " -PIXEL_SCALE " + str(pix2sec) + " -SEEING_FWHM " + str(fwhm)   
    callSE += " -c %s" % (PATH_TO_SCRIPT+'/'+se_file)   
    callSE += " -VERBOSE_TYPE=QUIET"

    subprocess.call(callSE, shell=True)
    
    # Remove copied files
    for file in files:
      os.remove(file)
    print('Done!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Isophote analysis")
    parser.add_argument("inputImage", help="Input fits image")
    parser.add_argument("ZeroPoint", nargs='?', const=1, help="Input Zero point in [mag arcsec^-2] (optional)",type=float,default=20.0)
    parser.add_argument("Gain", nargs='?', const=1, help="Input gain in [e-/ADU] (optional)",type=float,default=4.0)
    parser.add_argument("Scale", nargs='?', const=1, help="Input scale in [arcsec/pix] (optional)",type=float,default=1.0)
    parser.add_argument("FWHM", nargs='?', const=1, help="Input PSF FWHM (optional)",type=str,default=1.0)
    parser.add_argument("seFile", nargs='?', const=1, help="Input sextractor input file (optional)",type=str,default='default.sex') 

    args = parser.parse_args()

    input_image = args.inputImage
    m0 = args.ZeroPoint
    gain = float(args.Gain)
    pix2sec = float(args.Scale)
    FWHM = float(args.FWHM)
    seFile = args.seFile

    run_sextr(input_image,m0,gain,pix2sec,FWHM,se_file=seFile)   
