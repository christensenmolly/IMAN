#! /usr/bin/env python

import sys
import subprocess
import os
import shutil
import argparse
from os.path import exists
from math import hypot, pi, log10
import math
import numpy  as np
from numpy import ma, mean, std, zeros_like, where, arange, exp, zeros
from numpy import copy as npcopy
from numpy import sum as npsum
from pylab import plot, show, savefig, xlabel, ylabel, vlines, clf
import astropy.io.fits as pyfits
from scipy.odr.odrpack import *
import warnings


LOCAL_DIR = "/imp/rotate"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]
sys.path.append(os.path.join(IMAN_DIR, 'sextractor'))

import run_SExtractor

warnings.filterwarnings("ignore")

FNULL = open(os.devnull, 'w')





def get_galaxy_params(fitsFile, xCenField, yCenField):
    """ Function finds object nearest to the center
    of the field."""
    hdu = pyfits.open(fitsFile)
    ySize, xSize = hdu[0].data.shape
    hdu.close()

    minArea = 100  # Minimal area of interested objects [pix^2]
    minCenterDist = 1e10
    for line in open("field.cat"):
        if line.startswith("#"):
            continue
        params = line.split()
        N = int(params[0])
        xCen = float(params[1])
        yCen = float(params[2])
        kron = float(params[8])
        ellA = kron * float(params[4])
        ellB = kron * float(params[5])
        PA = float(params[6])
        ellArea = pi * ellA * ellB
        if ellArea > minArea:
            centerDist = hypot(xCen-xCenField, yCen-yCenField)
            if centerDist < minCenterDist:
                minCenterDist = centerDist
                galN = N
                galXCen = xCen
                galYCen = yCen
                galEllA = ellA
                galEllB = ellB
                galPA = PA
    return galN, galXCen, galYCen, galEllA, galEllB, galPA






def main(input_images, xc, yc, verbosity=True):
    # 1. Launching of the SE for background subtracted images
    # to obtain some geometric parameters of the galaxy
    if verbosity: print('Computing the position angle to rotate the galaxy...')
    if verbosity: print('Rotation is done counterclockwise from the x-axis')
    galaxy_image = input_images[0]
    if len(input_images)>1:
        mask_image = input_images[1]
    else:
        mask_image = None


    run_SExtractor.call_SE(galaxy_image, snr=None, min_pix=None, sextr_dir=None, sextr_setup='default.sex', sextr_param='default.param', output_cat='field.cat', checkimage_type='SEGMENTATION',checkimage_name='segm.fits', sextr_add_string=None,determine_sky=False, m0=28.0, GAIN=10., pix2sec=1., fwhm=4.)

    galN, xCen, yCen, ellA, ellB, ellPA = get_galaxy_params(galaxy_image, xc, yc)

    
        
    # Use max pixels method
    optPosAng = ellPA
    if verbosity: print("SE Position angle = %1.1f" % (optPosAng))
    os.remove('segm.fits')
    os.remove('field.cat')
    return optPosAng
