#!/usr/bin/python
# -*- coding:  cp1251 -*-
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
import re
import glob
from astropy.io import fits as pyfits
import collections
import warnings
from types import SimpleNamespace
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like

from astropy.io import fits
from photutils.isophote import IsophoteList
from photutils.isophote import build_ellipse_model
from astropy.stats import sigma_clipped_stats
import argparse

warnings.filterwarnings("ignore")
tmp_out = sys.stdout
FNULL = open(os.devnull, 'w')

LOCAL_DIR = "/imp/psf"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'iraf_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/misc'))
sys.path.append(os.path.join(IMAN_DIR, 'Ellipse_photometry'))
mpl.style.use('classic')

import iraf_ellipse
import arithm_operations
import azimProfile
import compare_iraf_and_ellipse
import build_iraf_model
import norm_image
from scipy.signal import savgol_filter
import fit_star


def ds9_to_np(x):
    '''
    Function to convert ds9 pixel coordinates to numpy format, i.e
    0.5 pix is 0, 1.5 pix is 1
    '''
    return int(round(x)) - 1

def get_ellipse_coords_at_sma(x0, y0, sma, pa):
    pa = np.radians(pa)
    x = sma*math.cos(pa) + x0
    y = sma*math.sin(pa) + y0
    return ds9_to_np(x), ds9_to_np(y)


def ellipse_mask(X, Y, xc, yc, ellA, ellB, ellPA, xSize, ySize):
        cen = PPoint(xc, yc)
        cospa = cos(radians(ellPA))
        sinpa = sin(radians(ellPA))
        # Check if whole ellipse is inside of the image
        # and obtain size of the ellipse in xy plane
        xMax = 0.0
        yMax = 0.0
        xMin = 1e10
        yMin = 1e10
        for e in linspace(0, 4*pi, 1000):
            cose = cos(e)
            sine = sin(e)
            x = xc + ellA * cose * cospa - ellB * sine * sinpa
            y = yc + ellB * sine * cospa + ellA * cose * sinpa
            if x > xMax:
                xMax = x
            if y > yMax:
                yMax = y
            if x < xMin:
                xMin = x
            if y < yMin:
                yMin = y
        xMin = max(0, int(round(xMin)))
        xMax = min(xSize, int(round(xMax)))
        yMin = max(0, int(round(yMin)))
        yMax = min(ySize, int(round(yMax)))
        focusR = (ellA ** 2.0 - ellB ** 2.0) ** 0.5
        focus10 = PPoint(xc + focusR, yc)  # Unrotated
        focus20 = PPoint(xc - focusR, yc)  #
        focus1 = rot_point(focus10, cen, radians(ellPA))
        focus2 = rot_point(focus20, cen, radians(ellPA))
        # Find pixels inside of the ellipse
        dEll = 2 * ellA

        if True:
            for x in range(xMin, xMax+1):
               for y in range(yMin, yMax+1):
                    dFocus1 = hypot(x-focus1.x, y-focus1.y)
                    dFocus2 = hypot(x-focus2.x, y-focus2.y)
                    dPoint = dFocus1 + dFocus2
                    if dPoint < dEll and x-1==X and y-1==Y:
                      return True
        return False
    
class PPoint:
    def __init__(self, x, y):
        self.x = x
        self.y = y


def rot_point(p, orig, angle):
    """ performs rotation of a point (x,y) around point (x0,y0)
        by angle radians"""
    x1 = cos(angle)*(p.x-orig.x) - sin(angle)*(p.y-orig.y) + orig.x
    y1 = sin(angle)*(p.x-orig.x) + cos(angle)*(p.y-orig.y) + orig.y
    return PPoint(x1, y1)


def main(core_psf, extended_psf, R_core_max, hdu_core=0, hdu_ext=0, output_image='extended_psf.fits'):
    norm_image.main(core_psf, 'tmp_core.fits')
    os.rename('tmp_core.fits', core_psf)

    norm_image.main(extended_psf, 'tmp_ext.fits')
    os.rename('tmp_ext.fits', extended_psf)    
    
    mag,FWHM,beta,q,PA = fit_star.main(core_psf, window='gauss', remove_model=False)
    
    
    azimProfile.main(core_psf, output_model='psf_model.fits', azim_tab='azim_model.txt', mask_image=None, xcen=-1, ycen=-1, ell=1.-q, posang=PA-90., sma_min=-1, sma_max=-1, step=1., sigma_sky=None, sigma_cal=None, outside_frame=False)
    
    a,inten,inten_err = np.loadtxt('azim_model.txt', usecols=[0,1,2], dtype=float, unpack=True, skiprows = 1, delimiter='\t')
    if np.max(a)<R_core_max:
        print('WARNING! R_core_max is larger than the maximum radius of the fitting. We set R_core_max as max(Radius)!')
        R_core_max = np.max(a)
    
    core_psf = 'psf_model.fits'
    
    R_core_max = int(R_core_max)
    
    # Read in the input PSFs
    hdulist_core = pyfits.open(core_psf)
    data_core = hdulist_core[hdu_core].data
    ny_core,nx_core = np.shape(data_core)
    xc_core = int(nx_core/2. + 0.5) - 1
    yc_core = int(ny_core/2. + 0.5) - 1

    
    hdulist_ext = pyfits.open(extended_psf)
    data_ext = hdulist_ext[hdu_ext].data # WARNING
    header_ext = hdulist_ext[hdu_ext].header # WARNING
    ny_ext,nx_ext = np.shape(data_ext)
    xc_ext = int(nx_ext/2. + 0.5) - 1
    yc_ext = int(ny_ext/2. + 0.5) - 1
    
    
    # Normalize core at R_core_max
    x_norm,y_norm = get_ellipse_coords_at_sma(nx_core/2., ny_core/2., R_core_max, PA-90.)

    
    I_R_core_max = data_core[y_norm,x_norm]
    data_core = data_core / I_R_core_max

    # Normalize ext at R_core_max
    x_norm,y_norm = get_ellipse_coords_at_sma(nx_ext/2., ny_ext/2., R_core_max, PA-90.)
    I_R_ext_max = data_ext[y_norm,x_norm]
    data_ext = data_ext / I_R_ext_max
    
    
    
    
    for y_ext in range(yc_ext-R_core_max, yc_ext+R_core_max+1):
        for x_ext in range(xc_ext-R_core_max, xc_ext+R_core_max+1):
            if math.sqrt((xc_ext-x_ext)**2+(yc_ext-y_ext)**2)<=R_core_max:
                i = x_ext - xc_ext
                j = y_ext - yc_ext
                x_core = xc_core + i
                y_core = yc_core + j
                if ellipse_mask(x_ext, y_ext, xc_ext, yc_ext, R_core_max, R_core_max*q, PA-90., nx_ext, ny_ext):
                    data_ext[y_ext,x_ext] = data_core[y_core,x_core]
    
    
    data_ext = np.nan_to_num(data_ext)
    data_ext = data_ext / np.nansum(data_ext)
    
    
    outHDU = pyfits.PrimaryHDU(data_ext)#, header=header_ext)
    outHDU.writeto(output_image, overwrite=True)  
    os.remove('azim_model.txt')


#main('psf.fits', 'extended_psf.fits', 10, hdu_core=0, hdu_ext=0, output_image='extended_psf_new_core.fits')

#main('psf_r.fits', '/home/amosenko/MyLocalDatabase/PSFs/SDSS/psf-filter-r-v0.4.fits', 9, output_image='extended_psf.fits', hdu_ext=1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Masking")
    parser.add_argument("core_psf", help="Input core fits image")
    parser.add_argument("extended_psf", help="Input extended fits image")
    parser.add_argument("R_core_max", help="Input maximum radius for the core")    
    parser.add_argument("output_image", help="Output image")

    args = parser.parse_args()

    core_psf = args.core_psf
    extended_psf = args.extended_psf
    R_core_max = float(args.R_core_max)
    output_image = args.output_image
    
    
    
    main(core_psf, extended_psf, R_core_max, output_image=output_image)