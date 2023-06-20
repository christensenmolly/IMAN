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








def main(core_psf, extended_psf, R_core_max, hdu_core=0, hdu_ext=1, output_image='extended_psf.fits', verbosity=True, Radius=None):
    if verbosity: print('Creating extended PSF...')
    azimProfile.main(core_psf, output_model='psf_model.fits', azim_tab='azim_model.txt', mask_image=None, xcen=-1, ycen=-1, ell=0., posang=0., sma_min=-1, sma_max=-1, step=1., sigma_sky=None, sigma_cal=None, outside_frame=False, verbosity=verbosity)

    a,inten,inten_err = np.loadtxt('azim_model.txt', usecols=[0,1,2], dtype=float, unpack=True, skiprows = 1, delimiter='\t')
    if np.max(a)<R_core_max:
        if verbosity: print('WARNING! R_core_max is larger than the maximum radius of the fitting. We set R_core_max as max(Radius)!')
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
    I_R_core_max = np.mean([data_core[yc_core,xc_core-R_core_max],data_core[yc_core+R_core_max,xc_core],data_core[yc_core,xc_core+R_core_max],data_core[yc_core-R_core_max,xc_core]])
    data_core = data_core / I_R_core_max

    # Normalize ext at R_core_max
    I_R_ext_max = np.mean([data_ext[yc_ext,xc_ext-R_core_max],data_ext[yc_ext+R_core_max,xc_ext],data_ext[yc_ext,xc_ext+R_core_max],data_ext[yc_ext-R_core_max,xc_ext]])
    data_ext = data_ext / I_R_ext_max    
    
    for y_ext in range(yc_ext-R_core_max, yc_ext+R_core_max+1):
        for x_ext in range(xc_ext-R_core_max, xc_ext+R_core_max+1):
            if math.sqrt((xc_ext-x_ext)**2+(yc_ext-y_ext)**2)<=R_core_max:
                i = x_ext - xc_ext
                j = y_ext - yc_ext
                x_core = xc_core + i
                y_core = yc_core + j
                data_ext[y_ext,x_ext] = data_core[y_core,x_core]
    data_ext = np.nan_to_num(data_ext)
    data_ext = data_ext / np.nansum(data_ext)
    
    if Radius is not None:
        final_data = data_ext[yc_ext-int(Radius):yc_ext+int(Radius),xc_ext-int(Radius):xc_ext+int(Radius)]
        final_data = final_data / np.nansum(final_data)
    else:
        final_data = data_ext
        
    outHDU = pyfits.PrimaryHDU(final_data)#, header=header_ext)
    outHDU.writeto(output_image, overwrite=True)  
    os.remove('azim_model.txt')
    if verbosity: print('Done!')

#main('psf_r.fits', '/home/amosenko/MyLocalDatabase/PSFs/SDSS/psf-filter-r-v0.4.fits', 9, output_image='extended_psf.fits', hdu_ext=1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Masking")
    parser.add_argument("core_psf", help="Input core fits image")
    parser.add_argument("extended_psf", help="Input extended fits image")
    parser.add_argument("R_core_max", help="Input maximum radius for the core")    
    parser.add_argument("output_image", help="Output image")
    parser.add_argument("--R", nargs='?', const=1, help="Optional: Radius of the output PSF", type=float, default=None) 

    args = parser.parse_args()

    core_psf = args.core_psf
    extended_psf = args.extended_psf
    R_core_max = float(args.R_core_max)
    output_image = args.output_image
    Radius = args.R
    
    
    
    main(core_psf, extended_psf, R_core_max, output_image=output_image, Radius=Radius)