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








def main(core_psf, extended_psf_profile, R_core_max, output_model='extended_psf.fits', core_fit=True):
    # Read in the input PSFs
    hdulist_core = pyfits.open(core_psf)
    data_core = hdulist_core[0].data
    nx_core,ny_core = data_core.shape[1], data_core.shape[0]
    
    #hdulist_ext = pyfits.open(extended_psf)
    #data_ext = hdulist_ext[0].data
    #nx_ext,ny_ext = data_ext.shape[1], data_ext.shape[0]


    
    
    azimProfile.main(core_psf, output_model=None, azim_tab='azim_model_core.txt', mask_image=None, xcen=-1, ycen=-1, ell=0., posang=0., sma_min=-1, sma_max=R_core_max, step=1., sigma_sky=None, sigma_cal=None, outside_frame=False)
    
    sma_azim_core,inten_azim_core,inten_azim_err_core = np.loadtxt('azim_model_core.txt', usecols=[0,1,2], dtype=float, unpack=True, skiprows = 1, delimiter='\t')
    #inten_azim_core = inten_azim_core/np.max(inten_azim_core)
    #inten_azim_err_core = inten_azim_err_core/np.max(inten_azim_core)
    

    sma_azim_ext,inten_azim_ext,inten_azim_err_ext = np.loadtxt(extended_psf_profile, usecols=[0,1,2], dtype=float, unpack=True, skiprows = 1, delimiter='\t')

    y_size, x_size = 2*np.max(sma_azim_ext),2*np.max(sma_azim_ext)
    image_x_center, image_y_center = x_size/2, y_size/2
    #object_x_center = xcen if xcen > 0 else image_x_center
    #object_y_center = ycen if ycen > 0 else image_y_center
    

    #inten_azim_ext = inten_azim_ext/np.max(inten_azim_ext)
    #inten_azim_err_ext = inten_azim_err_ext/np.max(inten_azim_ext)
    
    new_sma = []; new_inten = []
    R_max = np.max(sma_azim_core)
    I_max = inten_azim_core[list(sma_azim_core).index(R_max)]
    
    R = 0.
    k = 0
    while R<R_max and k<len(sma_azim_core):
        new_sma.append(sma_azim_core[k])
        new_inten.append(inten_azim_core[k]/I_max)
        k = k+1
    
    
    I_max = inten_azim_ext[list(sma_azim_ext).index(R_max)]
    for k in range(len(sma_azim_ext)):
        if sma_azim_ext[k]>=R_max+1:
            new_sma.append(sma_azim_ext[k])
            new_inten.append(inten_azim_ext[k]/I_max)            

    ellipses = IsophoteList([])
    #print(new_sma)
    #exit()
    for k in range(len(new_sma)):
            ellipse = SimpleNamespace(x0=image_x_center, y0=image_y_center, sma=new_sma[k], intens=new_inten[k],
                                  eps=0., pa=np.radians(0.), grad=0, a3=0, b3=0, a4=0, b4=0, fit_ok=False)
            ellipses.append(ellipse)
    
    model = azimProfile.build_ellipse_model((int(y_size), int(x_size)), ellipses)
    fits.PrimaryHDU(data=model).writeto(output_model, overwrite=True)    
    
    
main('psf_r.fits', 'azim_model.txt', 13, output_model='extended_psf.fits', core_fit=True)    
    
    
    
    
    
   