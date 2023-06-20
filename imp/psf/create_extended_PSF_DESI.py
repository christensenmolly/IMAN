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
from photutils.isophote import IsophoteList
from photutils.isophote import build_ellipse_model
from types import SimpleNamespace


def main(core_psf_image, output_model, radius_arcsec=100., scale=0.262):
    azimProfile.main(core_psf_image, output_model=None, azim_tab='azim_model_core.txt', mask_image=None, xcen=-1, ycen=-1, ell=0., posang=0., sma_min=-1, sma_max=28, step=1., sigma_sky=None, sigma_cal=None, outside_frame=False, center_model=False, stop_where_negative=False, verbosity=True, linear=True)
    
    sma,inten,inten_err = np.loadtxt('azim_model_core.txt', usecols=[0,1,2], dtype=float, unpack=True, skiprows = 1, delimiter='\t')
    
    radius_max = int(radius_arcsec/scale)
    
    radius = np.arange(1., radius_max) # in pixels
    
    inten_wings = 0.00033 * radius**(-2)
    
    ind = list(inten).index(np.nanmin(inten))

    inten_wings = inten_wings * inten[ind]/inten_wings[list(radius).index(sma[ind])]
    
    extended_inten = [ ]
    for k in range(len(inten_wings)):
        try:
            extended_inten.append(inten[k])
        except:
            extended_inten.append(inten_wings[k])

    '''
    fig = plt.figure()
    ax = fig.add_subplot()
    plt.plot(sma, inten, 'o')
    plt.plot(radius, inten_wings, '*')
    plt.plot(radius, extended_inten, '-')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.show()
    '''
    
    ellipses = IsophoteList([])
    
    for k in range(len(extended_inten)):
        ellipse = SimpleNamespace(x0=radius_max, y0=radius_max, sma=radius[k], intens=extended_inten[k],
                                  eps=0., pa=0., grad=0, a3=0, b3=0, a4=0, b4=0, fit_ok=False)
        ellipses.append(ellipse)

    model = build_ellipse_model((int(radius_max*2.+1),int(radius_max*2.+1)), ellipses)
    pyfits.PrimaryHDU(model).writeto(output_model, overwrite=True)


main('psf_r.fits', 'psf_extended.fits')

def combine_two_profiles(core_ellipse, wings_ellipse, R_intersect, output_image, rmax, m0, pix2sec):

    f = open('psf_extended_for_paper.txt', 'w')
    f.write('sma[pix]\tinten[DN]\n')
    for k in range(len(sma_wings)):
       f.write('%.3f\t%.8f\n' % (sma[k],inten[k]))
    f.close()


def main(core_psf, wings_psf, R_intersect, core_fit=True):
    # Read in the input PSF 
    hdulist_core = pyfits.open(core_psf)
    data_core = hdulist_core[0].data
    nx,ny = data_core.shape[1], data_core.shape[0]
    
    # Determine the center:
    if nx%2==False:
      xc = nx/2. + 0.5
    else:
      xc = nx/2. + 1.

    if ny%2==False:
      yc = ny/2. + 0.5
    else:
      yc = ny/2. + 1.
      

    iraf_ellipse.main_ell(core_psf, 14., 14., ellip=0.05, pa=0., sma=5., m0=20., pix2sec=1.0, step=1.0, minsma=0., maxsma=21., outp_format='png', ell_file='ellipse_star.txt', fits_mask=None, fflag=1.0, olthresh=0.0, linear='yes')
    shutil.move('azim_aver.png','azim_aver_star.png')


    if core_fit:

        # Fit the core_psf with a Moffat function

        # Create input galfit file:
        f = open('galfit.inp', "w") 
        sys.stdout = f  
        header(core_psf, 'none', None, 'none', 20., generalScaleValue=1., sampling=1)
        moffat_function(xc,yc,20. - 2.5*log10(np.sum(data_core)),4.,3.,1.0,0.)
        sys.stdout = tmp_out
        f.close()

        subprocess.call("galfit %s" % ("galfit.inp"), shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)

        xc,yc,mtot,FWHM,beta,q,PA = read_galfit_file('galfit.01')
    
        # Create input galfit file:
        f = open('galfit_core_model.inp', "w") 
        sys.stdout = f  
        header_model()
        moffat_function(251.,251.,mtot,FWHM,beta,q,PA)
        sys.stdout = tmp_out
        f.close()

        subprocess.call("galfit %s" % ("galfit_core_model.inp"), shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)
        
        iraf_ellipse.main_ell('psf_core_model.fits', 251., 251., ellip=0.05, pa=0., sma=100., m0=20., pix2sec=1.0, step=1.0, minsma=0., maxsma=300., outp_format='png', ell_file='ellipse_core.txt', fits_mask=None, fflag=1.0, olthresh=0.0, linear='yes')
        shutil.move('azim_aver.png','azim_aver_core.png')
        #exit()

