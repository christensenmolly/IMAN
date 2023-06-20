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

FNULL = open(os.devnull, 'w')
tmp_out = sys.stdout

def Header(file_image, file_out, bad_pix_mask, m0, pix2sec, nx, ny):
    print("\n===============================================================================")
    print("# IMAGE and GALFIT CONTROL PARAMETERS")
    print("A) %s                # Input data image (FITS file)" % (file_image))
    print("B) %s                    # Output data image block" % (file_out))
    print("C) none                # Sigma image name (made from data if blank or none)" )
    print("D) none                # Input PSF image and (optional) diffusion kernel")
    print("E) 1                   # PSF fine sampling factor relative to data" )
    print("F) %s                # Bad pixel mask (FITS image or ASCII coord list)" % (bad_pix_mask))
    print("G) none                # File with parameter constraints (ASCII file)" )
    print("H) %i    %i   %i    %i   # Image region to fit (xmin xmax ymin ymax)" % (1,nx,1,ny))
    print("I) 0    0        # Size of the convolution box (x y)")
    print("J) %.3f              # Magnitude photometric zeropoint" % (m0))
    print("K) %.3f    %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (pix2sec,pix2sec))
    print("O) regular             # Display type (regular, curses, both)")
    print("P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n")

    print("# INITIAL FITTING PARAMETERS")
    print("#")
    print("#   For object type, the allowed functions are:" )
    print("#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat," )
    print("#       ferrer, powsersic, sky, and isophote." )
    print("#"  )
    print("#   Hidden parameters will only appear when they're specified:")
    print("#       C0 (diskyness/boxyness)," )
    print("#       Fn (n=integer, Azimuthal Fourier Modes),")
    print("#       R0-R10 (PA rotation, for creating spiral structures).")
    print("#" )
    print("# -----------------------------------------------------------------------------")
    print("#   par)    par value(s)    fit toggle(s)    # parameter description" )
    print("# -----------------------------------------------------------------------------\n")

def moffat(xc, yc, mtot, fwhm, ell, PA):
     print("Moffat function")
     print("0) moffat           # object type")
     print("1) %.3f  %.3f  1 1  # position x, y        [pixel]" % (xc,yc))
     print("3) %.3f       1       # total magnitude" % (mtot)     )
     print("4) %.3f       1       #   FWHM               [pixels]" % (fwhm))
     print("5) %.3f        1       # powerlaw" % (beta) )
     print("9) %.3f      1       # axis ratio (b/a)" % (1.-ell))
     print("10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA))
     print("Z) 0                  # leave in [1] or subtract [0] this comp from data?")

def gauss(xc, yc, mtot, fwhm, ell, PA):
     print("Gaussian function")
     print("0) gaussian           # object type")
     print("1) %.3f  %.3f  1 1  # position x, y        [pixel]" % (xc,yc))
     print("3) %.3f       1       # total magnitude" % (mtot)     )
     print("4) %.3f       1       #   FWHM               [pixels]" % (fwhm))
     print("9) %.3f      1       # axis ratio (b/a)" % (1.-ell))
     print("10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA))
     print("Z) 0                  # leave in [1] or subtract [0] this comp from data?")


def read_galfit_file(galfit_file):
    beta = float('nan')
    mag = float('nan') 
    FWHM = float('nan') 
    q = float('nan') 
    PA = float('nan') 
    
    if os.path.exists(galfit_file):
        fff = open(galfit_file, 'r')
        lines = fff.readlines()
        
        for line in lines:
            if '3)' in line and '#  Integrated magnitude' in line:
                mag = float(line.split()[1])    
            if '4)' in line and '#     FWHM       [pix]' in line:
                FWHM = float(line.split()[1])
            if '5)' in line and '#  powerlaw' in line:
                beta = float(line.split()[1])
            if '9)' in line and '#  Axis ratio (b/a)' in line:
                q = float(line.split()[1])    
            if '10)' in line and '#  Position angle (PA) [deg: Up=0, Left=90]' in line:
                PA = float(line.split()[1])    
        fff.close()
    return mag,FWHM,beta,q,PA
    

def remove_files(files):
    for file in files:
        if os.path.exists(file):
            os.remove(file)


def main(input_image, window='gauss', remove_model=False):
        hdulist_core = pyfits.open(input_image)
        data_core = hdulist_core[0].data
        ny,nx = np.shape(data_core)

        pix2sec = 1.
        m0 = 20.
        
        f = open(r"modelIN.txt", "w") 
        sys.stdout = f
        Header(input_image, 'psf_model_galf.fits', 'none', m0, pix2sec, nx, ny)
        if window=='moffat':    moffat(nx/2.,ny/2.,m0-2.5*math.log10(np.sum(data_core)),nx/6.,0.0,0.)
        if window=='gauss':    gauss(nx/2.,ny/2.,m0-2.5*math.log10(np.sum(data_core)),nx/6.,0.0,0.)
        sys.stdout = tmp_out
        f.close()
        #os.chmod(r"modelIN.txt",0777)

        ## GALFIT RUNING
        subprocess.call("galfit modelIN.txt", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
        mag,FWHM,beta,q,PA = read_galfit_file('galfit.01')
        os.remove('galfit.01')
        os.remove('fit.log')
        remove_files(['galfit.01','fit.log','modelIN.txt'])
        if remove_model:
            remove_files(['psf_model_galf.fits'])
        
        return mag,FWHM,beta,q,PA
    

if __name__ == "__main__":
    input_image = sys.argv[1]
    window = sys.argv[2]
    main(input_image, window=window)
    
        
    