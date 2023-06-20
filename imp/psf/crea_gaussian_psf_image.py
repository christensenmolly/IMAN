#!/usr/bin/python

#import sys
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
import shutil
import subprocess
import random
from numpy import fft
import astropy.io.fits as pyfits
import re
import os

from photutils import datasets
from astropy.stats import sigma_clipped_stats


FNULL = open(os.devnull, 'w')
tmp_out = sys.stdout

import argparse


def normalize(input_file, output_file):
    hdulist1 = pyfits.open(input_file)
    data = hdulist1[0].data
    header = hdulist1[0].header
    new_data = data/np.sum(data)


    hdu = pyfits.PrimaryHDU(new_data, header)
    hdu.writeto(output_file, overwrite=True)


def main(function, fwhm_psf, ell, PA, beta, psf_size, output_fits='psf.fits'):
    nx_psf = psf_size
    ny_psf = psf_size
    if nx_psf%2==0:
        xc_psf = int(nx_psf/2. + 1)
    else:
        xc_psf = int(nx_psf/2. + 0.5)
    if ny_psf%2==0:
        yc_psf = int(ny_psf/2. + 1)
    else:
        yc_psf = int(ny_psf/2. + 0.5)


    

    f = open(r"modelPSF.txt", "w") 
    sys.stdout = f
    print("\n===============================================================================")
    print("# IMAGE and GALFIT CONTROL PARAMETERS")
    print("A) none                # Input data image (FITS file)")
    print("B) %s         # Output data image block" % (output_fits))
    print("C) none                # Sigma image name (made from data if blank or none)" )
    print("D) none                # Input PSF image and (optional) diffusion kernel")
    print("E) 1                   # PSF fine sampling factor relative to data" )
    print("F) none                # Bad pixel mask (FITS image or ASCII coord list)")
    print("G) none                # File with parameter constraints (ASCII file)" )
    print("H) 1    %i   1    %i   # Image region to fit (xmin xmax ymin ymax)" % (nx_psf, ny_psf))
    print("I) %.3f    %.3f        # Size of the convolution box (x y)" % (0, 0))
    print("J) %.3f              # Magnitude photometric zeropoint" % (25.))
    print("K) %.3f    %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (1.0,1.0))
    print("O) regular             # Display type (regular, curses, both)")
    print("P) 1                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n")

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

    if function=='gauss':
        print("# Gaussian function\n")
        print("0) gaussian           # object type")
        print("1) %.3f  %.3f  0 0   # position x, y        [pixel]" % (xc_psf,yc_psf))
        print("3) %.3f       0        # total magnitude" % (10.0)    ) 
        print("4) %.3f       0        #   FWHM               [pixels]" % (fwhm_psf))
        print("9) %.3f        0       # axis ratio (b/a)" % (1.-ell)  )
        print("10) %.3f         0       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA))
        print("Z) 0                  # leave in [1] or subtract [0] this comp from data")
        print("\n================================================================================")

    if function=='moffat':
        print("# Moffat function\n")
        print("0) moffat           # object type")
        print("1) %.3f  %.3f  0 0   # position x, y        [pixel]" % (xc_psf,yc_psf))
        print("3) %.3f       0        # total magnitude" % (10.0)     )
        print("4) %.3f       0        #   FWHM               [pixels]" % (fwhm_psf))
        print("5) %.3f        0       # powerlaw" % (beta))
        print("9) %.3f        0       # axis ratio (b/a)" % (1.-ell)  )
        print("10) %.3f         0       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA))
        print("Z) 0                  # leave in [1] or subtract [0] this comp from data")
        print("\n================================================================================"    )      

    sys.stdout = tmp_out
    f.close()

    subprocess.call("galfit modelPSF.txt", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    
    normalize(output_fits, output_fits)
    print('Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Sampling")  
    parser.add_argument("fwhm", help="FWHM [pix]", type=float, default=1.)
    parser.add_argument("--o", help="Optional: Output image",type=str,default='psf.fits')
    parser.add_argument("--w", help="Optional: FUnction: gauss or moffat", type=str, default='gauss')
    parser.add_argument("--ell", help="Optional: Ellipticity", type=float, default=0.)
    parser.add_argument("--PA", help="Optional: Position angle [Degrees: Up=0, Left=90]", type=float, default=90.)
    parser.add_argument("--beta", help="Optional: Beta (power) for the Moffat function", type=float, default=2.)
    parser.add_argument("--size", help="Optional: Size of the output psf image [nx=ny, pix]. Must be odd!", type=float, default=51.)
     
    args = parser.parse_args()


    output_image = args.o
    fwhm_psf = args.fwhm
    function = args.w
    ell = args.ell
    PA = args.PA
    beta = args.beta
    psf_size = args.size

    main(function, fwhm_psf, ell, PA, beta, psf_size, output_fits=output_image)
