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

def plot_azim(sma, inten, m0, pix2sec, rmax=None, I_norm=None):
        
        
        for kk in range(len(sma)):
            yy = np.array(sma[0:kk+1])*np.array(inten[0:kk+1])
            L = 2.*math.pi*np.trapz(yy, dx = 1.)*(1.-0.)
        inten = inten/L
        fig = figure(1, figsize=(6,5))
        ax = fig.add_subplot(111)
        mag = np.log10(inten)#+ 5.*np.log10(pix2sec)
        fsize=15

        color = 'black'
        symb = 'o'
        if rmax is None:
            rmax = np.max(sma)
        
        
        
        

        ax.plot(sma*pix2sec, mag, '-',color=color,markeredgecolor=color, markersize=2, lw=3)
        #ax.plot(sma*pix2sec, mag, symb,color=color,markeredgecolor=color, markersize=2)
        mean_mag, median_mag, std_mag = sigma_clipped_stats(mag[~np.isnan(mag)], sigma=3.0, iters=5)
        
        if rmax is None:
            min_mag = min(mag[~np.isnan(mag)])-std_mag
        else:
            min_mag = min(mag[np.where(sma<=rmax)[0]])-std_mag
            
        #ax.set_ylim(max(mag[~np.isnan(mag)])+std_mag,min(mag[~np.isnan(mag)])-std_mag)
        ax.set_ylim(min_mag,max(mag[~np.isnan(mag)])+std_mag)
        ax.set_xlim(0., rmax*pix2sec)

        for tick in ax.yaxis.get_major_ticks():
              tick.label.set_fontsize(fsize-2)   

        for tick in ax.xaxis.get_major_ticks():
              tick.label.set_fontsize(fsize-2)  

        ax.set_xlabel(r'$r$ (arcsec)', fontsize=fsize)
        ax.set_ylabel(r'Normalized Flux', fontsize=fsize)
 
        plt.savefig("PSF_extended.png", transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
        plt.clf()
        plt.close()

    

def read_galfit_file(galfit_file):
    # Function to parse galfit file
    fff = open(galfit_file, 'r')
    lines = fff.readlines()

    for line in lines:
        #print line
        if '1)' in line and '#  Position x, y' in line:
            xc = float(line.split()[1])
            yc = float(line.split()[2])
        if '3)' in line and '#  Integrated magnitude' in line:
            mtot = float(line.split()[1])
        if '4)' in line and '#     FWHM       [pix]' in line:
            FWHM = float(line.split()[1])
        if '5)' in line and '#  Moffat powerlaw index' in line:
            beta = float(line.split()[1])
        if '9)' in line and '#  Axis ratio (b/a)' in line:
            q = float(line.split()[1])
        if '10)' in line and '#  Position angle (PA) [deg: Up=0, Left=90]' in line:
            PA = float(line.split()[1])            
    
    return xc,yc,mtot,FWHM,beta,q,PA

def header(file_Image, file_Sigma, file_Psf, file_Mask, magZP, generalScaleValue=1., sampling=1):
        # Function to create header for galfit file
        hdulist = pyfits.open(file_Image)
        image = hdulist[0].data    
        ySize, xSize = image.shape
        if file_Psf is not None:
            hdulist = pyfits.open(file_Psf)
            image = hdulist[0].data    
            ySize_psf, xSize_psf = image.shape
        
            #psf_box = min([min([ySize_psf, xSize_psf]),ySize])
            psf_box = min([ySize_psf, xSize_psf])
        else:
            psf_box = 50
        print("\n===============================================================================")
        print("# IMAGE and GALFIT CONTROL PARAMETERS")
        print("A) %s                  # Input data image (FITS file)" % (file_Image))
        print("B) %s                  # Output data image block" % ('model.fits'))
        print("C) %s                # Sigma image name (made from data if blank or none)" % (str(file_Sigma)))
        print("D) %s                # Input PSF image and (optional) diffusion kernel" % (str(file_Psf)))
        print("E) %i                   # PSF fine sampling factor relative to data" % (sampling))
        print("F) %s                # Bad pixel mask (FITS image or ASCII coord list)" % (str(file_Mask)))
        print("G) %s                # File with parameter constraints (ASCII file)" % ('none'))
        print("H) %i    %i   %i    %i   # Image region to fit (xmin xmax ymin ymax)" % (1, xSize, 1, ySize))
        print("I) %i    %i              # Size of the convolution box (x y)" % (psf_box,psf_box))
        print("J) %.3f              # Magnitude photometric zeropoint" % (magZP))
        print("K) %.3f  %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (generalScaleValue,generalScaleValue))
        print("O) regular             # Display type (regular, curses, both)")
        print("P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n")

        print("# INITIAL FITTING PARAMETERS")
        print("#")
        print("#   For object type, the allowed functions are:")
        print("#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,")
        print("#       ferrer, powsersic, sky, and isophote.") 
        print("#")  
        print("#   Hidden parameters will only appear when they're specified:")
        print("#       C0 (diskyness/boxyness),") 
        print("#       Fn (n=integer, Azimuthal Fourier Modes),")
        print("#       R0-R10 (PA rotation, for creating spiral structures).")
        print("#") 
        print("# -----------------------------------------------------------------------------")
        print("#   par)    par value(s)    fit toggle(s)    # parameter description") 
        print("# -----------------------------------------------------------------------------\n")


def header_model():
        print("\n===============================================================================")
        print("# IMAGE and GALFIT CONTROL PARAMETERS")
        print("A) none                  # Input data image (FITS file)")
        print("B) psf_core_model.fits                  # Output data image block")
        print("C) none                # Sigma image name (made from data if blank or none)")
        print("D) none                # Input PSF image and (optional) diffusion kernel")
        print("E) 1                   # PSF fine sampling factor relative to data")
        print("F) none                # Bad pixel mask (FITS image or ASCII coord list)")
        print("G) none                # File with parameter constraints (ASCII file)")
        print("H) 1    501   1    501   # Image region to fit (xmin xmax ymin ymax)")
        print("I) 50    50              # Size of the convolution box (x y)")
        print("J) 20.0              # Magnitude photometric zeropoint")
        print("K) 1.000  1.000        # Plate scale (dx dy)    [arcsec per pixel]")
        print("O) regular             # Display type (regular, curses, both)")
        print("P) 1                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n")

        print("# INITIAL FITTING PARAMETERS")
        print("#")
        print("#   For object type, the allowed functions are:")
        print("#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,")
        print("#       ferrer, powsersic, sky, and isophote.") 
        print("#")  
        print("#   Hidden parameters will only appear when they're specified:")
        print("#       C0 (diskyness/boxyness),") 
        print("#       Fn (n=integer, Azimuthal Fourier Modes),")
        print("#       R0-R10 (PA rotation, for creating spiral structures).")
        print("#") 
        print("# -----------------------------------------------------------------------------")
        print("#   par)    par value(s)    fit toggle(s)    # parameter description") 
        print("# -----------------------------------------------------------------------------\n")

def moffat_function(xc,yc,mtot,FWHM,beta,q,PA):
    s = '''
 0) moffat             # object type
 1) %f  %f 1 1   # position x, y        [pixel]
 3) %f       1       # total magnitude     
 4) %f        1       #   FWHM               [Pixels]
 5) %f        1       # powerlaw      
 9) %f        0       # axis ratio (b/a)   
10) %f         0       # position angle (PA)  [Degrees: Up=0, Left=90]
 Z) 0                  #  Skip this model in output image?  (yes=1, no=0)
 ''' % (xc,yc,mtot,FWHM,beta,q,PA)
    print(s)


def read_ellipse_file(iraf_ell_file):
    sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = np.loadtxt(iraf_ell_file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 5, dtype='str')

    for k in range(len(sma)):
                    if sma[k]=='INDEF': sma[k]=0
                    if inten[k]=='INDEF': inten[k]=0
                    if inten_err[k]=='INDEF': inten_err[k]=0
                    if ell[k]=='INDEF': ell[k]=0
                    if errell[k]=='INDEF': errell[k]=0
                    if PA[k]=='INDEF': PA[k]=0
                    if errPA[k]=='INDEF': errPA[k]=0
                    if x0[k]=='INDEF': x0[k]=0
                    if y0[k]=='INDEF': y0[k]=0
                    if B4[k]=='INDEF': B4[k]=0
                    if errB4[k]=='INDEF': errB4[k]=0
    sma = np.array(sma,dtype='float')
    inten = np.array(inten,dtype='float')
    inten_err = np.array(inten_err,dtype='float')
    ell = np.array(ell,dtype='float')
    errell = np.array(errell,dtype='float')
    PA = np.array(PA,dtype='float')
    errPA = np.array(errPA,dtype='float')
    x0 = np.array(x0,dtype='float')
    y0 = np.array(y0,dtype='float')
    B4 = np.array(B4,dtype='float')
    errB4 = np.array(errB4,dtype='float')
    return sma, inten, inten_err, ell, errell, PA, errPA, x0, y0, B4, errB4


def combine_two_profiles(core_ellipse, wings_ellipse, R_intersect, output_image, rmax, m0, pix2sec):
    sma_core, inten_core, inten_err_core, ell_core, errell_core, PA_core, errPA_core, x0_core, y0_core, B4_core, errB4_core = read_ellipse_file(core_ellipse)
    
    sma_wings, inten_wings, inten_err_wings, ell_wings, errell_wings, PA_wings, errPA_wings, x0_wings, y0_wings, B4_wings, errB4_wings = read_ellipse_file(wings_ellipse)
    
    #inten_wings = savgol_filter(inten_wings, 21, 3)   #### CHECK THIS

    sma = []
    inten = []
    inten_err = []
    ell = []
    errell = []
    PA = []
    errPA = []
    x0 = []
    y0 = []
    B4 = []
    errB4 = []
    
    for k in range(len(sma_core)):
        if sma_core[k]<=R_intersect:
            sma.append(sma_core[k])
            inten.append(inten_core[k])
            inten_err.append(inten_err_core[k])
            ell.append(ell_core[k])
            errell.append(errell_core[k])
            PA.append(PA_core[k])
            errPA.append(errPA_core[k])
            x0.append(x0_core[k])
            y0.append(y0_core[k])
            B4.append(B4_core[k])
            errB4.append(errB4_core[k])
        else:
            break
    inten = list(np.array(inten, float)/inten[-1])
    #print(inten)
    #exit()
    
    for k in range(len(sma_wings)):
        if sma_wings[k]<=R_intersect:
            inten_norm = inten_wings[k]
            
        if sma_wings[k]>R_intersect:
            sma.append(sma_wings[k])
            inten.append(inten_wings[k]/inten_norm)
            inten_err.append(inten_err_wings[k])
            ell.append(ell_wings[k])
            errell.append(errell_wings[k])
            PA.append(PA_wings[k])
            errPA.append(errPA_wings[k])
            x0.append(x0_wings[k])
            y0.append(y0_wings[k])
            B4.append(B4_wings[k])
            errB4.append(errB4_wings[k])
    sma = np.array(sma,dtype='float')
    inten = np.array(inten,dtype='float')
    ell = np.array(ell,dtype='float')
    PA = np.array(PA,dtype='float')

    #plt.plot(sma, 28. -2.5*np.log10(inten), 'o')
    #plt.show()
    #exit()

    #'''
    f = open('psf_extended_for_paper.txt', 'w')
    f.write('sma[pix]\tinten[DN]\n')
    for k in range(len(sma_wings)):
       f.write('%.3f\t%.8f\n' % (sma[k],inten[k]))
    f.close()
    #exit()
    #'''

    plot_azim(sma, inten, m0, pix2sec, rmax)
    #exit()
    build_iraf_model.build_model(sma, inten, ell, PA, output_image, input_image=None, rmax=rmax, center=True, fix_ellipses=True)    
    norm_image.main(output_image, output_image)

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

#main('psf.fits', 'psf1.fits', 30., core_fit=False)

hdulist_core = pyfits.open('psf.fits')
data_core = hdulist_core[0].data
nx,ny = data_core.shape[1], data_core.shape[0]

#
#exit()
xc = 1203.81
yc = 399.61
m0 = 28.8607
pix2sec = 1.67
arithm_operations.main('new-image.fits', 4.577+2.32, 'add', 'new-image_sky_sub.fits')
#iraf_ellipse.main_ell('new-image_sky_sub.fits', xc, yc, ellip=0.05, pa=0., sma=100., m0=20., pix2sec=1.0, step=1.0, minsma=0., maxsma=700., outp_format='png', ell_file='ellipse_wings.txt', fits_mask='mask.fits', fflag=1.0, olthresh=0.0, linear='yes')
#shutil.move('azim_aver.png','azim_aver_wings.png')
#combine_two_profiles('ellipse_core.txt', 'ellipse_wings.txt', 20., 'final_ext_psf.fits', 600., m0, pix2sec)
combine_two_profiles('ellipse_star.txt', 'ellipse_wings.txt', 15., 'final_ext_psf_star.fits', 635., m0, pix2sec)




#azimProfile.main('new-image_sky_sub.fits', output_model='axim_model.fits', azim_tab='azim_model.txt', mask_image='mask.fits', xcen=xc, ycen=yc, ell=0., posang=0., sma_min=-1, sma_max=10., step=1., sigma_sky=None, sigma_cal=None, outside_frame=True)

#compare_iraf_and_ellipse.main('ellipse.txt', 'azim_model.txt', pix2sec=1.0, m0=28.0)

#build_iraf_model.main(input_image=None, iraf_ell_file='ellipse.txt', output_model='ellipse.fits', center=True, rmax=600., fix_ellipses=True)




