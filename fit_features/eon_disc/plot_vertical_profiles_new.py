#! /usr/bin/env python
import numpy as np
import gzip
import shutil
from joblib import Parallel, delayed
import astropy.io.fits as pyfits
import sys
import os
import shutil
import time
import subprocess
import glob
import math
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import warnings
import random
from scipy import signal
from scipy.special import kn
from scipy import stats
from math import *
warnings.filterwarnings("ignore")
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

import fit_eon_vertical
FNULL = open(os.devnull, 'w')
path = ''#'/Users/mosenkov/MEGA/HERON/ESO240-G011/'
imfit_path = path#'/home/amosenko/MEGA/MyPrograms/imfit-1.6.1/'
galfit_path = path
fitscopy = path

fsize = 18


def ds9_to_np(x):
    '''
    Function to convert ds9 pixel coordinates to numpy format, i.e
    0.5 pix is 0, 1.5 pix is 1
    '''
    return int(round(x)) - 1

def run_gulf_single_good(psf_image, r, xc, yc, nx, ny, mu0d, z0, h, M0):
    if psf_image is None:
        psf_image = 'none'
        n_psf=1
    else:
        hdu = pyfits.open(psf_image)
        image_psf = hdu[0].data 
        ny_psf,nx_psf =  np.shape(image_psf)
        n_psf = max([ny_psf,nx_psf,ny])
        
    s1="""

================================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) none         # Input data image (FITS file)
B) model_slice.fits          # Output data image block
C) none          # Sigma image name (made from data if blank or "none") 
D) %s            # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data 
F) none           # Bad pixel mask (FITS image or ASCII coord list)
G) none                # File with parameter constraints (ASCII file) 
H) 1  %i  1  %i  # Image region to fit (xmin xmax ymin ymax)
I) %i    %i          # Size of the convolution box (x y)
J) %f              # Magnitude photometric zeropoint 
K) 1.0  1.0        # Plate scale (dx dy)   [arcsec per pixel]
O) regular             # Display type (regular, curses, both)
P) 1                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

# INITIAL FITTING PARAMETERS
#
#   For component type, the allowed functions are: 
#       sersic, expdisk, edgedisk, devauc, king, nuker, psf, 
#       gaussian, moffat, ferrer, and sky. 
#  
#   Hidden parameters will only appear when they're specified:
#       Bn (n=integer, Bending Modes).
#       C0 (diskyness/boxyness), 
#       Fn (n=integer, Azimuthal Fourier Modes).
#       R0-R10 (coordinate rotation, for creating spiral structures).
#       To, Ti, T0-T10 (truncation function).
# 
# ------------------------------------------------------------------------------
#   par)    par value(s)    fit toggle(s)    # parameter description 
# ------------------------------------------------------------------------------
""" % (psf_image, nx, ny, nx, ny, M0)

    s2="""
# Component number: 1
 0) edgedisk               #  Component type
 1) %f %f %i %i  #  Position x, y
 3) %f     1          #     Mu(0)   [mag/arcsec^2]
 4) %f      %i          #  h_s (disk scale-height)   [pix]
 5) %f      %i          #  R_s (disk scale-length)   [pix]
10) 90.0     0          #  Position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)
 ================================================================================
""" % (xc, yc, 0, 0, mu0d, z0, 0, h, 0)
        
    s = s1 + s2
    
    f = open("galf_slice.txt", 'w')
    f.write(s)
    f.close()
    
    subprocess.call('%sgalfit %s' % (galfit_path, 'galf_slice.txt'), shell=True, stdout=FNULL)
    os.remove('galf_slice.txt')
    
    imageHDU = pyfits.open('model_slice.fits')[0]
    image = imageHDU.data
    os.remove('model_slice.fits')
    #exit()
    return np.array(image[:,ds9_to_np(xc+r)]) # TODO: Averaging as in fit_eon_vertical!!!!!!



def run_gulf_single(psf_image, yc, ny, mu0d, z0, h, M0):
    if psf_image is None:
        psf_image = 'none'
        n_psf=1
    else:
        hdu = pyfits.open(psf_image)
        image_psf = hdu[0].data 
        ny_psf,nx_psf =  np.shape(image_psf)
        n_psf = max([ny_psf,nx_psf,ny])

    s1="""

================================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) none         # Input data image (FITS file)
B) model_line.fits          # Output data image block
C) none          # Sigma image name (made from data if blank or "none") 
D) %s            # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data 
F) none           # Bad pixel mask (FITS image or ASCII coord list)
G) none                # File with parameter constraints (ASCII file) 
H) 1  %i  1  %i  # Image region to fit (xmin xmax ymin ymax)
I) %i    %i          # Size of the convolution box (x y)
J) %f              # Magnitude photometric zeropoint 
K) 1.0  1.0        # Plate scale (dx dy)   [arcsec per pixel]
O) regular             # Display type (regular, curses, both)
P) 1                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

# INITIAL FITTING PARAMETERS
#
#   For component type, the allowed functions are: 
#       sersic, expdisk, edgedisk, devauc, king, nuker, psf, 
#       gaussian, moffat, ferrer, and sky. 
#  
#   Hidden parameters will only appear when they're specified:
#       Bn (n=integer, Bending Modes).
#       C0 (diskyness/boxyness), 
#       Fn (n=integer, Azimuthal Fourier Modes).
#       R0-R10 (coordinate rotation, for creating spiral structures).
#       To, Ti, T0-T10 (truncation function).
# 
# ------------------------------------------------------------------------------
#   par)    par value(s)    fit toggle(s)    # parameter description 
# ------------------------------------------------------------------------------
""" % (psf_image, 1, ny, n_psf, n_psf, M0)

    s2="""
# Component number: 1
 0) edgedisk               #  Component type
 1) %f %f %i %i  #  Position x, y
 3) %f     1          #     Mu(0)   [mag/arcsec^2]
 4) %f      %i          #  h_s (disk scale-height)   [pix]
 5) %f      %i          #  R_s (disk scale-length)   [pix]
10) 90.0     0          #  Position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)
 ================================================================================
""" % (0.5, yc, 0, 0, mu0d, z0, 0, h, 0)
        
    s = s1 + s2
    
    f = open("galf_line.txt", 'w')
    f.write(s)
    f.close()
    
    subprocess.call('%sgalfit %s' % (galfit_path, 'galf_line.txt'), shell=True, stdout=FNULL)
    os.remove('galf_line.txt')
    
    imageHDU = pyfits.open('model_line.fits')[0]
    image = imageHDU.data
    os.remove('model_line.fits')
    return np.array(image[:,0])



def run_gulf_double(psf_image, yc, ny, mu0d1, z01, h1, mu0d2, z02, h2, M0):
    if psf_image is None:
        psf_image = 'none'
        n_psf=1
    else:
        hdu = pyfits.open(psf_image)
        image_psf = hdu[0].data 
        ny_psf,nx_psf =  np.shape(image_psf)
        n_psf = max([ny_psf,nx_psf,ny])

    s1="""

================================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) none         # Input data image (FITS file)
B) model_line.fits          # Output data image block
C) none          # Sigma image name (made from data if blank or "none") 
D) %s            # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data 
F) none           # Bad pixel mask (FITS image or ASCII coord list)
G) none                # File with parameter constraints (ASCII file) 
H) 1  %i  1  %i  # Image region to fit (xmin xmax ymin ymax)
I) %i    %i          # Size of the convolution box (x y)
J) %f              # Magnitude photometric zeropoint 
K) 1.0  1.0        # Plate scale (dx dy)   [arcsec per pixel]
O) regular             # Display type (regular, curses, both)
P) 1                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

# INITIAL FITTING PARAMETERS
#
#   For component type, the allowed functions are: 
#       sersic, expdisk, edgedisk, devauc, king, nuker, psf, 
#       gaussian, moffat, ferrer, and sky. 
#  
#   Hidden parameters will only appear when they're specified:
#       Bn (n=integer, Bending Modes).
#       C0 (diskyness/boxyness), 
#       Fn (n=integer, Azimuthal Fourier Modes).
#       R0-R10 (coordinate rotation, for creating spiral structures).
#       To, Ti, T0-T10 (truncation function).
# 
# ------------------------------------------------------------------------------
#   par)    par value(s)    fit toggle(s)    # parameter description 
# ------------------------------------------------------------------------------
""" % (psf_image, 1, ny, n_psf, n_psf, M0)

    s2="""
# Component number: 1
 0) edgedisk               #  Component type
 1) %f %f %i %i  #  Position x, y
 3) %f     1          #     Mu(0)   [mag/arcsec^2]
 4) %f      %i          #  h_s (disk scale-height)   [pix]
 5) %f      %i          #  R_s (disk scale-length)   [pix]
10) 90.0     0          #  Position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)

# Component number: 2
 0) edgedisk               #  Component type
 1) %f %f %i %i  #  Position x, y
 3) %f     1          #     Mu(0)   [mag/arcsec^2]
 4) %f      %i          #  h_s (disk scale-height)   [pix]
 5) %f      %i          #  R_s (disk scale-length)   [pix]
10) 90.0     0          #  Position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)
 ================================================================================
""" % (0.5, yc, 0, 0, mu0d1, z01, 0, h1, 0, 0.5, yc, 0, 0, mu0d2, z02, 0, h2, 0)
        
    s = s1 + s2
    
    f = open("galf_line.txt", 'w')
    f.write(s)
    f.close()
    
    subprocess.call('%sgalfit %s' % (galfit_path, 'galf_line.txt'), shell=True, stdout=FNULL)
    os.remove('galf_line.txt')
    
    imageHDU = pyfits.open('model_line.fits')[0]
    image = imageHDU.data
    os.remove('model_line.fits')
    return np.array(image[:,0])






def disk_edge_soph(z, I0d, z0, z_c, n):
    #*** For edge-on disk SB in mag/arcsec^2 (along z-axis). n - index of this law. ***
    return I0d * (1.0 / np.cosh(n*np.fabs(z-z_c)/(2.*z0)) ) **(2./n)

def read_profiles(profiles_file, mode):
    f=open(profiles_file, 'r')
    lines = f.read()
    cuts = lines.split('#END\n')

    if mode=='single':
        radius = []; Z = []; I = []; I_MOD = []
        for cut in cuts:
            if cut!='':
                cut = cut.split('\n')
                #print(cut[0],'here')
                radius.append(float( cut[0].split('#RADIUS:')[-1] ))
                zz = []; II = []; I_mod = []
                #print(cut)
                for k in range(1, len(cut)):
                    if cut[k]!='':
                        #print(cut[k])
                        z,i,i_mod = cut[k].split('\t')
                        zz.append(float(z))
                        II.append(float(i))
                        I_mod.append(float(i_mod))
                Z.append(zz)
                I.append(II)
                I_MOD.append(I_mod)
        return np.array(radius),np.array(Z),np.array(I),np.array(I_MOD)

    if mode=='double1' or mode=='double':
        radius = []; Z = []; I = []; I_MOD = []; I_1 = []; I_2 = []
        for cut in cuts:
            if cut!='':
                cut = cut.split('\n')
                #print(cut[0],'here')
                radius.append(float( cut[0].split('#RADIUS:')[-1] ))
                zz = []; II = []; I_mod = []; I_11 = []; I_22 = []
                #print(cut)
                for k in range(1, len(cut)):
                    if cut[k]!='':
                        #print(cut[k])
                        z,i,i_mod,i_1,i_2 = cut[k].split('\t')
                        zz.append(float(z))
                        II.append(float(i))
                        I_mod.append(float(i_mod))
                        I_11.append(float(i_1))
                        I_22.append(float(i_2))
                Z.append(zz)
                I.append(II)
                I_MOD.append(I_mod)
                I_1.append(I_11)
                I_2.append(I_22)
        return np.array(radius),np.array(Z),np.array(I),np.array(I_MOD),np.array(I_1),np.array(I_2)    

def main(profiles_file, units = ['arcsec',0.75], mode='single'):
    
    if mode=='single':
        radius,Z,I,I_mod = read_profiles(profiles_file, mode)
    if mode=='double1':
        radius,Z,I,I_mod,I_1,I_2 = read_profiles(profiles_file, mode)
    if mode=='double':
        radius,Z,I,I_mod,I_1,I_2 = read_profiles(profiles_file, mode)

    
    
    
    
    #R,L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err,L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err,L03,L03_err,N3,N3_err,Z03,Z03_err = np.loadtxt(results_file, usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26],dtype=float, unpack=True,skiprows=1,delimiter='\t')

    fig = plt.figure(0, figsize=(5, 10))
    ax=plt.subplot()
    
    M0 = 28.
    min_mag = []
    max_mag = []
    max_z = []
    for k in range(len(radius)):
        mag = M0 - np.log10(I[k]) 
        
        ax.plot(Z[k], mag,'o', markersize=9,color='gray',markeredgecolor='gray')

        if mode=='single':
            mag_single = M0 - np.log10(I_mod[k]) 
            ax.plot(Z[k], mag_single, '-', color='red')
        
        if mode=='double1' or mode=='double':           
            mag_1 = M0 - np.log10(I_1[k])
            mag_2 = M0 - np.log10(I_2[k])     
            mag_sum = M0 - np.log10(I_mod[k])     
            ax.plot(Z[k], mag_sum, '-', color='black')
            ax.plot(Z[k], mag_1, '--', color='blue')            
            ax.plot(Z[k], mag_2, '--', color='red')
            if np.isnan(mag_2[0]):
                ax.plot(Z[k], mag_sum, '--', color='red')

        if mode=='tripple':
            # TODO:
                
            mag_1 = M0 - np.log10(I_1)
            mag_2 = M0 - np.log10(I_2)     
            mag_3 = M0 - np.log10(I_3)  
            mag_sum = M0 - np.log10(I_sum) 
            ax.plot(Z[k], mag_sum, '-', color='black')
            ax.plot(Z[k], mag_1, '--', color='blue')            
            ax.plot(Z[k], mag_2, '--', color='red')
            ax.plot(Z[k], mag_3, '--', color='lime')
            
        #min_mag.append(np.nanmin(mag))
        #max_mag.append(np.nanmax(mag))
        mmag = []
        for i in range(len(mag)):
            if not np.isinf(mag[i]) and not np.isnan(mag[i]):
                mmag.append(mag[i])     
        
        
        min_mag.append(np.min(mmag))
        max_mag.append(np.max(mmag))  
        max_z.append(np.nanmax(np.fabs(Z[k])))
        M0 = M0 - 2
    
    Min_mag = []
    Max_mag = []
    for k in range(len(min_mag)):
        if not np.isinf(min_mag[k]):
            Min_mag.append(min_mag[k])

    for k in range(len(max_mag)):
        if not np.isinf(max_mag[k]):
            Max_mag.append(max_mag[k])
            
    #print(np.nanmax(Max_mag)+1, np.nanmin(Min_mag)-1)
    ##exit()
    #print(max_mag, Min_mag)
    ax.set_xlim(0., np.nanmax(max_z))
    ax.set_ylim(np.nanmax(Max_mag)+1, np.nanmin(Min_mag)-1)
    
    ax.set_xlabel(r' $z$ (%s)' % (units[0]), fontsize=fsize)
    ax.set_ylabel(r' $\mu$ (mag\,arcsec$^{-2}$)', fontsize=fsize)
    
    plt.savefig('z_porfiles_%s.eps' % (mode), transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.clf()
    plt.close()

# Test_eon/one_disc_galaxy
#main('vertical_profiles_1_single.txt', units = ['pix',1.0], mode='single')
#main('vertical_profiles_1_double.txt', units = ['pix',1.0], mode='double')
