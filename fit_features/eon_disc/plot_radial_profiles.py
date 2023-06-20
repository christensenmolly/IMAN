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
from scipy.special import kn
from scipy import signal

warnings.filterwarnings("ignore")
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


import fit_eon_radial

fsize = 18


def disk_edge_soph(z, r, I0d, z0, h, x_c, z_c, n):
    #*** For edge-on disk SB in mag/arcsec^2 (along z-axis). n - index of this law. ***
    return I0d * (np.fabs(r-x_c)/h)*kn(1,np.fabs(r-x_c)/h)  * (1.0 / np.cosh(n*np.fabs(z-z_c)/(2.*z0)) )**(2./n)

def read_profiles(profiles_file):
    f=open(profiles_file, 'r')
    lines = f.read()
    cuts = lines.split('#END\n')

    Z = []; R = []; I = []
    for cut in cuts:
      if cut!='':
        cut = cut.split('\n')
        Z.append(float( cut[0].split('#Z:')[-1] ))
        rr = []; II = []

        for k in range(1, len(cut)):
            if cut[k]!='':
                r,i = cut[k].split('\t')
                rr.append(float(r))
                II.append(float(i))
        R.append(rr)
        I.append(II)
    return np.array(Z),np.array(R),np.array(I)
    

def main(profiles_file, results_file, xc, psf_image = None, units = ['arcsec',0.75], mode='triple'):
    if psf_image is not None:
        kernel = fit_eon_radial.psf_slice(psf_image, axis='xaxis')
    else:
        kernel = None
        
    Z,radius,I = read_profiles(profiles_file)
    if I[0][-1]>I[0][0]:
        radius = np.fabs(radius - np.max(radius)) # WARNING: For LEFT!!!!!!!!!!!!

          
    Z,L0,L0_err,h,h_err,L01,L01_err,h1,h1_err,Rbr1,L02,L02_err,h2,h2_err,Rbr2,L03,L03_err,h3,h3_err = np.loadtxt(results_file, usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],dtype=float, unpack=True,skiprows=1,delimiter='\t')


    fig = plt.figure(0, figsize=(5, 10))
    ax=plt.subplot()
    
    M0 = 28.
    min_mag = []
    max_mag = []
    max_r = []
    for k in range(len(Z)):
        mag = M0 - np.log10(I[k])       
        ax.plot(radius[k], mag,'o', markersize=9,color='gray',markeredgecolor='gray')

        if mode=='single':
            I_single = disk_edge_soph(0., radius[k], 2.*h[k]*L0[k], 2.5, h[k], xc, 0., 1.)
            if kernel is not None:
                signal.fftconvolve(I_single, kernel, mode='same')
            
            #print(I_single)
            #print(I[k])
            #exit()
            mag_single = M0 - np.log10(I_single)
            ax.plot(radius[k], mag_single, '-', color='red')
        '''
        else:
            I_1 = disk_edge_soph(Z[k], 2.*10.*L01[k], Z01[k], Y0[k], N1[k])
            I_2 = disk_edge_soph(Z[k], 2.*10.*L02[k], Z02[k], Y0[k], N2[k])
            I_sum = I_1 + I_2
            mag_1 = M0 - np.log10(I_1)
            mag_2 = M0 - np.log10(I_2)     
            mag_sum = M0 - np.log10(I_sum)     
            ax.plot(Z[k]-Y0[k], mag_sum, '-', color='black')
            ax.plot(Z[k]-Y0[k], mag_1, '--', color='blue')            
            ax.plot(Z[k]-Y0[k], mag_2, '--', color='red')
        '''    
        min_mag.append(np.nanmin(mag))
        max_mag.append(np.nanmax(mag))
        max_r.append(np.nanmax(np.fabs(radius[k])))
        M0 = M0 - 2

    Min_mag = []
    Max_mag = []
    for k in range(len(min_mag)):
        if not np.isinf(min_mag[k]):
            Min_mag.append(min_mag[k])

    for k in range(len(max_mag)):
        if not np.isinf(max_mag[k]):
            Max_mag.append(max_mag[k])
    
    ax.set_xlim(0, np.nanmax(max_r))
    ax.set_ylim(np.nanmax(Max_mag)+1, np.nanmin(Min_mag)-1)
    
    ax.set_xlabel(r' $R$ (%s)' % (units[0]), fontsize=fsize)
    ax.set_ylabel(r' $\mu$ (mag\,arcsec$^{-2}$)', fontsize=fsize)
    
    plt.savefig('r_porfiles.eps', transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.clf()
    plt.close()

#main('radial_profiles_1.txt', 'radial_fits_1.dat', 651., psf_image='psf.fits', mode='single') # NGC3628
#main('radial_profiles_1.txt', 'radial_fits_1.dat', 366., psf_image='psf.fits',  mode='single') # NGC4302
#main('radial_profiles_1.txt', 'radial_fits_1.dat', 772, psf_image='psf.fits',  mode='single') # NGC891
#main('radial_profiles_1.txt', 'radial_fits_1.dat', 874, psf_image=None,  mode='single') # NGC4302 model (no psf)
#main('radial_profiles_1.txt', 'radial_fits_1.dat', 251, psf_image='psf.fits',  mode='single') # NGC4302 model (no psf)