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

    Z = []; R = []; I = []; I_MOD = []
    for cut in cuts:
      if cut!='':
        cut = cut.split('\n')
        Z.append(float( cut[0].split('#Z:')[-1] ))
        rr = []; II = []; I_mod = []

        for k in range(1, len(cut)):
            if cut[k]!='':
                r,i,ii = cut[k].split('\t')
                rr.append(float(r))
                II.append(float(i))
                I_mod.append(float(ii))
        R.append(rr)
        I.append(II)
        I_MOD.append(I_mod)
    return np.array(Z),np.array(R),np.array(I),np.array(I_MOD)
    

def main(profiles_file, units = ['arcsec',0.75], mode='triple',mag_difference = 0.5):
  if len(profiles_file)==1:     
    Z,radius,I,I_mod = read_profiles(profiles_file[0])
    if I[0][-1]>I[0][0]:
        radius = np.fabs(radius - np.max(radius)) # WARNING: For LEFT!!!!!!!!!!!!


    fig = plt.figure(0, figsize=(5, 10))
    ax=plt.subplot()
    
    M0 = 28.
    min_mag = []
    max_mag = []
    max_r = []
    min_r = []
    for k in range(len(Z)):
        mag = M0 - np.log10(I[k])       
        ax.plot(radius[k], mag,'o', markersize=5,color='gray',markeredgecolor='gray')

        if mode=='single':
            mag_single = M0 - np.log10(I_mod[k])
            ax.plot(radius[k], mag_single, '-', color='red')
 
        min_mag.append(np.nanmin(mag))
        max_mag.append(np.nanmax(mag))
        max_r.append(np.nanmax(np.fabs(radius[k])))
        min_r.append(np.nanmin(np.fabs(radius[k])))
        M0 = M0 - mag_difference

    Min_mag = []
    Max_mag = []
    for k in range(len(min_mag)):
        if not np.isinf(min_mag[k]):
            Min_mag.append(min_mag[k])

    for k in range(len(max_mag)):
        if not np.isinf(max_mag[k]):
            Max_mag.append(max_mag[k])
    
    ax.set_xlim(np.nanmin(min_r), np.nanmax(max_r))
    ax.set_ylim(np.nanmax(Max_mag)+mag_difference, np.nanmin(Min_mag)-mag_difference)
  else:     
    Z_in,radius_in,I_in,I_mod_in = read_profiles(profiles_file[0])
    Z_out,radius_out,I_out,I_mod_out = read_profiles(profiles_file[1])
    
    if I_in[0][-1]>I_in[0][0]:
        radius_in = np.fabs(radius_in - np.max(radius_in)) # WARNING: For LEFT!!!!!!!!!!!!

    if I_out[0][-1]>I_out[0][0]:
        radius_out = np.fabs(radius_out - np.max(radius_out)) # WARNING: For LEFT!!!!!!!!!!!!

    fig = plt.figure(0, figsize=(5, 10))
    ax=plt.subplot()
    
    M0 = 28.
    min_mag = []
    max_mag = []
    max_r = []
    min_r = []
    Z_done = []
    for k in range(len(Z_in)):
      if Z_in[k] not in Z_done:
        mag_in = M0 - np.log10(I_in[k])       
        ax.plot(radius_in[k], mag_in,'o', markersize=5,color='gray',markeredgecolor='gray')

        if mode=='single':
            mag_single_in = M0 - np.log10(I_mod_in[k])
            ax.plot(radius_in[k], mag_single_in, '-', color='red')
 
      
        
        
        if Z_in[k] in Z_out:
            ind = list(Z_out).index(Z_in[k])
            mag_out = M0 - np.log10(I_out[ind])       
            ax.plot(radius_out[ind], mag_out,'o', markersize=5,color='gray',markeredgecolor='gray')

            if mode=='single':
                mag_single_out = M0 - np.log10(I_mod_out[ind])
                ax.plot(radius_out[ind], mag_single_out, '-', color='blue')
            
        
        min_mag.append(np.nanmin([np.nanmin(mag_in),np.nanmin(mag_out)]))
        max_mag.append(np.nanmax([np.nanmin(mag_in),np.nanmax(mag_out)]))      
        max_r.append(np.nanmax([np.nanmax(np.fabs(radius_in[k])),np.nanmax(np.fabs(radius_out[ind]))]))
        min_r.append(np.nanmin([np.nanmin(np.fabs(radius_in[k])),np.nanmin(np.fabs(radius_out[ind]))]))
        M0 = M0 - mag_difference
        Z_done.append(Z_in[k])
    
    Min_mag = []
    Max_mag = []
    for k in range(len(min_mag)):
        if not np.isinf(min_mag[k]):
            Min_mag.append(min_mag[k])

    for k in range(len(max_mag)):
        if not np.isinf(max_mag[k]):
            Max_mag.append(max_mag[k])
    
    ax.set_xlim(np.nanmin(min_r), np.nanmax(max_r))
    ax.set_ylim(np.nanmax(Max_mag)+mag_difference, np.nanmin(Min_mag)-mag_difference)






    
    ax.set_xlabel(r' $R$ (%s)' % (units[0]), fontsize=fsize)
    ax.set_ylabel(r' $\mu$ (mag\,arcsec$^{-2}$)', fontsize=fsize)
    
    plt.savefig('r_porfiles.eps', transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.clf()
    plt.close()



# Test_eon/one_disc_galaxy
#main('radial_profiles_1_single.txt',  mode='single')




