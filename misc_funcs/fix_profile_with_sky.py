#!/usr/bin/python

from astropy.io import fits as pyfits
import numpy as np
import math
import itertools
from scipy import ndimage
import sys
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import NullFormatter
from itertools import product
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import os
import subprocess
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats
#from photutils import data_properties, properties_table
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from astropy import wcs
from astroquery import ned

def create_profile(a, inten, m0, pix2sec, inten_fixed=None):
    #pix2sec = 0.396
    #m0 = 28.294
    mag = m0 - 2.5*np.log10(inten)+ 5.*np.log10(pix2sec)
    

    fig = plt.figure(1,figsize=(5,5))
    ax = fig.add_subplot(111)
    fsize=16        
    rmax = max(a)+5

    mean_mag, median_mag, std_mag = sigma_clipped_stats(mag, sigma=3.0)
    ax.set_ylim(32.,min(mag)-std_mag)
    ax.set_xlim(0.,rmax*pix2sec)


    ax.set_xlabel(r'$r$ (arcsec)', fontsize=fsize)
    ax.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
       

    
    ax.plot(a*pix2sec, mag, 'o',color='gray',markeredgecolor='gray', markersize=3)
    if inten_fixed is not None:
        mag_fixed = m0 - 2.5*np.log10(inten_fixed)+ 5.*np.log10(pix2sec)
        ax.plot(a*pix2sec, mag_fixed, 'X',color='red',markeredgecolor='red', markersize=3)
    plt.grid()
    plt.show()
    plt.clf()
    plt.close()



def main(r, I, m0, pix2sec):
    create_profile(r, I, m0, pix2sec, inten_fixed=None)
    s = input("Input measured SB level and estimated SB level (using commas):")
    #print(s)
    s = s.split(',')
    SB_inp = s[0]
    SB_est = s[1]
    #print(SB_inp,SB_est)
    #exit()
    delta = (10**(0.4*m0)) * ( 10**(-0.4*float(SB_est)) - 10**(-0.4*float(SB_inp)) ) * (pix2sec**2)
    new_I = I + delta
    create_profile(r, I, m0, pix2sec, inten_fixed=new_I)
    return float(SB_inp),float(SB_est),delta

'''
#f = open('fix_sky.dat','w')
#f.close()
for galaxy in ['NGC278','NGC247','NGC470','NGC509','NGC524','NGC525','NGC661','NGC1023','NGC1084','NGC1400','NGC1407','NGC2549','NGC3031','NGC3379','NGC3521','NGC3623','NGC3627','NGC4406','NGC4449','NGC4594','NGC4621','NGC4649','NGC4697','NGC5055','NGC5236','NGC5576','NGC5719','NGC5806','NGC5811','NGC5813','NGC5866','NGC6340','NGC6643','NGC7217','NGC467','M51','NGC3032','NGC3079','NGC3384','NGC4754','NGC5814','NGC7463','NGC7465','NGC7743','NGC205','NGC5577']:
  if galaxy == 'NGC5577':
    f = open('fix_sky.dat','a')
    a,inten,inten_err = np.loadtxt('/home/amosenko/Toshiba_1/CurrentWork/Rich/Paper_sample/%s/ellipse_phot/azim_model.txt' % (galaxy), usecols=[0,1,2], dtype=float, unpack=True, skiprows = 1, delimiter='\t')
    SB_inp,SB_est,delta = main(a, inten)
    print('%s: %f' % (galaxy, delta))
    f.write('%s\t%.2f\t%.2f\t%f\n' % (galaxy, SB_inp,SB_est,delta))
    f.close()
'''

ell_file = 'ellipse.txt'
if True:            
            sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = np.loadtxt(ell_file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 5, dtype='str')

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


SB_inp,SB_est,delta = main(sma, inten, 28.294, 0.396)
print('delta=%f' % (delta))