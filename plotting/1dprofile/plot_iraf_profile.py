#!/usr/bin/python
# -*- coding:  cp1251 -*-


# Import standard modules
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
import argparse
import warnings

from astropy.stats import sigma_clipped_stats

def smoothTriangle(data, data_err, degree, dropVals=False):
    #https://plot.ly/python/smoothing/
    triangle=np.array(list(range(degree)) + [degree] + list(range(degree)[::-1])) + 1
    
    smoothed=[]

    for i in range(degree, len(data) - degree * 2):
        point=data[i:i + len(triangle)] * triangle 
        smoothed.append(sum(point)/sum(triangle))
    if dropVals:
        return smoothed
    smoothed=[smoothed[0]]*int(degree + degree/2) + smoothed
    while len(smoothed) < len(data):
        smoothed.append(smoothed[-1])
    return smoothed

def create_average_profile(r, I, I_err, bin_box=10):
    #I = smooth(I,bin_box)
    I = smoothTriangle(I, I_err, bin_box)
    #I = regr_fit(r,I)
    return r, I

def find_lim_radius(x,y,bin_box):
    r_lim = []
    for k in range(len(x)-1):
        if fabs(y[k]-y[k+1])/y[k]<0.01 and k>bin_box+2:
            if len(r_lim)<=5:
                r_lim.append(x[k])
            else:
                break
    if r_lim==[]:
        r_lim = max(x)
    print('Radii with SB within 1 per cent:',r_lim)
    return max(r_lim)
            
def plot_azim(ax,a,inten,inten_err,pix2sec,m0,rmax,outp_format='eps',xmax=None,ymin=None,ymax=None,text=None, bin_box=3,plot_add=False):
        mag = m0 - 2.5*log10(inten)+ 5.*log10(pix2sec)
        mag_err = fabs((2.5/log(10.0)) * inten_err/inten)
        fsize=16

        if True:
           color = 'gray'
           color_bar = 'black'
           symb = 'o'

        ax.errorbar(a*pix2sec,mag,yerr=mag_err,fmt=symb,color=color,ecolor=color, markersize=5)
        ax.plot(a*pix2sec, mag, symb,color=color,markeredgecolor=color, markersize=3)
        mean_mag, median_mag, std_mag = sigma_clipped_stats(mag, sigma=3.0)
        if plot_add:
            r,I = create_average_profile(a, inten, inten_err, bin_box=bin_box)
            ax.plot(r*pix2sec, m0 - 2.5*log10(I)+ 5.*log10(pix2sec), '*',color='red',markeredgecolor='red', markersize=3)
            r_lim = find_lim_radius(r,I,bin_box)
            ax.axvline(x=r_lim*pix2sec)
            print('Maximum radius [pix]: %.1f' % (r_lim) )
        else:
            r_lim = float('nan')
            
        if ymin is None:  
           ax.set_ylim(max(mag)+std_mag,min(mag)-std_mag)
        else:
           ax.set_ylim(ymax,ymin)
        ax.set_xlim(0.,rmax*pix2sec)



        ax.set_xlabel(r'$r$ (arcsec)', fontsize=fsize)
        ax.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
        if text is not None:
           ax.text(0.7, 0.9, text, fontsize=fsize, color='black',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline')



        if True:  
           plt.savefig("azim_aver."+outp_format, transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
           plt.clf()
           plt.close()
        return r_lim
    
def main(ellipse_file,pix2sec,m0,plot_add=False, bin_box=3):
        sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = loadtxt(ellipse_file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 5, dtype='str')

        
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
    
        fig = figure(1,figsize=(5,5))
        ax = fig.add_subplot(111)
        r_lim = plot_azim(ax,sma,inten,inten_err,pix2sec,m0,rmax=max(sma),outp_format='eps',xmax=None,ymin=None,ymax=None,text=None,bin_box=bin_box,plot_add=plot_add)
        return r_lim

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot azimuthally averaged profile")
    parser.add_argument("ellipse_file", help="Input IRAF/ELLIPSE file with results")
    parser.add_argument("--scale", help="Pixel scale in arcsec/pix", type=float, default=1.)
    parser.add_argument("--m0", help="Zero-point", type=float, default=28.)
    parser.add_argument("--plot_add", help="Overplot averaged profile and the maximum radius", action="store_true", default=False)
    
    args = parser.parse_args()

    ellipse_file = args.ellipse_file
    m0 = args.m0
    pix2sec = args.scale
    plot_add = args.plot_add

    main(ellipse_file,pix2sec,m0,plot_add)