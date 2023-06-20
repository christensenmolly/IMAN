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

axes_fontsize = 20
axes_label_fontsize = 0.8*axes_fontsize

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
            
def plot_azim(ax,a,inten,inten_err,pix2sec,m0,rmax,output_image='azim_aver.png',xmax=None,ymin=None,ymax=None,text=None, bin_box=3, plot_add=False, ylabel=r'$\mu$ (mag arcsec$^{-2}$)', xlabel=r'$r$ (arcsec)'):
        if m0 is not None:
            mag = m0 - 2.5*log10(inten)+ 5.*log10(pix2sec)
            mag_err = fabs((2.5/log(10.0)) * inten_err/inten)
        else:
            mag = inten
            mag_err = inten_err          
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
        
        if m0 is not None:
            if ymin is None:  
                ax.set_ylim(max(mag)+std_mag,min(mag)-std_mag)
            else:
                ax.set_ylim(ymax,ymin)
        else:
            if ymin is None:  
                ax.set_ylim(min(mag)-std_mag,max(mag)+std_mag)
            else:
                ax.set_ylim(ymin,ymax)            
        
        
        ax.set_xlim(0.,rmax*pix2sec)



        ax.set_xlabel(xlabel, fontsize=fsize)
        ax.set_ylabel(ylabel, fontsize=fsize)
        if text is not None:
           ax.text(0.7, 0.9, text, fontsize=fsize, color='black',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline')

        ax.minorticks_on()
        ax.tick_params(direction='in', length=3, which='minor', top=True, right=True, labelsize=axes_label_fontsize)
        ax.tick_params(direction='in', length=6, which='major', top=True, right=True, labelsize=axes_label_fontsize)


        if True:  
           plt.savefig(output_image, transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
           plt.clf()
           plt.close()
        return r_lim
    
def main(azim_model_file,pix2sec,m0,plot_add=False, bin_box=3, ylabel=r'$\mu$ (mag arcsec$^{-2}$)', xlabel=r'$r$ (arcsec)', text=None, output_image="azim_aver.png"):
    a,inten,inten_err = np.loadtxt(azim_model_file, usecols=[0,1,2], dtype=float, unpack=True, skiprows = 1, delimiter='\t')
    
    fig = figure(1,figsize=(5,5))
    ax = fig.add_subplot(111)
    r_lim = plot_azim(ax,a,inten,inten_err,pix2sec,m0,rmax=max(a),output_image=output_image,xmax=None,ymin=None,ymax=None,text=text,bin_box=bin_box,plot_add=plot_add, xlabel=xlabel, ylabel=ylabel)
    return r_lim

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot azimuthally averaged profile")
    parser.add_argument("azim_model", help="Input azimuthally averaged profile (txt file)")
    parser.add_argument("--scale", help="Pixel scale in arcsec/pix", type=float, default=1.)
    parser.add_argument("--m0", help="Zero-point", type=float, default=None)
    parser.add_argument("--plot_add", help="Overplot averaged profile and the maximum radius", action="store_true", default=False)
    parser.add_argument("--xlabel", help="X label", type=str, default='$r$ (arcsec)')
    parser.add_argument("--ylabel", help="Y label", type=str, default='$\mu$ (mag arcsec$^{-2}$)')
    parser.add_argument("--text", help="Text in plot", type=str, default=None)
    parser.add_argument("--o", help="Output image", type=str, default="azim_aver.png")
    args = parser.parse_args()

    azim_model = args.azim_model
    m0 = args.m0
    pix2sec = args.scale
    plot_add = args.plot_add
    xlabel = args.xlabel
    ylabel = args.ylabel
    text = args.text
    output_image = args.o

    main(azim_model,pix2sec,m0,plot_add, xlabel=xlabel, ylabel=ylabel, text=text, output_image=output_image)