#!/usr/bin/python
# DESCRIPTION:
# Function to compare IRAF/ELLIPSE results and results from my function to create azimuthally averaged
# profile with fixed PA and ellipticity.
# MINIMAL USAGE: python  compare_iraf_and_ellipse.py [iraf_output] [ellipse_output]
# EXAMPLE:  python compare_iraf_and_ellipse.py ellipse.txt azim_average.txt

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from numpy import *
from pylab import *
import os
import shutil
import astropy.io.fits as pyfits
from astropy.stats import sigma_clipped_stats
import argparse


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

tmp_out = sys.stdout


fsize=16



#def plot_azim(ax,a,inten,inten_err,pix2sec,m0,rmax,overplot,outp_format,comp,xmax=None,ymin=None,ymax=None):
def plot_azim(ax, sma, inten, inten_err, sma_azim,inten_azim,inten_azim_err, pix2sec=1., m0=28.):
        mag = m0 - 2.5*log10(inten)+ 5.*log10(pix2sec)
        mag_err = fabs((2.5/log(10.0)) * inten_err/inten)

        mag_azim = m0 - 2.5*log10(inten_azim)+ 5.*log10(pix2sec)
        mag_azim_err = fabs((2.5/log(10.0)) * inten_azim_err/inten_azim)


        ax.errorbar(sma*pix2sec, mag, yerr=mag_err, fmt='o', color='black', ecolor='black', markersize=5, label='IRAF')
        ax.plot(sma*pix2sec, mag, 'o', color='black', markeredgecolor='black', markersize=3)
        
        ax.errorbar(sma_azim*pix2sec, mag_azim, yerr=mag_azim_err, fmt='o', color='red', ecolor='red', markersize=5, label='Fixed')
        ax.plot(sma_azim*pix2sec, mag_azim, 'o', color='red', markeredgecolor='red', markersize=3)        
        
        
        mean_mag, median_mag, std_mag = sigma_clipped_stats(mag, sigma=3.0, iters=5)
        mean_mag_azim, median_mag_azim, std_mag_azim = sigma_clipped_stats(mag_azim, sigma=3.0, iters=5)

        ax.set_xlim(0.,max([np.max(sma*pix2sec), np.max(sma_azim*pix2sec)]))

        ax.set_ylim( max([max(mag)+std_mag, max(mag_azim)+std_mag_azim]), min([min(mag)-std_mag, min(mag_azim)-std_mag_azim]))

                

        ax.set_xlabel(r'$r$ (arcsec)', fontsize=fsize)
        ax.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
        
        plt.legend(loc=0, borderaxespad=0.,fontsize=10,numpoints=1,markerscale=2)
        
        
 
        plt.savefig("compar_profiles.png", transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
        plt.clf()
        plt.close()





def main(iraf_ell_file, azim_model_file, pix2sec=1.0, m0=28.0):
    sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = loadtxt(iraf_ell_file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 5, dtype='str')

    sma_azim,inten_azim,inten_azim_err = np.loadtxt(azim_model_file, usecols=[0,1,2], dtype=float, unpack=True, skiprows = 1, delimiter='\t')

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




    fig = figure(1,figsize=(5,5))
    ax = fig.add_subplot(111)    
    
    
    plot_azim(ax, sma, inten, inten_err, sma_azim,inten_azim,inten_azim_err, pix2sec, m0)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Comparison between IRAF/ELLIPSE and azimuthally averaged fitting with constant PA and ellipticity.")
    parser.add_argument("iraf_file", help="IRAF file with ELLIPSE fitting results")
    parser.add_argument("azim_file", help="File with average fitting results") 

    parser.add_argument("--ZP", nargs='?', const=1, help="Optional: Input the y-coordinate of the centre", type=float, default=28.) 
    parser.add_argument("--pix2sec", nargs='?', const=1, help="Optional: Input the y-coordinate of the centre", type=float, default=1.) 

    args = parser.parse_args()

    iraf_file = args.iraf_file
    azim_file = args.azim_file
    ZP = args.ZP
    pix2sec = args.pix2sec

    
    main(iraf_file, azim_file, pix2sec, ZP)
