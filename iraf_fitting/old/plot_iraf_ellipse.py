#!/usr/bin/python
# DESCRIPTION:
# STSDAS IRAF/ELLIPSE fitting. IRAF is required. Do not create your own login.cl where you launch this script. 
# MINIMAL USAGE: python  iraf_ellipse.py [input_image]
# EXAMPLE:  python iraf_ellipse.py galaxy_clean_galf.fits --xc 1122 --yc 1021 --maxsma 850 --ZP 29.5962 --pix2sec 0.833 --mask mask.fits

import sys
import math
import numpy as np
#from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import astropy.io.fits as pyfits
from astropy.stats import sigma_clipped_stats
import argparse
import glob
import warnings
warnings.filterwarnings("ignore")


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

tmp_out = sys.stdout


fsize=30
FNULL = open(os.devnull, 'w')



def plot_iso(ax1, ax2, ax3, output_image, PA, errPA, ell, errell, B4, errB4, a, pix2sec, PA_lim=None, ell_lim=None, B4_lim=None, Rlim=None):
        color = 'gray'
        color_bar = 'black'
        symb = 'o'  
        
        PA[PA<0.] = PA[PA<0.] + 180.
        
        if PA_lim is not None:
            [min_PA, max_PA] = PA_lim
        else:
            mean_errPA, median_errPA, std_errPA = sigma_clipped_stats(errPA, sigma=3.0)
            mean_PA, median_PA, std_PA = sigma_clipped_stats(PA, sigma=3.0)
            min_PA = np.min(fabs(PA)-errPA)*0.99
            max_PA = np.max(fabs(PA)+errPA)*1.01
            
            
        if ell_lim is not None:
            [min_ell,max_ell] = ell_lim
        else:
            min_ell = np.min(fabs(ell)-errell)*0.99
            max_ell = np.max(fabs(ell)+errell)*1.01
            
            
        
        if B4_lim is not None:
           [min_B4, max_B4] = B4_lim 
        else:
            mean_B4, median_B4, std_B4 = sigma_clipped_stats(B4, sigma=3.0)
            if fabs(median_B4-5*std_B4)<2. and fabs(median_B4+5*std_B4)<2.:
              min_B4 = median_B4 - 5*std_B4
              max_B4 = median_B4 + 5*std_B4
            else:
              min_B4 = -2
              max_B4 = 2


            

        ax1.errorbar(a*pix2sec, fabs(PA),yerr=errPA,fmt=symb,color=color_bar, markersize=3)
        ax1.plot(a*pix2sec, fabs(PA),symb,color=color,markeredgecolor=color, markersize=3)
        ax1.set_ylabel(r'PA (deg)', fontsize=fsize)
        ax1.set_ylim(float(min_PA), float(max_PA))


        ax2.errorbar(a*pix2sec, fabs(ell),yerr=errell,fmt=symb,color=color_bar, markersize=3)
        ax2.plot(a*pix2sec, fabs(ell),symb,color=color,markeredgecolor=color, markersize=3)
        ax2.set_ylabel(r'$\epsilon$ ', fontsize=fsize)
        ax2.set_ylim(float(min_ell), float(max_ell))


        ax3.errorbar(a*pix2sec, B4,yerr=errB4,fmt=symb,color=color_bar, markersize=3)        
        ax3.plot(a*pix2sec, B4, symb,color=color,markeredgecolor=color, markersize=3)
        ax3.set_xlabel(r'r (arcsec)', fontsize=fsize)
        ax3.set_ylabel(r'B$_4$ ', fontsize=fsize)
        ax3.set_ylim(float(min_B4),float(max_B4))
        ax3.axhline(y=0., ls='--')
        

        if Rlim is not None:
            [min_R,max_R] = Rlim
            ax1.set_xlim(float(min_R), float(max_R))            
            ax2.set_xlim(float(min_R), float(max_R))
            ax3.set_xlim(float(min_R), float(max_R))


  
        xticklabels = ax1.get_xticklabels()+ax2.get_xticklabels()
        setp(xticklabels, visible=False)

        for ax in [ax1,ax2,ax3]:
          ax.get_xaxis().set_tick_params(direction='in', width=1)
          ax.get_yaxis().set_tick_params(direction='in', width=1)
          for tick in ax.xaxis.get_major_ticks():
              tick.label.set_fontsize(fsize-6)
              ax.set_xlim(0.,np.max(a)*pix2sec)
          for tick in ax.yaxis.get_major_ticks():
              tick.label.set_fontsize(fsize-8)
        

        plt.savefig(output_image, transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
        plt.clf()
        plt.close() 




def plot_azim(ax,a,inten,inten_err,pix2sec,m0,rmax,overplot,outp_format,comp,xmax=None,ymin=None,ymax=None):
        # Function to plot azimuthally averaged profile
        mag = m0 - 2.5*log10(inten)+ 5.*log10(pix2sec)
        mag_err = fabs((2.5/log(10.0)) * inten_err/inten)


        if overplot==True:
          color = 'gray'
          color_bar = 'black'
          symb = 'o'
        else:
          color = 'red'
          color_bar = 'red'
          symb = '*'
          
        if comp==False:
          color = 'gray'
          color_bar = 'black'
          symb = 'o'  

        if comp==True:
          if overplot==True:
            label='galaxy'
          else:
            label='model'
        else:
            label=''

        ax.errorbar(a*pix2sec,mag,yerr=mag_err,fmt=symb,color=color,ecolor=color, markersize=5,label=label)
        ax.plot(a*pix2sec, mag, symb,color=color,markeredgecolor=color, markersize=3)
        mean_mag, median_mag, std_mag = sigma_clipped_stats(mag, sigma=3.0)

        if overplot==True:  
            ax.set_ylim(max(mag)+std_mag,min(mag)-std_mag)
            ax.set_xlim(0.,rmax*pix2sec)
        else:
            if xmax==None:
                ax.set_xlim(0.,rmax*pix2sec)
            else:
                ax.set_xlim(0.,xmax)
            if ymax==None:
                ax.set_ylim(max(mag)+std_mag,min(mag)-std_mag)
            else:
                ax.set_ylim(ymax,ymin)
                
            #mag_ylim = [max(mag)+std_mag,min(mag)-std_mag]

        ax.set_xlabel(r'$r$ (arcsec)', fontsize=fsize)
        ax.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
        
        if comp==True:
          plt.legend(loc=0, borderaxespad=0.,fontsize=10,numpoints=1,markerscale=2)
        
        ax.get_xaxis().set_tick_params(direction='in', width=1)
        ax.get_yaxis().set_tick_params(direction='in', width=1)
        if overplot==False:  
            plt.savefig("azim_aver."+outp_format, transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
            plt.clf()
            plt.close()
            return None,None,None
        else:
            return rmax*pix2sec, min(mag)-std_mag,max(mag)+std_mag


def main(ellipse_file, output_image_azim, output_image_ell, pix2sec, m0, PA_lim=None, ell_lim=None, B4_lim=None, Rlim=None):
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

        f= figure(0,figsize=(4,10))

        gs = gridspec.GridSpec(3, 1)#,width_ratios=[1,3])
        gs.update(left=0.25, right=3.00, hspace=0.0)

        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        ax3 = plt.subplot(gs[2])

        
        plot_iso(ax1, ax2, ax3, output_image_ell, PA, errPA, ell, errell, B4, errB4, sma, pix2sec, PA_lim=PA_lim, ell_lim=ell_lim, B4_lim=B4_lim, Rlim=Rlim)
        
        exit()
        
        fig = figure(1,figsize=(5,5))
        ax4 = fig.add_subplot(111)    
        plot_azim(ax4,sma,inten,inten_err,pix2sec,m0,float(maxsma),overplot,outp_format,comp)

        print('Done!')
        





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="IRAF/ELLIPSE fitting")
    parser.add_argument("ellipse_file", help="Input IRAF/ELLIPSE file with results")
    parser.add_argument("--oell", nargs='?', const=1, help="Optional: Output picture with Ellipse results", type=str, default='iraf.eps') 
    parser.add_argument("--oazi", nargs='?', const=1, help="Optional: Output picture with azimuthal averaging", type=str, default='azim_aver.eps') 
    parser.add_argument("--PAlim", nargs='?', const=1, help="Optional: PA limits separated by comma", type=str, default=None) 
    parser.add_argument("--Elim", nargs='?', const=1, help="Optional: Ellipticity limits separated by comma", type=str, default=None) 
    parser.add_argument("--B4lim", nargs='?', const=1, help="Optional: B4 limits separated by comma", type=str, default=None) 
    parser.add_argument("--Rlim", nargs='?', const=1, help="Optional: R limits separated by comma", type=str, default=None) 

    parser.add_argument("--ZP", nargs='?', const=1, help="Optional: Input the y-coordinate of the centre", type=float, default=28.) 
    parser.add_argument("--pix2sec", nargs='?', const=1, help="Optional: Input the y-coordinate of the centre", type=float, default=1.) 

    args = parser.parse_args()

    ellipse_file = args.ellipse_file
    output_image_azim = args.oazi
    output_image_ell = args.oell
    
    PA_lim = args.PAlim
    ell_lim = args.Elim
    B4_lim = args.B4lim
    Rlim = args.Rlim

    ZP = args.ZP
    pix2sec = args.pix2sec
    
    if PA_lim is not None:
        PA_lim = PA_lim.split(',')
        
    if ell_lim is not None:
        ell_lim = ell_lim.split(',')
        
    if B4_lim is not None:
        B4_lim = B4_lim.split(',')

    if Rlim is not None:
        Rlim = Rlim.split(',')
        
    main(ellipse_file, output_image_azim, output_image_ell, pix2sec=pix2sec, m0=ZP, PA_lim=PA_lim, ell_lim=ell_lim, B4_lim=B4_lim, Rlim=Rlim)
