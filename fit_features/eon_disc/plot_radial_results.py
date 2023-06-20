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
from scipy import stats
from math import *
warnings.filterwarnings("ignore")
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


fsize = 18


def find_regression(xx, yy):
        k1, b1, r_value, p_value, std_err = stats.linregress(xx,yy)
        k2, b2, r_value, p_value, std_err = stats.linregress(yy,xx)

        t =  (1.0-k1*k2)/(k1+k2)
        t2 =  t/(1.0+sqrt(1.0+t*t))
        k4 = (t2+k1)/(1.0-t2*k1)
        xc = ((b2+b1*k2)/(1.0-k1*k2))
        yc = k1*xc + b1
        b4 = yc - k4*xc
        return k1,b1,k2,b2,k4,b4,r_value,std_err
    
def RegrLine(xxx, yyy, xx1, xx2, yy1, yy2, i, sigma_outliers=3.):
        xx = []
        yy = []
        if xx1>xx2:
            xx1,xx2 = xx2,xx1
        if yy1>yy2:
            yy1,yy2 = yy2,yy1


        for k in range(0,len(xxx)):
            if xxx[k]>xx1 and xxx[k]<xx2 and yyy[k]>yy1 and yyy[k]<yy2:
                xx.append(xxx[k])
                yy.append(yyy[k])

        k1,b1,k2,b2,k4,b4,r_value,std_err = find_regression(xx,yy)

        if sigma_outliers is not None:
            XX = []; YY = []; Dist = []
            for k in range(len(xx)):
                Dist.append( abs(k4*xx[k]-yy[k]+b4) / ( sqrt(k4*k4+1.) ) )
                
            std_dist = np.std(Dist)    
            for k in range(len(xx)):                
                if Dist[k]<=sigma_outliers*std_dist:
                    XX.append(xx[k])
                    YY.append(yy[k])
            k1,b1,k2,b2,k4,b4,r_value,std_err = find_regression(XX,YY)    

        if i==1:
                return k1,b1,r_value
        if i==2:
                return 1./k2,-b2/k2,r_value
        if i==3:
                return k4,b4,r_value


def main(input_files, Ropt=None, units = ['arcsec',0.75], regression_lims = None, xlim = None, ylim = None):
    # [minx, maxx, miny, maxy] = regression_lims - limits to plot regression line (should be in pix!)
    
    # First file should be left
    # Second file should be right
    Z,L0,L0_err,h,h_err,L01,L01_err,h1,h1_err,Rbr1,L02,L02_err,h2,h2_err,Rbr2,L03,L03_err,h3,h3_err = np.loadtxt(input_files[0], usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],dtype=float, unpack=True,skiprows=1,delimiter='\t')
    Z = Z * units[1]
    h = h * units[1]
    h_err = h_err * units[1]

    h1 = h1 * units[1]
    h1_err = h1_err * units[1]
    Rbr1 = Rbr1 * units[1]

    h2 = h2 * units[1]
    h2_err = h2_err * units[1]
    Rbr2 = Rbr2 * units[1]

    h3 = h3 * units[1]
    h3_err = h3_err * units[1]
    #Rbr3 = Rbr3 * units[1]
    
    if len(input_files)==2:
        Z_r,L0_r,L0_err_r,h_r,h_err_r,L01,L01_err,h1,h1_err,Rbr1,L02,L02_err,h2,h2_err,Rbr2,L03,L03_err,h3,h3_err = np.loadtxt(input_files[1], usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],dtype=float, unpack=True,skiprows=1,delimiter='\t')
        Z_r = Z_r * units[1]
        h_r = h_r * units[1]
        h_err_r = h_err_r * units[1]

        h1_r = h1 * units[1]
        h1_err_r = h1_err * units[1]
        Rbr1_r = Rbr1 * units[1]

        h2_r = h2 * units[1]
        h2_err_r = h2_err * units[1]
        Rbr2_r = Rbr2 * units[1]

        h3_r = h3 * units[1]
        h3_err_r = h3_err * units[1]
        #Rbr3_r = Rbr3 * units[1]
    
    
    
    
    ########### PLOTTING: SINGLE DISK ########### 
    
    #### 1) h Vs Z
    fig = plt.figure(0, figsize=(6, 5))
    ax=plt.subplot()

    if not np.isnan(h_err[0]):
      if len(input_files)==2:
        ax.errorbar(Z, h, yerr=h_err,fmt='o',markersize=9,color='red',markeredgecolor='red', ecolor='red', capthick=2, label='Inner', zorder=0)
        ax.errorbar(Z_r, h_r, yerr=h_err_r,fmt='^',markersize=9,color='blue',markeredgecolor='blue', ecolor='blue', capthick=2, label='Outer', zorder=0)
      else:
        ax.errorbar(Z, h, yerr=h_err,fmt='o',markersize=9,color='red',markeredgecolor='red', ecolor='red', capthick=2, zorder=1)
    else:
      if len(input_files)==2:
        ax.plot(Z, h,'o', markersize=9,color='red',markeredgecolor='red', label='Inner', zorder=0)
        ax.plot(Z, h_r,'^', lw=2,color='blue', zorder=1)
      else:
        ax.plot(Z, h,'o', markersize=9,color='red',markeredgecolor='red', zorder=1)
    
    #if Rbreak!=0. and not math.isnan(Rbreak):
    #  plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')
    
    
    if regression_lims is not None:
        if len(input_files)==1:
            [minx, maxx, miny, maxy] = regression_lims[0]
            k,b,r_value = RegrLine(Z, h, minx* units[1], maxx* units[1], miny* units[1], maxy* units[1], 3)
            
            if b>=0.:
                print '\t Regression line: Y = %.3f * X + %.3f (r=%.2f)\n' % (k,b,r_value)
            else:
                print '\t Regression line: Y = %.3f * X - %.3f (r=%.2f)\n' % (k,abs(b),r_value)
        

            x_regr = np.arange(minx* units[1], maxx* units[1], (maxx-minx)* units[1]/100.)
            y_regr = k*x_regr + b
            ax.plot(x_regr, y_regr, '-.', lw=2, color='black', zorder=2)
        else:
            [minx, maxx, miny, maxy] = regression_lims[0]
            k,b,r_value = RegrLine(Z, h, minx* units[1], maxx* units[1], miny* units[1], maxy* units[1], 3)
            
            if b>=0.:
                print '\t Regression line: Y = %.3f * X + %.3f (r=%.2f)\n' % (k,b,r_value)
            else:
                print '\t Regression line: Y = %.3f * X - %.3f (r=%.2f)\n' % (k,abs(b),r_value)
        

            x_regr = np.arange(minx* units[1], maxx* units[1], (maxx-minx)* units[1]/100.)
            y_regr = k*x_regr + b
            ax.plot(x_regr, y_regr, '--', lw=2, color='red', zorder=2)




            [minx, maxx, miny, maxy] = regression_lims[1]
            k,b,r_value = RegrLine(Z_r, h_r, minx* units[1], maxx* units[1], miny* units[1], maxy* units[1], 3)
            
            if b>=0.:
                print '\t Regression line: Y = %.3f * X + %.3f (r=%.2f)\n' % (k,b,r_value)
            else:
                print '\t Regression line: Y = %.3f * X - %.3f (r=%.2f)\n' % (k,abs(b),r_value)
        

            x_regr = np.arange(minx* units[1], maxx* units[1], (maxx-minx)* units[1]/100.)
            y_regr = k*x_regr + b
            ax.plot(x_regr, y_regr, '-.', lw=2, color='blue', zorder=2)
            
            
            
   
    
    
    if len(input_files)==2:
        ax.set_ylim(0.5*min([min(h),min(h_r)]), 1.3*max([max(h),max(h_r)]))
        legend = ax.legend(loc=2, shadow=False, fontsize=10, numpoints=1)
    else:
        ax.set_ylim(0.5*min(h), 1.3*max(h))

    if xlim is not None:
        ax.set_xlim(xlim[0]* units[1], xlim[1]* units[1])

    if ylim is not None:
        ax.set_ylim(ylim[0]* units[1], ylim[1]* units[1])
                 
    ax.set_xlabel(r' $z$ (%s)' % (units[0]), fontsize=fsize)
    ax.set_ylabel(r' $h$ (%s) ' % (units[0]), fontsize=fsize)
    
    plt.savefig('h_Z.eps', transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.clf()
    plt.close()  







    '''
    #### 2) L0 Vs Z
    fig = plt.figure(1, figsize=(6, 5))
    ax=plt.subplot()
    L0 = np.log10(L0)
    
    if not np.isnan(L0_err[0]):
      if len(input_files)==2:
        ax.errorbar(Z, L0, yerr=L0_err,fmt='o',markersize=9,color='blue',markeredgecolor='blue', ecolor='blue', capthick=2, label='Left', zorder=0)
        ax.errorbar(Z_r, L0_r, yerr=L0_err_r,fmt='^',markersize=7,color='red',markeredgecolor='red', ecolor='red', capthick=2, label='Right', zorder=1)
      else:
        ax.errorbar(Z, L0, yerr=L0_err,fmt='o',markersize=9,color='red',markeredgecolor='red', ecolor='red', capthick=2, zorder=1)
    else:
      if len(input_files)==2:
        ax.plot(Z, L0,'o', markersize=9,color='blue',markeredgecolor='blue', label='Left', zorder=0)
        ax.plot(Z, L0_r,'^', markersize=7,color='red',markeredgecolor='red', label='Right', zorder=1)
      else:
        ax.plot(Z, L0,'o', markersize=9,color='red',markeredgecolor='red', zorder=1)
    
    #if Rbreak!=0. and not math.isnan(Rbreak):
    #  plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')
    

    if regression_lims is not None:
        [minx, maxx, miny, maxy] = regression_lims
        miny = -3.
        maxy = 1.
        k,b,r_value = RegrLine(Z, L0, minx* units[1], maxx* units[1], miny, maxy, 3)
        if b>=0.:
            print '\t Regression line: Y = %.3f * X + %.3f (r=%.2f)\n' % (k,b,r_value)
        else:
            print '\t Regression line: Y = %.3f * X - %.3f (r=%.2f)\n' % (k,abs(b),r_value)
    

        x_regr = np.arange(minx* units[1], maxx* units[1], (maxx-minx)* units[1]/100.)
        y_regr = k*x_regr + b
        ax.plot(x_regr, y_regr, '--', lw=3, color='black', zorder=2)

    
    
    if len(input_files)==2:
        ax.set_ylim(0.5*min([min(L0),min(L0_r)]), 1.3*max([max(L0),max(L0_r)]))
        legend = ax.legend(loc=1, shadow=False, fontsize=10, numpoints=1)
    else:
        ax.set_ylim(0.5*min(L0), 1.3*max(L0))

    if xlim is not None:
        ax.set_xlim(xlim[0]* units[1], xlim[1]* units[1])

    #if ylim is not None:
    #    ax.set_ylim(ylim[0], ylim[1])
    ax.set_ylim(-3., 1.) 
    #ax.set_yscale('log')
    ax.set_xlabel(r' $z$ (%s)' % (units[0]), fontsize=fsize)
    ax.set_ylabel(r' $\log L_0$ (DN)' , fontsize=fsize)
    
    plt.savefig('L0_Z.eps', transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.clf()
    plt.close()  
    '''

#input_files = ['radial_fits_1_single.dat']
#main(input_files, Ropt=None, units = ['pix',1.0]) 
    
    
    
    
    
    
    
    