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

def main(input_files, Ropt=None, mode='vertical', units = ['arcsec',0.75], regression_lims = None, xlim = None, ylim = None):
    # [minx, maxx, miny, maxy] = regression_lims - limits to plot regression line (should be in pix!)
    # First file should be left
    # Second file should be right
    R,L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err,L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err = np.loadtxt(input_files[0], usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],dtype=float, unpack=True,skiprows=1,delimiter='\t')
    R = R * units[1]
    Z0 = Z0 * units[1]
    Z0_err = Z0_err * units[1]
    Y0 = (Y0 - np.mean(Y0)) * units[1]
    Y0_err = Y0_err * units[1]
    
    Z01 = Z01 * units[1]
    Z01_err = Z01_err * units[1]
    Z02 = Z02 * units[1]
    Z02_err = Z02_err * units[1]

    
    if len(input_files)==2:
        R_r,L0_r,L0_err_r,Y0_r,Y0_err_r,N_r,N_err_r,Z0_r,Z0_err_r,L01_r,L01_err_r,N1_r,N1_err_r,Z01_r,Z01_err_r,L02_r,L02_err_r,N2_r,N2_err_r,Z02_r,Z02_err_r = np.loadtxt(input_files[1], usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],dtype=float, unpack=True,skiprows=1,delimiter='\t')        

        R_r = R_r * units[1]
        Z0_r = Z0_r * units[1]
        Z0_err_r = Z0_err_r * units[1]
        Y0_r = (Y0_r - np.mean(Y0_r)) * units[1]
        Y0_err_r = Y0_err_r * units[1]

        Z01_r = Z01_r * units[1]
        Z01_err_r = Z01_err_r * units[1]
        Z02_r = Z02_r * units[1]
        Z02_err_r = Z02_err_r * units[1]    
    
    
    
    
    ########### PLOTTING: SINGLE DISK ########### 
    
    #### 1) z0 Vs R
    fig = plt.figure(0, figsize=(6, 5))
    ax=plt.subplot()

    if not np.isnan(L0_err[0]):
      if len(input_files)==2:
        ax.errorbar(R, Z0, yerr=Z0_err,fmt='o',markersize=9,color='blue',markeredgecolor='blue', ecolor='black', capthick=2, label='Thin disc', zorder=0)
        ax.errorbar(R_r, Z0_r, yerr=Z0_err_r,fmt='^',markersize=7,color='red',markeredgecolor='red', ecolor='black', capthick=2, label='Thick disc', zorder=1)
      else:
        ax.errorbar(R, Z0, yerr=Z0_err,fmt='o',markersize=9,color='red',markeredgecolor='red', ecolor='red', capthick=2, zorder=1)
    else:
      if len(input_files)==2:
        ax.plot(R, Z0,'o', markersize=9,color='blue',markeredgecolor='blue', label='Thin disc', zorder=0)
        ax.plot(R_r, Z0_r,'^', markersize=7,color='red',markeredgecolor='red', label='Thick disc', zorder=1)
      else:
        ax.plot(R, Z0,'o', markersize=9,color='red',markeredgecolor='red', zorder=1)
    
    #if Rbreak!=0. and not math.isnan(Rbreak):
    #  plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')
    
    if len(input_files)==2:
        ax.set_ylim(0.5*min([min(Z0),min(Z0_r)]), 1.2*max([max(Z0),max(Z0_r)]))
        legend = ax.legend(loc=2, shadow=False, fontsize=10, numpoints=1)
    else:
        ax.set_ylim(0.5*min(Z0), 1.2*max(Z0))

    if regression_lims is not None:
        [minx, maxx, miny, maxy] = regression_lims
        k,b,r_value = RegrLine(R, Z0, minx* units[1], maxx* units[1], miny* units[1], maxy* units[1], 3)
        if b>=0.:
            print '\t Regression line: Y = %.3f * X + %.3f (r=%.2f)\n' % (k,b,r_value)
        else:
            print '\t Regression line: Y = %.3f * X - %.3f (r=%.2f)\n' % (k,abs(b),r_value)
    

        x_regr = np.arange(minx* units[1], maxx* units[1], (maxx-minx)* units[1]/100.)
        y_regr = k*x_regr + b
        ax.plot(x_regr, y_regr, '--', lw=3, color='black', zorder=2)
        
    
    if xlim is not None:
        ax.set_xlim(xlim[0]* units[1], xlim[1]* units[1])

    if ylim is not None:
        ax.set_ylim(ylim[0]* units[1], ylim[1]* units[1])
    
    #ax.axvline(x=146, ls='--')  # WARNING: NGC4302!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 
    ax.set_xlabel(r' $r$ (%s)' % (units[0]), fontsize=fsize)
    ax.set_ylabel(r' $z_0$ (%s) ' % (units[0]), fontsize=fsize)
    
    plt.savefig('z0_R.eps', transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.clf()
    plt.close()  
    


    #### 2) y0 Vs R
    fig = plt.figure(1, figsize=(6, 5))
    ax=plt.subplot()

    if not np.isnan(L0_err[0]):
      if len(input_files)==2:
        ax.errorbar(R, Y0, yerr=Y0_err,fmt='o',markersize=9,color='blue',markeredgecolor='blue', ecolor='black', capthick=2, label='Left')
        ax.errorbar(R_r, Y0_r, yerr=Y0_err_r,fmt='^',markersize=7,color='red',markeredgecolor='red', ecolor='black', capthick=2, label='Right')
      else:
        ax.errorbar(R, Y0, yerr=Y0_err,fmt='o',markersize=9,color='blue',markeredgecolor='blue', ecolor='black', capthick=2)
    else:
      if len(input_files)==2:
        ax.plot(R, Y0,'o', markersize=9,color='blue',markeredgecolor='blue', label='Left')
        ax.plot(R_r, Y0_r,'^', markersize=7,color='red',markeredgecolor='red', label='Right')
      else:
        ax.plot(R, Y0,'o', markersize=9,color='blue',markeredgecolor='blue')
    
    #if Rbreak!=0. and not math.isnan(Rbreak):
    #  plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')
    
    if len(input_files)==2:
        ax.set_ylim(0.5*min([min(Y0),max(Y0_r)]), 2.*max([max(Y0),max(Y0_r)]))
        legend = ax.legend(loc=1, shadow=False, fontsize=10, numpoints=1)
    else:
        ax.set_ylim(0.5*min(Y0), 2.*max(Y0))
                 
    ax.set_xlabel(r' $r$ (%s)' % (units[0]), fontsize=fsize)
    ax.set_ylabel(r' $y_0$ (%s) ' % (units[0]), fontsize=fsize)
    
    plt.savefig('y0_R.eps', transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.clf()
    plt.close()      


    

    #### 3) N Vs R
    fig = plt.figure(2, figsize=(6, 5))
    ax=plt.subplot()

    if not np.isnan(L0_err[0]):
      if len(input_files)==2:
        ax.errorbar(R, N, yerr=N_err,fmt='o',markersize=9,color='blue',markeredgecolor='blue', ecolor='black', capthick=2, label='Left')
        ax.errorbar(R_r, N_r, yerr=N_err_r,fmt='^',markersize=7,color='red',markeredgecolor='red', ecolor='black', capthick=2, label='Right')
      else:
        ax.errorbar(R, N, yerr=Z0_err,fmt='o',markersize=9,color='blue',markeredgecolor='blue', ecolor='black', capthick=2)
    else:
      if len(input_files)==2:
        ax.plot(R, N,'o', markersize=9,color='blue',markeredgecolor='blue', label='Left')
        ax.plot(R_r, N_r,'^', markersize=7,color='red',markeredgecolor='red', label='Right')
      else:
        ax.plot(R, N,'o', markersize=9,color='blue',markeredgecolor='blue')
    
    #if Rbreak!=0. and not math.isnan(Rbreak):
    #  plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')
    
    if len(input_files)==2:
        ax.set_ylim(0.1, 200)
        ax.set_yscale('log')
        legend = ax.legend(loc=1, shadow=False, fontsize=10, numpoints=1)
    else:
        ax.set_ylim(0.1, 200)
        ax.set_yscale('log')
                 
    ax.set_xlabel(r' $r$ (%s)' % (units[0]), fontsize=fsize)
    ax.set_ylabel(r' $n$ ', fontsize=fsize)
    
    plt.savefig('z0_n.eps', transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.clf()
    plt.close()     






    ########### PLOTTING: DOUBLE DISK ########### 
    
    #### 1) z0 Vs R
    fig = plt.figure(3, figsize=(6, 5))
    ax=plt.subplot()

    if not np.isnan(L0_err[0]):
      if len(input_files)==2:
        #### THIS IS IN USE
        ax.errorbar(R, Z01, yerr=Z01_err,fmt='^',markersize=9,color='blue',markeredgecolor='blue', ecolor='blue', capthick=2, label='Thin disc')
        #ax.plot(R_r, Z01_r,'--',lw=2,color='blue')#,markeredgecolor='blue', ecolor='blue', capthick=2, label='Right, thin')
        ax.errorbar(R, Z02, yerr=Z02_err,fmt='o',markersize=9,color='red',markeredgecolor='red', ecolor='red', capthick=2, label='Thick disc')
        #ax.plot(R_r, Z02_r,'--',lw=2,color='red')#,markeredgecolor='red', ecolor='red', capthick=2, label='Right, thick')                
      else:
        ax.errorbar(R, Z01, yerr=Z01_err,fmt='^',markersize=9,color='blue',markeredgecolor='blue', ecolor='blue', capthick=2, label='Thin disc')
        ax.errorbar(R, Z02, yerr=Z02_err,fmt='o',markersize=9,color='red',markeredgecolor='red', ecolor='red', capthick=2, label='Thick disc')        
    else:
      if len(input_files)==2:
        ax.plot(R, Z01,'^', markersize=9,color='blue',markeredgecolor='blue', label='Thin disc')
        ax.plot(R, Z02,'o', markersize=9,color='red',markeredgecolor='red', label='Thick disc')
        #ax.plot(R_r, Z01_r,'--', lw=2,color='blue',markeredgecolor='blue')#, label='Right, thin')
        #ax.plot(R_r, Z02_r,'--', lw=2,color='red',markeredgecolor='red')#, label='Right, thick')
      else:
        ax.plot(R, Z01,'^', markersize=9,color='blue',markeredgecolor='blue', label='Thin disc')
        ax.plot(R, Z02,'o', markersize=9,color='red',markeredgecolor='red', label='Thick disc')
    
    #if Rbreak!=0. and not math.isnan(Rbreak):
    #  plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')
    
    if len(input_files)==2:
        ax.set_ylim(0.5*min([min(Z01),min(Z02),min(Z01_r),min(Z02_r)]), 1.5*max([max(Z01),max(Z02),max(Z01_r),max(Z02_r)]))
        legend = ax.legend(loc=2, shadow=False, fontsize=10, numpoints=1)
    else:
        ax.set_ylim(0.5*min([min(Z01),min(Z02)]), 1.5*max([max(Z01),max(Z02)]))

    if regression_lims is not None:
        [minx, maxx, miny, maxy] = regression_lims
        k,b,r_value = RegrLine(R, Z02, minx* units[1], maxx* units[1], miny* units[1], maxy* units[1], 3)
        if b>=0.:
            print '\t Regression line: Y = %.3f * X + %.3f (r=%.2f)\n' % (k,b,r_value)
        else:
            print '\t Regression line: Y = %.3f * X - %.3f (r=%.2f)\n' % (k,abs(b),r_value)
    

        x_regr = np.arange(minx* units[1], maxx* units[1], (maxx-minx)* units[1]/100.)
        y_regr = k*x_regr + b
        ax.plot(x_regr, y_regr, '-.', lw=2, color='black', zorder=2)
        
        #'''
        # NGC3628
        #y_regr = 0.96*x_regr -463.5*units[1]   
        #ax.plot(x_regr, y_regr, '-.', lw=2, color='red', zorder=2)
        #'''
        
    #ax.set_ylim(0., 28.)
    if ylim is not None:
        ax.set_ylim(ylim[0]* units[1], ylim[1]* units[1])

    if xlim is not None:
        ax.set_xlim(xlim[0]* units[1], xlim[1]* units[1])
    
    ax.set_xlabel(r' $r$ (%s)' % (units[0]), fontsize=fsize)
    ax.set_ylabel(r' $z_0$ (%s) ' % (units[0]), fontsize=fsize)
    
    plt.savefig('z0_R_double.eps', transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.clf()
    plt.close()  
    

    #### 2) N Vs R
    fig = plt.figure(4, figsize=(6, 5))
    ax=plt.subplot()

    if not np.isnan(L0_err[0]):
      if len(input_files)==2:
        ax.errorbar(R, N1, yerr=N1_err,fmt='o',markersize=9,color='blue',markeredgecolor='blue', ecolor='black', capthick=2, label='Left, thin')
        ax.errorbar(R_r, N1_r, yerr=N1_err_r,fmt='^',markersize=7,color='blue',markeredgecolor='blue', ecolor='black', capthick=2, label='Right, thin')
        ax.errorbar(R, N2, yerr=N2_err,fmt='o',markersize=9,color='red',markeredgecolor='red', ecolor='black', capthick=2, label='Left, thick')
        ax.errorbar(R_r, N2_r, yerr=N2_err_r,fmt='^',markersize=7,color='red',markeredgecolor='red', ecolor='black', capthick=2, label='Right, thick')        
        
      else:
        ax.errorbar(R, N1, yerr=N1_err,fmt='o',markersize=9,color='blue',markeredgecolor='blue', ecolor='black', capthick=2, label='Thin')
        ax.errorbar(R, N2, yerr=N2_err,fmt='o',markersize=9,color='red',markeredgecolor='red', ecolor='black', capthick=2, label='Thick')        
    else:
      if len(input_files)==2:
        ax.plot(R, N1,'o', markersize=9,color='blue',markeredgecolor='blue', label='Left, thin')
        ax.plot(R, N2,'o', markersize=9,color='red',markeredgecolor='red', label='Left, thick')
        ax.plot(R_r, N1_r,'^', markersize=7,color='blue',markeredgecolor='blue', label='Right, thin')
        ax.plot(R_r, N2_r,'^', markersize=7,color='red',markeredgecolor='red', label='Right, thick')
      else:
        ax.plot(R, N1,'o', markersize=9,color='blue',markeredgecolor='blue', label='Thin')
        ax.plot(R, N2,'o', markersize=9,color='red',markeredgecolor='red', label='Thick')
    
    #if Rbreak!=0. and not math.isnan(Rbreak):
    #  plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')
    
    if len(input_files)==2:
        ax.set_ylim(0.1, 200)
        ax.set_yscale('log')
        legend = ax.legend(loc=1, shadow=False, fontsize=10, numpoints=1)
    else:
        ax.set_ylim(0.1, 200)
        ax.set_yscale('log')
                 
    ax.set_xlabel(r' $r$ (%s)' % (units[0]), fontsize=fsize)
    ax.set_ylabel(r' $n$', fontsize=fsize)
    
    plt.savefig('n_R_double.eps', transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.clf()
    plt.close()  


#input_files = ['vertical_fits_1_single.dat']
#input_files = ['vertical_fits_1_double.dat']
#main(input_files, Ropt=None, mode='vertical', units = ['pix',1.0])
    
    
    
    
    
    
    
    