#!/usr/bin/python
# Tested on /home/amosenko/Toshiba_1/CurrentWork/DustPedia_NEW/iraf_res/351_NGC4067

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
from mpl_toolkits.axes_grid1 import AxesGrid
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
import re
import glob
from matplotlib.colors import LogNorm
from astropy.io import fits as pyfits
import argparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import wcs
import warnings
warnings.filterwarnings("ignore")
from scipy.interpolate import splprep, splev
from scipy.interpolate import interp1d
from scipy import interpolate
import astropy
from scipy.optimize import least_squares

from sklearn import linear_model, datasets
from astropy.stats import sigma_clipped_stats
import collections
from matplotlib.path import Path
from itertools import product
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like


def show_complete(i,N):
  percent = 100. * float(i) / float(N)
  sys.stdout.write("\r%2d%%" % percent)
  sys.stdout.flush()

class PPoint:
    def __init__(self, x, y):
        self.x = x
        self.y = y    
        
def rot_point(p, orig, angle):
    """ performs rotation of a point (x,y) around point (x0,y0)
        by angle radians"""
    x1 = cos(angle)*(p.x-orig.x) - sin(angle)*(p.y-orig.y) + orig.x
    y1 = sin(angle)*(p.x-orig.x) + cos(angle)*(p.y-orig.y) + orig.y
    return PPoint(x1, y1)

# Main function, this part actually runs when routine is called
def ellipse_mask(cen, ellA, ellB, ellPA, inframe, xSize, ySize, img_inp=None,img_mask=None, outer = False):
        cospa = cos(radians(ellPA))
        sinpa = sin(radians(ellPA))
        # Check if whole ellipse is inside of the image
        # and obtain size of the ellipse in xy plane
        xMax = 0.0
        yMax = 0.0
        xMin = 1e10
        yMin = 1e10
        for e in linspace(0, 4*pi, 1000):
            cose = cos(e)
            sine = sin(e)
            x = cen.x + ellA * cose * cospa - ellB * sine * sinpa
            y = cen.y + ellB * sine * cospa + ellA * cose * sinpa
            if x > xMax:
                xMax = x
            if y > yMax:
                yMax = y
            if x < xMin:
                xMin = x
            if y < yMin:
                yMin = y
        xMin = max(0, int(round(xMin)))
        xMax = min(xSize, int(round(xMax)))
        yMin = max(0, int(round(yMin)))
        yMax = min(ySize, int(round(yMax)))
        focusR = (ellA ** 2.0 - ellB ** 2.0) ** 0.5
        focus10 = PPoint(cen.x + focusR, cen.y)  # Unrotated
        focus20 = PPoint(cen.x - focusR, cen.y)  #
        focus1 = rot_point(focus10, cen, radians(ellPA))
        focus2 = rot_point(focus20, cen, radians(ellPA))
        # Find pixels inside of the ellipse
        dEll = 2 * ellA
        
        if outer==True:
          xMin = 0
          xMax = xSize -1
          yMin = 0
          yMax = ySize -1

          if img_inp==None:
            for x in range(xMin, xMax+1):
              for y in range(yMin, yMax+1):
                    dFocus1 = hypot(x-focus1.x, y-focus1.y)
                    dFocus2 = hypot(x-focus2.x, y-focus2.y)
                    dPoint = dFocus1 + dFocus2
                    if dPoint >= dEll:
                      try:
                        inframe[y-1,x-1] = 0.0
                        img_mask[y-1,x-1] = 1
                      except:
                        zz=1
                    else:
                      try:
                        inframe[y-1,x-1] = inframe[y-1,x-1]
                        img_mask[y-1,x-1] = 0
                      except:
                        zz=1      
          else:
            for x in range(xMin, xMax+1):
              for y in range(yMin, yMax+1):
                    dFocus1 = hypot(x-focus1.x, y-focus1.y)
                    dFocus2 = hypot(x-focus2.x, y-focus2.y)
                    dPoint = dFocus1 + dFocus2
                    if dPoint >= dEll:
                      try:
                        inframe[y-1,x-1] = img_inp[y-1,x-1]
                      except:
                        zz=1

        else:
          if img_inp==None:
            for x in range(xMin, xMax+1):
              for y in range(yMin, yMax+1):
                    dFocus1 = hypot(x-focus1.x, y-focus1.y)
                    dFocus2 = hypot(x-focus2.x, y-focus2.y)
                    dPoint = dFocus1 + dFocus2
                    if dPoint < dEll:
                      try:
                        inframe[y-1,x-1] = 0.0
                        img_mask[y-1,x-1] = 1
                      except:
                        zz=1
          else:
            for x in range(xMin, xMax+1):
              for y in range(yMin, yMax+1):
                    dFocus1 = hypot(x-focus1.x, y-focus1.y)
                    dFocus2 = hypot(x-focus2.x, y-focus2.y)
                    dPoint = dFocus1 + dFocus2
                    if dPoint < dEll:
                      try:
                        inframe[y-1,x-1] = img_inp[y-1,x-1]
                      except:
                        zz=1
        return inframe
    
def grad_mag(radius, mag, max_radius):
    grad = []
    new_mag = []
    new_rad = []
    #max_radius = np.max(radius)
    for k in range(1,len(radius)-1):
      if radius[k]>max_radius/2.:
        g = (mag[k+1] - mag[k-1])/(2.*(radius[k+1]-radius[k-1]))
        if g<=0.:
            
            if new_mag==[]:
                grad.append( g )
                new_rad.append(radius[k])
                new_mag.append(mag[k])
            else:
                if abs(new_mag[-1]-mag[k])<1.:
                    grad.append( g )
                    new_rad.append(radius[k])
                    new_mag.append(mag[k])            


    return np.array(new_mag, float),np.array(grad, float), np.array(new_rad, float)

def func(k, x, y):
    return k[0]*x+k[1] - y

def fit_linear(x, y):
    k0 = np.ones(2)
    res_robust = least_squares(func, k0, loss='soft_l1', f_scale=0.1, args=(x, y))

    x_new = np.linspace(min(x), max(x), 300)
    y_new = func(res_robust.x, x_new, 0.)
    intercept = res_robust.x[1]
    return x_new,y_new,intercept


def fit_sklern(x, y):
    #print(x,y)
    # Robustly fit linear model with RANSAC algorithm
    ransac = linear_model.RANSACRegressor()
    ransac.fit(x, y)
    intercept = ransac.estimator_.intercept_ #coef_  
    
    inlier_mask = ransac.inlier_mask_
    outlier_mask = np.logical_not(inlier_mask)
    line_x = np.linspace(min(x), max(x), 300)[:, np.newaxis]
    try:
        line_y_ransac = ransac.predict(line_x)
    except:
        nsamples, nx, ny = line_x.shape
        line_x = line_x.reshape((nsamples,nx*ny))
        line_y_ransac = ransac.predict(line_x)
    return line_x,line_y_ransac,intercept[0]


def find_reff_simple(radius, mag, tot_mag):
    reff = float('nan')
    for k in range(len(radius)-1):
        if mag[k]>=tot_mag+2.5*log10(2.) and mag[k+1]<tot_mag+2.5*log10(2.):
            #reff = radius[k]
            r_l = radius[k]
            r_r = radius[k+1]
            mag_l = mag[k]
            mag_r = mag[k+1]
            break
    
    reff = (tot_mag+2.5*log10(2.)-mag_l)*(r_r-r_l)/(mag_r-mag_l) + r_l
    return reff
    


def find_eff_radius(radius, mag, tot_mag):
    reff = [float('nan')]
    s = sorted(zip(radius, mag))
    radius, mag = map(list, zip(*s))

    radius = np.array(radius)
    mag = np.array(mag)
    #plt.plot(radius, mag, 'o')
    #plt.show()
    #exit()
    again = False
    try:
        yToFind = tot_mag+2.5*log10(2.)
        yreduced = np.array(mag) - yToFind
        freduced = interpolate.UnivariateSpline(radius, yreduced, s=0)
        reff = freduced.roots()    
    except:
        again=True
    
    if list(reff)==[]:
        again=True

    if again==True:
        new_mag = []; new_radius = []
        new_mag = [mag[0]]
        new_radius = [radius[0]]
        for k in range(1,len(radius)):
            if abs(new_mag[-1]-mag[k])<1.:
                new_mag.append(mag[k])
                new_radius.append(radius[k])
        
        
        radius = np.array(new_radius)
        mag = np.array(new_mag)


    
        yToFind = tot_mag+2.5*log10(2.)
        yreduced = np.array(mag) - yToFind

        try:
            freduced = interpolate.UnivariateSpline(radius, yreduced, s=0)
            reff = freduced.roots()     
            #print reff
        except:
            try:
                reff = [find_reff_simple(radius, mag, tot_mag)]
            except:
                reff = []
    #plt.show()
    #if list(reff) == []:
    #        reff = [find_reff_simple(radius, mag, tot_mag)]       
    #print find_reff_simple(radius, mag, tot_mag)        
    if list(reff)==[]:
        reff = [float('nan')]
    #print reff
    return reff[0]    



def find_sigma_instr(image_data, N_ap, sigma_data=None):
    if sigma_data is None:
        #print(image_data)
        #exit()
        mean, median, stddev = sigma_clipped_stats(image_data, mask=None, mask_value=0.0, sigma=10.0, sigma_lower=None, sigma_upper=None, maxiters=1)
        
        sigma_poisson = stddev / math.sqrt(N_ap)
        sigma_instr = N_ap*sigma_poisson
    else:        
        sigma_instr = math.sqrt( np.sum(sigma_data**2) )
    
    return sigma_instr

def find_sigma_sky(rms_sky, N_ap):
    return rms_sky/math.sqrt(N_ap)#*N_ap#/math.sqrt(N_ap) ###TODO:????


def find_error(galaxy_apperture, tot_mag, m0, input_image, rms_sky, sigma_image=None, sigma_cal=0., sigma_cal_units='%%', output_units='mag', sigma_apert=0., sigma_apert_units='mag'):
    # galaxy_apperture = [x0, y0, sma, ell, PA]
    # tot_mag - total magnitude, in mags
    # m0 - Zero point in mags
    # input_image - fits file 
    # rms_sky - rms in ADU
    # sigma_image - fits file
    # sigma_cal - callibration error
    # sigma_cal_units - units of callibration error (ADU or %) 
    # output_units - mag, ADU or %
    
    hdulist_i = pyfits.open(input_image)
    image_data = hdulist_i[0].data


    [x0, y0, sma, ell, PA] = galaxy_apperture
    #sma = sma/3.
    PA = PA + 90.
    cen = PPoint(x0,y0)
    ySize, xSize = np.shape(image_data)

    image_data_masked = ellipse_mask(cen, sma, sma*(1.-ell), PA, image_data, xSize, ySize, img_inp=None,img_mask=None, outer = True)
    
    #hdu = pyfits.PrimaryHDU(image_data_masked)
    #hdu.writeto('tmp.fits', clobber=True)    
    #exit()
    image_data = image_data_masked

    #print(image_data)
    #exit()    
    if sigma_image is not None:
        hdulist_s = pyfits.open(sigma_image)
        sigma_data = hdulist_s[0].data
        sigma_data_masked = ellipse_mask(cen, sma, sma*(1.-ell), PA, sigma_data, xSize, ySize, img_inp=None,img_mask=None, outer = True)
        sigma_data = sigma_data_masked
    else:
        sigma_data = None

    N_ap = np.count_nonzero(image_data)
    sigma_instr = find_sigma_instr(image_data, N_ap, sigma_data=sigma_data)
    sigma_sky = find_sigma_sky(rms_sky, N_ap)
    
    if sigma_cal_units=='%%':
        sigma_cal = 10**(0.4*(m0-tot_mag)) * sigma_cal/100.
    
    if sigma_cal_units=='mag':
        sigma_tot = math.sqrt( sigma_instr**2 + sigma_sky**2) 
        sigma_tot_m = (2.5*log10(math.e) / ( 10**(0.4*(m0-tot_mag)) ) ) * (sigma_tot) + sigma_apert + sigma_cal
        
    if sigma_cal_units=='ADU' or sigma_cal_units=='%%':
        sigma_tot = math.sqrt( sigma_instr**2 + sigma_sky**2 + sigma_cal**2 )    
        sigma_tot_m = (2.5*log10(math.e) / ( 10**(0.4*(m0-tot_mag)) ) ) * (sigma_tot) + sigma_apert
    
    sigma_tot_i = sigma_tot_m*(10**(0.4*(m0-tot_mag))) / (2.5*log10(math.e))
    

    if output_units=='mag':
        return sigma_tot_m
    elif output_units=='%%':
        return sigma_tot_i*100./(10**(0.4*(m0-tot_mag)))
    else:
        return sigma_tot_i

def smoothTriangle(data, data_err, degree, dropVals=False):
    #https://plot.ly/python/smoothing/
    triangle=np.array(list(range(degree)) + [degree] + list(range(degree)[::-1])) + 1
    
    smoothed=[]; smoothed_err=[]

    for i in range(degree, len(data) - degree * 2):
        point=data[i:i + len(triangle)] * triangle
        point_err=data_err[i:i + len(triangle)] * triangle
        smoothed.append(sum(point)/sum(triangle))
        smoothed_err.append(sum(point_err)/sum(triangle))
    if dropVals:
        return smoothed, smoothed_err
    smoothed=[smoothed[0]]*int(degree + degree/2) + smoothed
    smoothed_err=[smoothed_err[0]]*int(degree + degree/2) + smoothed_err
    while len(smoothed) < len(data):
        smoothed.append(smoothed[-1])

    while len(smoothed_err) < len(data):
        smoothed_err.append(smoothed_err[-1])
    
    return smoothed,smoothed_err

def create_average_profile(r, I, I_err, bin_box=10):
    if bin_box>1:
        #I = smooth(I,bin_box)
        I,I_err = smoothTriangle(I, I_err, bin_box)
        #I = regr_fit(r,I)
    return r, I, I_err

def calc_scale(D):
    return D*1000./206265. # kpc per arcsec

def calc_Mag(mag, Dist, Aext):
        return mag - 5.*np.log10(Dist) - Aext - 25.


def get_mag_by_curve(r, mag, R_lim):
    f = interpolate.interp1d(r, mag)
    m_lim = f(R_lim)
    return float(m_lim)


def create_growth_curve(a, mag, pix2sec, output_file, R_lim, mag_lim, text):
        fig = figure(2,figsize=(5,5))
        ax = fig.add_subplot(111)
    
        fsize=16

        if True:
            color = 'gray'
            color_bar = 'black'
            symb = 'o'


        ax.plot(a*pix2sec, mag, symb,color=color,markeredgecolor=color, markersize=3)

        ax.axhline(y=mag_lim, ls='--')
        ax.axvline(x=R_lim*pix2sec, ls='-.')


        ax.set_ylim(max(mag),min(mag)-1.)
        ax.set_xlim(0.,max(a)*pix2sec)


        ax.set_xlabel(r'$r$ (arcsec)', fontsize=fsize)
        ax.set_ylabel(r'$mag$', fontsize=fsize)
        if text is not None:
            ax.text(0.7, 0.9, text, fontsize=fsize, color='black',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline')

        fig.savefig(output_file, transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
        plt.clf()
        plt.close()        

def find_min_radius(func, r, SB_lim, m0,pix2sec,rmin=10):

        for R in arange(rmin,max(r),0.1):
            #print(R)
            if m0 - 2.5*log10(func(R))+ 5.*log10(pix2sec)<SB_lim:
                Rmax = R
            else:
                break
        return Rmax


def create_profile_with_averaging(a, inten, inten_err, a_av, inten_av, m0, pix2sec, R_lim, R_lim_err, R_min, mu_min, output_file, text, SB_lim):
        fig = figure(1,figsize=(5,5))
        ax = fig.add_subplot(111)
    
        mag = m0 - 2.5*log10(inten)+ 5.*log10(pix2sec)
        mag_err = fabs((2.5/log(10.0)) * inten_err/inten)
        
        mag_av = m0 - 2.5*log10(inten_av)+ 5.*log10(pix2sec)
        ymin = None
        ymax = None
        rmax = max(a)
        fsize=16

        if True:
          color = 'gray'
          color_bar = 'black'
          symb = 'o'

        '''
        set_interp = sp.interpolate.PchipInterpolator(a*pix2sec, mag, axis=0, extrapolate=True)
        x_new = []
        y_new = []
        for k in range(len(a)):
            ss = set_interp(a[k]*pix2sec)
            if not np.isnan(ss):
                y_new.append(ss)
                x_new.append(a[k]*pix2sec)
        '''
        
        
        '''
        # Overplot PSF fo NGC2903: WARNING -use this to create Fig. in paper
        file_ext_psf = '/home/amosenko/MyCurrentWork/Edge_on_HERON/PSF_HERON/NEW_PSF/psf_extended_for_paper.txt'
        sma_psf, inten_psf = np.loadtxt(file_ext_psf, usecols=[0,1], unpack=True, skiprows = 1, dtype=float)
        ax.plot(sma_psf*pix2sec,  m0 - 2.5*log10(inten_psf*inten[0]/inten_psf[0])+ 5.*log10(pix2sec), '^', color='lime', markersize=5, zorder=0)
        '''
        
        
        
        ax.errorbar(a*pix2sec,mag,yerr=mag_err,fmt=symb,color=color,ecolor=color, markersize=5, zorder=1)
        ax.plot(a*pix2sec, mag, symb,color=color,markeredgecolor=color, markersize=3, zorder=2)

        #ax.plot(x_new, y_new, color='red', lw=2)

        #ax.plot(a_av*pix2sec, mag_av, color='red', lw=2, zorder=3) #### COMMENTED!!!!!!!!!!!!
        ax.axhline(y=SB_lim, ls='-', zorder=4)
        ax.axvline(x=R_lim*pix2sec, ls='-', zorder=5)
        ax.axvspan((R_lim-R_lim_err)*pix2sec, (R_lim+R_lim_err)*pix2sec, alpha=0.5, color='blue', zorder=6)

        ax.axhline(y=mu_min, ls='-.', lw=1, color='red', zorder=7)
        ax.axvline(x=R_min*pix2sec, ls='-.', lw=1, color='red', zorder=8)
        
        mean_mag, median_mag, std_mag = sigma_clipped_stats(mag, sigma=3.0)

        if ymin is None:  
            #ax.set_ylim(max(mag)+std_mag,min(mag)-std_mag)
            ax.set_ylim(32.,min(mag)-std_mag)
        else:
            ax.set_ylim(ymax,ymin)
        ax.set_xlim(0.,1.02*max([R_lim,R_min])*pix2sec)


        ax.set_xlabel(r'$r$ (arcsec)', fontsize=fsize)
        ax.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
        if text is not None:
           ax.text(0.5, 0.9, text, fontsize=fsize, color='black',transform=ax.transAxes, horizontalalignment='center',verticalalignment='baseline')
        
        ax.get_xaxis().set_tick_params(direction='in', width=1)
        ax.get_yaxis().set_tick_params(direction='in', width=1)        
        #plt.show()
        #exit()
        plt.savefig(output_file, transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
        plt.clf()
        plt.close()   

def find_minSB_func(r, f_int, f_int_err, sigma=2., rmin=10):
        for R in arange(rmin, max(r), 0.1):
            if f_int(R)<sigma*f_int_err(R):
                min_SB=f_int(R)
                Rmin = R
                break
            else:
                min_SB=f_int(R)
                Rmin = R                

        return min_SB,Rmin

def simulate_profile(sma, inten, sigma_sum, sky_med_unc):
    new_inten = []
    for k in range(len(sma)):
        new_inten.append(random.gauss(inten[k]+random.gauss(0., sky_med_unc), sigma_sum[k]))
    
    return new_inten

def find_min_ellipticity(sma, ell, R_lim):
    ELL = []
    for i in range(len(sma)):
        if sma[i]>R_lim-5. and sma[i]<R_lim+5.:
            ELL.append(ell[i])
    lim_ell = np.median(ELL)        

    return lim_ell


def curve_of_growth(ellipse_file, m0, pix2sec=1, bin_box=5, SB_lim=28., rmin=10., sky_std=0., sigma_cal_mag=None, sky_med_unc=0., Aext=0., find_Rlim_err=False):
              sky_std = sky_std * 10**(0.4*Aext)
              sky_med_unc = sky_med_unc * 10**(0.4*Aext)
              
              try:
                    sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4,TFLUX_E,GRAD,NDATA = loadtxt(ellipse_file, usecols=[1,2,3,6,7,8,9,10,12,33,34,21,15,35], unpack=True, skiprows = 6, dtype='str')
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
                                    if TFLUX_E[k]=='INDEF': TFLUX_E[k]=0
                                    if GRAD[k]=='INDEF': GRAD[k]=0
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
                    TFLUX_E = np.array(TFLUX_E,dtype='float')
                    GRAD_I = np.array(GRAD,dtype='float')
                    NDATA = np.array(NDATA,dtype='float')
                    
                    mag = m0 - 2.5*log10(TFLUX_E) - Aext
                    inten = inten * 10**(0.4*Aext)
                    inten_err = inten_err * 10**(0.4*Aext)
              except:
                    # Read azimuthally-averaged profile:
                    sma,inten,inten_err = np.loadtxt(ellipse_file, usecols=[0,1,2], dtype=float, unpack=True, skiprows = 1, delimiter='\t')
                    inten = inten * 10**(0.4*Aext)
                    inten_err = inten_err * 10**(0.4*Aext)
                    
                    # Convert SB profile to the growth curve
                    mag = []
                    smaa = []
                    intenn = []
                    intenn_err = []
                    
                    for kk in range(len(sma)):
                        yy = np.array(sma[0:kk+1])*np.array(inten[0:kk+1])
                        L = 2.*math.pi*np.trapz(yy, dx = 1.)*(1.-ell)
                        if L!=0.:
                            mag.append( m0 - 2.5*log10(L) )
                            intenn.append(inten[kk])
                            intenn_err.append(inten_err[kk])
                    mag = np.array(mag)
                    sma = np.array(smaa)
                    inten = np.array(intenn)
                    inten_err = np.array(intenn_err)

              # Calculate full errors:
              sigma_instr = inten_err
              sigma_sky = sky_std/sqrt(NDATA)
              
              if sigma_cal_mag is not None:
                 sigma_cal = sigma_cal_mag*inten / (2.5*log10(math.e))
              else:
                 sigma_cal = 0. 

              sigma_sum = sqrt( sigma_instr**2 + sigma_sky**2 + sigma_cal**2 )
              
                  
              # Averaging profile:
              sma_av, inten_av, inten_av_err = create_average_profile(sma, inten, sigma_sum, bin_box=bin_box)              
              
              
              # Interpolate averaged profile:
              f_int = interpolate.interp1d(sma_av, inten_av)
              f_int_err = interpolate.interp1d(sma_av, inten_av_err)              
              
              R_lim = find_min_radius(f_int, sma_av, SB_lim, m0, pix2sec, rmin)
              ell_lim = find_min_ellipticity(sma, ell, R_lim)
              
              #'''
              # Monte-Carlo modelling to estimate the error on R_lim:
              R_lim_gauss = []
              if find_Rlim_err:
                    print('Monte Carlo simulations to find error on Rlim...')
                    for k in range(1000):
                            show_complete(k+1, 1000) # Progress status
                            inten_gauss = simulate_profile(sma, inten, sigma_sum, sky_med_unc)

                            # Averaging profile:
                            sma_av_gauss, inten_av_gauss, inten_av_gauss_err = create_average_profile(sma, inten_gauss, sigma_sum, bin_box=bin_box)              
                            
                            
                            # Interpolate averaged profile:
                            f_int_gauss = interpolate.interp1d(sma_av_gauss, inten_av_gauss)
                            #f_int_gauss_err = interpolate.interp1d(sma_av_gauss, inten_av_gauss_err)              
                            try:
                                R_lim_gauss.append(find_min_radius(f_int_gauss, sma_av_gauss, SB_lim, m0, pix2sec, rmin)) 
                            except:
                                z=1
                    print('\n')
                    R_lim_err = np.std(R_lim_gauss)
              else:
                  R_lim_err = float('nan')
              #'''
              #R_lim_err = 18.96 ###TODO: WARNING!!!! This is for NGC2903
              
              # Find magnitude mag_lim at the R_lim
              mag_lim = get_mag_by_curve(np.array(sma), mag, R_lim)

              # Find minimal SB level (where I < sigma*Ierr) and corresponding radius:
              I_min, R_min = find_minSB_func(sma_av, f_int, f_int_err, sigma=2., rmin=rmin)
              mu_min = m0 - 2.5*math.log10(I_min/(pix2sec**2))


              # Plot curve of growth:
              create_growth_curve(np.array(sma), mag, pix2sec, 'growth_curve.eps', R_lim, mag_lim , text=None)
      
              mag_grad,grad,rad = grad_mag(sma, mag, R_lim/2.)
              #plt.plot(grad, np.log10(rad),'o')

              grad = np.array(grad)
              mag_grad = np.array(mag_grad)
              #plt.plot(grad, mag_grad,'o')
              #plt.show()
              #exit()
              grad_new_new,mag_new_new,mag_asympt = fit_sklern(grad.reshape((len(grad),1)), mag_grad.reshape((len(mag_grad),1)))
              
              grad_new_new_new,rad_new_new_new,R_asympt = fit_sklern(grad.reshape((len(grad),1)), np.log10(rad).reshape((len(rad),1)))
              R_asympt = 10**(R_asympt)
              #plt.plot(grad_new_new_new, rad_new_new_new)
              #plt.show()

              
              eff_radius = find_eff_radius(sma, mag, mag_asympt)
              
              
              x0 = np.median(x0)
              y0 = np.median(y0)
              ell = np.median(ell)
              PA = np.median(PA)

              
              # Plot profiles:
              create_profile_with_averaging(sma, inten, sigma_sum, sma_av, inten_av, m0, pix2sec, R_lim, R_lim_err, R_min, mu_min, 'azim_aver_%.1f.png' % (SB_lim), '', SB_lim)
              
              return R_min, mu_min, R_lim, R_lim_err, mag_lim, R_asympt, mag_asympt, eff_radius, [x0,y0,ell,PA], ell_lim






def convert_AB_to_Jy(mag):
    return 10**(0.4*(8.90-mag))

def convert_Ierr_to_magerr():
    z=1

def convert_magerr_to_Ierr(tot_mag_I, magerr):
    return magerr*(tot_mag_I / (2.5*log10(math.e)))

def main(input_image, ellipse_file, m0, pix2sec, rms_sky=0., apert_cor=0., sigma_image=None, sigma_cal=0., sigma_apert=0., bin_box=5., SB_lim=28., rmin=10., sky_med_unc=0., Aext=0.):
    print('Estimating general parameters...')
    res = collections.OrderedDict()

    R_min, mu_min, R_lim, R_lim_err, mag_lim, R_asympt, mag_asympt, eff_radius, [x0,y0,ell,PA], ell_lim = curve_of_growth(ellipse_file, m0, pix2sec=pix2sec, bin_box=bin_box, SB_lim=SB_lim, rmin=rmin, sky_std=rms_sky, sigma_cal_mag=sigma_cal, sky_med_unc=sky_med_unc, Aext=Aext)
    
    mag_lim = mag_lim - apert_cor
    mag_asympt = mag_asympt - apert_cor
    
    
    res['min_radius'] = R_min
    res['min_SB'] = mu_min   
    res['lim_radius'] = R_lim
    res['lim_radius_err'] = R_lim_err    
    res['lim_mag'] = mag_lim
    res['asympt_radius'] = R_asympt
    res['asympt_mag'] = mag_asympt
    res['eff_radius'] = eff_radius
    res['galaxy_aperture'] = [x0,y0,ell,PA]
    res['lim_ell'] = ell_lim
    
    # Estimate error of the mag_lim:
    galaxy_apperture = [x0, y0, R_lim, ell, PA]
    mag_lim_err = find_error(galaxy_apperture, mag_lim, m0, input_image, rms_sky, sigma_image=sigma_image, sigma_cal=sigma_cal, sigma_cal_units='mag', output_units='mag', sigma_apert=0., sigma_apert_units='mag')

    res['lim_mag_err'] = mag_lim_err

    galaxy_apperture = [x0, y0, R_asympt, ell, PA]
    mag_asympt_err = find_error(galaxy_apperture, mag_asympt, m0, input_image, rms_sky, sigma_image=sigma_image, sigma_cal=sigma_cal, sigma_cal_units='mag', output_units='mag', sigma_apert=0., sigma_apert_units='mag')

    res['asympt_mag_err'] = mag_asympt_err
    
    print(res['lim_mag'],res['lim_mag_err'])
    print(res['asympt_mag'],res['asympt_mag_err'])
    return res


# NGC509
#main('galaxy_r.fits', 'ellipse.txt', 28.2487, 0.396, rms_sky=6.146, apert_cor=0., sigma_image=None, sigma_cal=0.01, sigma_apert=0., bin_box=5, SB_lim=26., rmin=10., sky_med_unc=0.4672, Aext=0.106) 


# NGC4452
#main('galaxy_r.fits', 'ellipse.txt', 28.1883, 0.396, rms_sky=4.407, apert_cor=0., sigma_image=None, sigma_cal=0.01, sigma_apert=0., bin_box=5, SB_lim=26., rmin=10., sky_med_unc=0.5345, Aext=0.068) 
#main('composed_model.fits', 'ellipse.txt', 28.098, 0.396, rms_sky=2.06, apert_cor=0., sigma_image=None, sigma_cal=0.01, sigma_apert=0., bin_box=5, SB_lim=25., rmin=10., sky_med_unc=0.5345, Aext=0.068) 
#main('composed_model.fits', 'ellipse.txt', 28.098, 0.396, rms_sky=2.06, apert_cor=0., sigma_image=None, sigma_cal=0.01, sigma_apert=0., bin_box=5, SB_lim=28., rmin=10., sky_med_unc=0.5345, Aext=0.068) 

# NGC4469
#main('galaxy_r.fits', 'ellipse.txt', 28.2838, 0.396, rms_sky=4.9524, apert_cor=0., sigma_image=None, sigma_cal=0.01, sigma_apert=0., bin_box=5, SB_lim=26., rmin=10., sky_med_unc=0.397, Aext=0.046) 
#main('combined_fixed.fits', 'ellipse.txt', 28.294, 0.396, rms_sky=2.06, apert_cor=0., sigma_image=None, sigma_cal=0.01, sigma_apert=0., bin_box=5, SB_lim=25., rmin=10., sky_med_unc=0.397, Aext=0.046)
#main('combined_fixed.fits', 'ellipse.txt', 28.294, 0.396, rms_sky=2.06, apert_cor=0., sigma_image=None, sigma_cal=0.01, sigma_apert=0., bin_box=5, SB_lim=28., rmin=10., sky_med_unc=0.397, Aext=0.046) 
 
'''' 
def calc_apertures(input_image, outer_ellipse, step):
    [xc,yc,sma_max,smb_max,pa] = outer_ellipse
    q = smb_max/sma_max
    sma = 1.
    smb = sma*q
    fluxes = []
    radii = []
    N = 1
    while sma<=sma_max:
        print(N)
        ellipse = "ellipse(%f,%f,%f,%f,%f)" % (xc,yc,sma,smb,pa)
        command = "~/MEGA/MOL_A/Pipeline_scripts/funtools-master/funcnts tmp.fits \"%s\" > res.txt" % (ellipse)
        subprocess.call(command, shell=True)
        sma = sma + step
        smb = q*sma
        
        f = open('res.txt', 'r')
        lines = f.readlines()
        for k in range(len(lines)):
            if '# source_data' in lines[k]:
                tot_flux = float(lines[k+3].split()[1])
                fluxes.append(tot_flux)
                radii.append(sma)
        os.remove('res.txt')
        N = N + 1
'''







    
'''
    if Distance is not None:
        Scale = calc_scale(Distance)
    else:
        Scale = float('nan')


    R_lim_kpc = R_lim * pix2sec * Scale
    R_asympt_kpc = R_asympt_kpc * pix2sec * Scale
    eff_radius_kpc = eff_radius * pix2sec * Scale
    
    
    Mag_lim = calc_Mag(mag_lim, float(distance[k]), A_r)


              # Absolute magnitude:
              #Mag_lim = calc_Mag(mag_lim, float(distance[k]), A_r)



              #R_lim_kpc = R_lim*pix2sec[k]*calc_scale(float(distance[k]))


    sigma_tot = find_error(galaxy_apperture, tot_mag, m0, input_image, rms_sky, sigma_image=sigma_image, sigma_cal=sigma_cal, sigma_cal_units='%%', output_units='mag', sigma_apert=sigma_apert, sigma_apert_units='mag')
    res['apparent_mag_err'] = sigma_tot
    
    #mags = []
    #mags.append(tot_mag)
    tot_mag = tot_mag + to_AB
    res['apparent_mag_AB'] = tot_mag    
    
    #mags.append(tot_mag)    
    tot_mag = tot_mag - apert_cor
    res['apparent_mag_AB_apert_cor'] = tot_mag      
    
    
    #mags.append(tot_mag)    
    tot_mag = tot_mag - apert_cor - Aext
    #mags.append(tot_mag)   
    res['apparent_mag_AB_apert_cor_ext'] = tot_mag        
    
    tot_mag = tot_mag - apert_cor - Kcorr
    #mags.append(tot_mag)  
    res['apparent_mag_AB_apert_cor_ext_Kcorr'] = tot_mag     
    
    tot_mag_Jy = convert_AB_to_Jy(tot_mag)
    res['apparent_mag_AB_apert_cor_ext_Kcorr_Jy'] = tot_mag_Jy   


    tot_mag_Jy_err =convert_magerr_to_Ierr(tot_mag_Jy, sigma_tot)
    res['apparent_mag_AB_apert_cor_ext_Kcorr_Jy_err'] = tot_mag_Jy_err   

    #mags.append(tot_mag)    
    tot_mag = tot_mag - 5.*(log10(Distance*10**6)-1.)
    res['abs_mag'] = tot_mag   

    #mags.append(tot_mag)       
    

    #print res
'''










'''
outer_ellipse = [413.6651,341.65102,78.949773,40.424773,50]
input_image = '/home/amosenko/Toshiba_1/CurrentWork/DustPedia_NEW/results_sersic/350_NGC4051/W1/galaxy_clean_galf.fits'
step = 1.
calc_apertures(input_image, outer_ellipse, step)
exit()        
 
input_image = '/home/amosenko/Toshiba_1/CurrentWork/DustPedia_NEW/results_sersic/350_NGC4051/W1/galaxy_clean_galf.fits'
rms_sky = 0.059
sigma_image = '/home/amosenko/Toshiba_1/CurrentWork/DustPedia_NEW/results_sersic/350_NGC4051/W1/sigma.fits'
sigma_cal = 2.9
sigma_apert = 0.#-0.034
apert_cor = -0.034
Aext = 0.
Kcorr = 0.
to_AB = 2.699
main('ellipse.txt', 20.5, 1.,  Aext, Kcorr, apert_cor, input_image, rms_sky, sigma_image, sigma_cal, sigma_apert, to_AB)
'''