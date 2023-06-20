#!/usr/bin/python
# Module to create initial guess for edge-on galaxies
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **        © Astronomical Observatory, Ghent University         **
# *****************************************************************

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
import warnings
from scipy import special
import argparse
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import scipy.optimize as opt
from scipy.integrate import quad
from scipy.odr.odrpack import *
from scipy import optimize

# Import astro packages
import pyfits

#*** Colour fonts ***
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''


def help_func():
      '''
      The script to find bulge and disc domination radii and estimate breaks (truncations) of the disc (rough double disc fitting).
      INPUT:
      input_image - image with the centered edge-on galaxy oriented horizontally. However, you can optionally specify the center of the galaxy.
      mask_image - image with the masks where masked objects have intenisty >0. This mask will be taken into account during the analysis.
      mag_level - level in [mag/arcsec^2] of the outer isophote. All pixels with m>mag_level will be ignored!
      Radius - radius of the outer isophote.
      zmin - all pixels in |z|<=zmin will not be taken into account.
      m0 - magnitude Zero Point in [mag/arcsec^2].
      pix2sec - pixel size in [arcsec/pix].
	Optional:
	xc,zc - center of the galaxy. Otherwise the center of the frame will be used.
	
      OUTPUT:
      The disctionary with the mean values (averaged for the left and right side of the galaxy):
      PARS["FIRST_DISC"] = [h [pix],m0d [mag/arcsec^2]]
      PARS["SECOND_DISC"] = [h [pix],m0d [mag/arcsec^2]]
      PARS["BREAK_RADIUS"] = [R_left [pix],R_right [pix]]
      PARS["BULGE_DOM"] = [Bulge dominated radius along x-axis [pix],Bulge dominated radius along z-axis [pix]]
      For the single disc:
      PARS["DISC_DOM"] = [h [pix],m0d [mag/arcsec^2], Disc dominated radius [pix]]
      '''
      print help_func.__doc__

def check_sign(lst):
    if lst[0] > 0.:
	  dif_sign = 1
	  for k in range(1,len(lst),1):
	    if lst[k]<0.:
	      dif_sign = 0
	  if dif_sign == 1:
	    return True
	  else:
	    return False
    else:
	  dif_sign = 1
	  for k in range(1,len(lst),1):
	    if lst[k]>0.:
	      dif_sign = 0
	  if dif_sign == 1:
	    return True
	  else:
	    return False

def disk_edge_soph(z,I0d,n,z0,z_c):
	      #*** For edge-on disk SB in mag/arcsec^2 (along z-axis). n - index of this law. ***
	      # B[0] = I0d
	      # B[1] = n
	      # B[2] = z0 
	      # B[3] = z_c
	      
	      return I0d * (1.0 / np.cosh(n*fabs(z-z_c)/(2.*z0)) ) **(2./n)
	    
def disk_edge(z,I0d,z0,z_c):
	      #*** For edge-on disk SB in mag/arcsec^2 (along z-axis).  Sech^2-law. ***
	      # B[0] = I0d
	      # B[1] = z0 
	      # B[2] = z_c

	      return I0d * (1.0 / np.cosh(2.*fabs(z-z_c)/(2.*z0)) ) **(2./2.)
	    
def disk_edge_exp(z,I0d,z0,z_c):
	      #*** For edge-on disk SB in mag/arcsec^2 (along z-axis). n - index of this law. ***
	      # B[0] = I0d
	      # B[2] = z0 
	      # B[3] = z_c
	      
	      return I0d * exp(-fabs(z-z_c)/(2.*z0))

def chi2_func(ref,mod):
  return sum( (ref-mod)**2)

def disk_edge_r(B,x):
	#*** For edge-on disk SB in mag/arcsec^2 (along r-axis). ***
	# B[0] = I0d
	# B[1] = h
	# B[2] = x_c
	x = np.array(x)
	return B[0] * fabs(x-B[2])/B[1] * sp.special.kn(1,fabs(x-B[2])/B[1])  

def fit_eon_radial_prof(r, Int, I0d, h):
	def disk_edge_r(B,x):
		#*** For edge-on disk SB in mag/arcsec^2 (along r-axis). ***
		# B[0] = I0d
		# B[1] = h
		# B[2] = x_c
		x = np.array(x)
		return B[0] * fabs(x-B[2])/B[1] * sp.special.kn(1,fabs(x-B[2])/B[1])  
      	dataToFit = RealData(r, Int)
	mo = Model(disk_edge_r)
	fitt = ODR(dataToFit, mo, [I0d/1.9,h/1.3,0.])
	fitt.set_job()
	fit_result = fitt.run()
	H = fit_result.beta[1]
	I0D = fit_result.beta[0]
	'''
	plt.clf()
	plt.close()
	plt.plot(r,Int,'o',r,disk_edge_r([I0D,H,fit_result.beta[2]],r))
	'''
	return H,I0D

def dust_lane_extrema(I,z):
	    #INPUT:
	    #I [DN], z[pix] - centered
	    
	    # 0) Smooth the data
	    #len(z)
	    I_sav = savgol_filter(I, 11, 3)


	    # 1) Find local minima and maxima:
	    
	    # Split array: left and right parts
	    # Left
	    z_left = []
	    I_left = []
	    z_right = []
	    I_right = [] 
	    for k in range(len(z)):
	      if z[k]<0.:
		z_left.append(z[k])
		I_left.append(I_sav[k])
	      if z[k]>0.:
		z_right.append(z[k])
		I_right.append(I_sav[k])
	    
	    I_left_poly = np.polyfit(z_left, I_left, 11)
	    f_left = np.poly1d(I_left_poly)
	    
	    z_left_dust = opt.fminbound(-f_left, min(z_left),max(z_left))
	    #print max_z_left

	    
	    I_right_poly = np.polyfit(z_right, I_right, 11)
	    f_right = np.poly1d(I_right_poly)
	    
	    z_right_dust = opt.fminbound(-f_right, min(z_right),max(z_right))
	    #print max_z_right 
	    
	    z_mid = []
	    I_mid = []
	    
	    for k in range(len(z)):
	      if z[k]<z_right_dust and z[k]>z_left_dust:
		z_mid.append(z[k])
		I_mid.append(I_sav[k])   

	    I_mid_poly = np.polyfit(z_mid, I_mid, 11)
	    f_mid = np.poly1d(I_mid_poly)
	    
	    
	    z_cen_dust = opt.fminbound(f_mid, min(z_mid),max(z_mid))
	    #print min_z_mid 	    
	    
	    #plt.plot(z,I,'o',z_left,f_left(z_left),z_right,f_right(z_right),z_mid,f_mid(z_mid))    
	    #plt.show()         
  
	    return z_cen_dust,z_left_dust,z_right_dust


def bulge(B, x):
	return B[0]*np.exp( -(1.9992*B[2] - 0.3271)*((fabs(x)/B[1])**(1./B[2]) - 1.) )
'''
def disk_edge(B,z):
	#*** For edge-on disk SB in mag/arcsec^2 (along z-axis).  Sech^2-law. ***
	# B[0] = I0d
	# B[1] = z0 
	# B[2] = z_c
	z = np.array(z)	
	return B[0] * ( 1. / (np.cosh(fabs(z-B[2])/B[1]))**2   )
'''

def local_maxima(xval, yval):
    xval = np.asarray(xval)
    yval = np.asarray(yval)

    sort_idx = np.argsort(xval)
    yval = yval[sort_idx]
    gradient = np.diff(yval)
    maxima = np.diff((gradient > 0).view(np.int8))
    return np.concatenate((([0],) if gradient[0] < 0 else ()) +
                          (np.where(maxima == -1)[0] + 1,) +
                          (([len(yval)-1],) if gradient[-1] > 0 else ()))

'''
def crea_reg(x,y_up,y_down):
    coords = str()
    for k in range(0,len(x),1):
	coords = coords + str(x[k]) + ',' + str(y_up) + ','

    for k in range(len(x)-1,-1,-1):
      if k<=0:
	coords = coords + str(x[k]) + ',' + str(y_down)
      else:
	coords = coords + str(x[k]) + ',' + str(y_down) + ','

    fout = open("dust_lane.reg", "w")
    fout.write("polygon(%s)\n" % (coords))
    fout.close()  
'''

def get_slice(fitsName,maskName, xOrig, yOrig, posang, layer=0):
    """
    Function gets slice in fits file along specified line.
    Parameters:
        fitsName -- name of fits file.
        maskName -- name of the mask file.
        xOrig, yOrig -- coordinates of reference point
        posang [degrees] -- position angle of line. posang=0 -> vertical line (North-South slice).
                            positive posang values are for counterclockwise rotation, i.e. slice
                            with posang=45 is from upper-left corner till buttom right
        layer -- number of hdu in multilayered images
    """
    xOrig, yOrig = yOrig, xOrig
    hdu = pyfits.open(fitsName)
    data = hdu[layer].data
    
    hdu_mask = pyfits.open(maskName)
    data_mask = hdu_mask[layer].data    
    
    xSize, ySize = data.shape
    if xOrig==0:
      xOrig = xSize/2.
    if yOrig==0:
      yOrig = ySize/2.      
    rArray = []
    iArray = []
    posang += 90
    # tan of 90 and 270 degrees is infinity, so we have to avoid this nombers
    if (not (89.5 <= posang <= 90.5)) and (not (269.5 <= posang <= 270.5)):
        m = -tan(radians(posang))
        xbegin = xOrig - yOrig/m
        xend = min((ySize-yOrig)/m + xOrig, xSize)
        if xend < xbegin:
            xend, xbegin = xbegin, xend
        if xend > xSize:
            xend = xSize
        if xbegin < 0:
            xbegin = 0.0
        for x in linspace(xbegin, xend, xSize):
            y = m * (x-xOrig)+yOrig
            if (y<0) or (y>ySize-2) or (x>xSize-2) or (x<0):
                continue
            fx, ix = modf(x)
            fy, iy = modf(y)
            ix = int(ix)
            iy = int(iy)
            I = ((1.0-fx)*(1.0-fy)*data[ix, iy] + fx*(1.0-fy)*data[ix+1, iy] +
                fy*(1.0-fx)*data[ix, iy+1] + fx*fy*data[ix+1, iy+1])
            I_mask = ((1.0-fx)*(1.0-fy)*data_mask[ix, iy] + fx*(1.0-fy)*data_mask[ix+1, iy] +
                fy*(1.0-fx)*data_mask[ix, iy+1] + fx*fy*data_mask[ix+1, iy+1])    
	    if I_mask==0.:
	      rArray.append(copysign(hypot(x-xOrig, y-yOrig), xOrig-x))
	      iArray.append(I)
    else: # if posang is near 90 or 270 degrees then x is constant
        for y in arange(0, ySize-1):
            fx, ix = modf(xOrig)
            fy, iy = modf(y)
            ix = int(ix)
            iy = int(iy)
            I = ((1.0-fx)*(1.0-fy)*data[ix, iy] + fx*(1.0-fy)*data[ix+1, iy] +
                fy*(1.0-fx)*data[ix, iy+1] + fx*fy*data[ix+1, iy+1])
	    I_mask = ((1.0-fx)*(1.0-fy)*data_mask[ix, iy] + fx*(1.0-fy)*data_mask[ix+1, iy] +
                fy*(1.0-fx)*data_mask[ix, iy+1] + fx*fy*data_mask[ix+1, iy+1])
	    if I_mask==0.:    
	      rArray.append(y-yOrig)
	      iArray.append(I)
    rArray = np.array(rArray,float)
    iArray = np.array(iArray,float)
    return rArray, iArray

def rat(n1,n2):
  return max([n1,n2])/min([n1,n2])
  
def results(pars,zc):
  INCL = []
  DUST = []
  FIT_SOPH = []
  FIT_HYP = []
  FIT_EXP = []
  GAL_LUM = []
  GAL_LUM_soph = []
  GAL_LUM_hyp = []
  GAL_LUM_exp = []
  GAL_LUM = []
  for Pars in pars:
      for key in Pars.keys():
	if key=='INCL':
	  if not np.isnan(float(Pars['FIT_EXP'][0])):
	    INCL.append(Pars[key][0])
	if key=='DUST':
	    DUST.append([Pars[key][0],Pars[key][1],Pars[key][2],Pars[key][3]])
	if key=='FIT_SOPH':
	    if not np.isnan(float(Pars[key][0])):
	      FIT_SOPH.append([Pars[key][0],Pars[key][1],Pars[key][2],Pars[key][3],Pars[key][4],Pars[key][5],Pars[key][6]])
	      GAL_LUM_soph.append(Pars['GAL_LUM'][0])
	if key=='FIT_HYP':
	    if not np.isnan(float(Pars[key][0])):
	      FIT_HYP.append([Pars[key][0],Pars[key][1],Pars[key][2],Pars[key][3],Pars[key][4],Pars[key][5]])
	      GAL_LUM_hyp.append(Pars['GAL_LUM'][0])
	if key=='FIT_EXP':
	    if not np.isnan(float(Pars[key][0])):
	      FIT_EXP.append([Pars[key][0],Pars[key][1],Pars[key][2],Pars[key][3],Pars[key][4],Pars[key][5]])
	      GAL_LUM_exp.append(Pars['GAL_LUM'][0])
	if key=='GAL_LUM':
	    GAL_LUM.append(Pars[key][0])

  FIT_SOPH = np.array(FIT_SOPH)
  FIT_SOPH =FIT_SOPH.transpose()

  FIT_HYP = np.array(FIT_HYP)
  FIT_HYP =FIT_HYP.transpose()
  
  FIT_EXP = np.array(FIT_EXP)
  FIT_EXP =FIT_EXP.transpose()

  Incl = median(INCL); Incl_err = std(INCL)
  z0_soph = median(FIT_SOPH[2]); z0_soph_err = std(FIT_SOPH[2])
  n_soph = median(FIT_SOPH[1]); n_soph_err = std(FIT_SOPH[1])
  L_soph = sum(FIT_SOPH[4])
  A_soph = median(FIT_SOPH[5])
  chi2_soph = median(FIT_SOPH[6])
  L_gal_soph = sum(GAL_LUM_soph)
     
  z0_hyp = median(FIT_HYP[1]); z0_hyp_err = std(FIT_HYP[1])
  L_hyp = sum(FIT_HYP[3])
  A_hyp = median(FIT_HYP[4])
  chi2_hyp = median(FIT_HYP[5])
  L_gal_hyp = sum(GAL_LUM_hyp)  
 
  z0_exp = median(FIT_EXP[1]); z0_exp_err = std(FIT_EXP[1])
  L_exp = sum(FIT_EXP[3])
  A_exp = median(FIT_EXP[4])
  chi2_exp = median(FIT_EXP[5])
  L_gal_exp = sum(GAL_LUM_exp)
  
  
  L_gal = sum(GAL_LUM)

	    
  min_exp_z0  = z0_exp/2. 
  max_exp_z0 = z0_exp*2.
  min_hyp_z0 = z0_hyp/2.
  max_hyp_z0 = z0_hyp*2.
  min_soph_z0 = z0_soph/2.
  max_soph_z0 = z0_soph*2.
  
	    
  INCL = []
  DUST = []
  FIT_SOPH = []
  FIT_HYP = []
  FIT_EXP = []
  GAL_LUM = []
  GAL_LUM_soph = []
  GAL_LUM_hyp = []
  GAL_LUM_exp = []
  GAL_LUM = []
  for Pars in pars:
      for key in Pars.keys():
	if key=='INCL':
	  if not np.isnan(float(Pars['FIT_EXP'][0])) and float(Pars['FIT_EXP'][1])>min_exp_z0 and float(Pars['FIT_EXP'][1])<max_exp_z0 and not np.isnan(float(Pars[key][0])):
	    INCL.append(Pars[key][0])
	if key=='DUST':
	    DUST.append([Pars[key][0],Pars[key][1],Pars[key][2],Pars[key][3]])
	if key=='FIT_SOPH':
	    if not np.isnan(float(Pars[key][0])) and float(Pars[key][2])>min_soph_z0 and float(Pars[key][2])<max_soph_z0:
	      FIT_SOPH.append([Pars[key][0],Pars[key][1],Pars[key][2],Pars[key][3],Pars[key][4],Pars[key][5],Pars[key][6],Pars[key][7]])
	      GAL_LUM_soph.append(Pars['GAL_LUM'][0])
	if key=='FIT_HYP':
	    if not np.isnan(float(Pars[key][0])) and float(Pars[key][1])>min_hyp_z0 and float(Pars[key][1])<max_hyp_z0:
	      FIT_HYP.append([Pars[key][0],Pars[key][1],Pars[key][2],Pars[key][3],Pars[key][4],Pars[key][5],Pars[key][6]])
	      GAL_LUM_hyp.append(Pars['GAL_LUM'][0])
	if key=='FIT_EXP':
	    if not np.isnan(float(Pars[key][0])) and float(Pars[key][1])>min_exp_z0 and float(Pars[key][1])<max_exp_z0:
	      FIT_EXP.append([Pars[key][0],Pars[key][1],Pars[key][2],Pars[key][3],Pars[key][4],Pars[key][5],Pars[key][6]])
	      GAL_LUM_exp.append(Pars['GAL_LUM'][0])
	if key=='GAL_LUM':
	    GAL_LUM.append(Pars[key][0])
	       
	    
	    
	    
	    
	    
  FIT_SOPH = np.array(FIT_SOPH)
  FIT_SOPH =FIT_SOPH.transpose()

  FIT_HYP = np.array(FIT_HYP)
  FIT_HYP =FIT_HYP.transpose()
  
  FIT_EXP = np.array(FIT_EXP)
  FIT_EXP =FIT_EXP.transpose()

  DUST = np.array(DUST)
  DUST =DUST.transpose()
  
  y_coord_up = zc + median(DUST[0])+median(DUST[3])/2.
  y_coord_down = zc + median(DUST[0])-median(DUST[3])/2.


  Incl = median(INCL); Incl_err = std(INCL)
  z0_soph = median(FIT_SOPH[2]); z0_soph_err = std(FIT_SOPH[2])
  n_soph = median(FIT_SOPH[1]); n_soph_err = std(FIT_SOPH[1])
  L_soph = sum(FIT_SOPH[4])
  A_soph = median(FIT_SOPH[5])
  chi2_soph = median(FIT_SOPH[6])
  L_gal_soph = sum(GAL_LUM_soph)
     
  z0_hyp = median(FIT_HYP[1]); z0_hyp_err = std(FIT_HYP[1])
  L_hyp = sum(FIT_HYP[3])
  A_hyp = median(FIT_HYP[4])
  chi2_hyp = median(FIT_HYP[5])
  L_gal_hyp = sum(GAL_LUM_hyp)  
 
  z0_exp = median(FIT_EXP[1]); z0_exp_err = std(FIT_EXP[1])
  L_exp = sum(FIT_EXP[3])
  A_exp = median(FIT_EXP[4])
  chi2_exp = median(FIT_EXP[5])
  L_gal_exp = sum(GAL_LUM_exp)
  
  
  L_gal = sum(GAL_LUM)
  
  print bcolors.HEADER+"Rough estimation of dust lane and stellar disc parameters: "+ bcolors.ENDC
  print bcolors.OKGREEN+"inclination [deg]: "+ bcolors.ENDC+ "%.2f +/- %.2f" % (Incl,Incl_err)
  print bcolors.OKGREEN+"z0 [pix] (soph,hyp,exp): "+ bcolors.ENDC+ "%.2f +/- %.2f,\t%.2f +/- %.2f,\t%.2f +/- %.2f" % (z0_soph,z0_soph_err,z0_hyp,z0_hyp_err,z0_exp,z0_exp_err)  
  print bcolors.OKGREEN+"n_soph : "+ bcolors.ENDC+ "%.2f +/- %.2f" % (n_soph,n_soph_err)  
  print bcolors.OKGREEN+"A_soph (total) : "+ bcolors.ENDC+ "%.2f" % (L_soph/L_gal_soph)
  print bcolors.OKGREEN+"A_soph (mean) : "+ bcolors.ENDC+ "%.2f" % (A_soph)  
  print bcolors.OKGREEN+"A_hyp (total): "+ bcolors.ENDC+ "%.2f" % (L_hyp/L_gal_hyp)
  print bcolors.OKGREEN+"A_hyp (mean) : "+ bcolors.ENDC+ "%.2f" % (A_hyp)  
  print bcolors.OKGREEN+"A_exp (total): "+ bcolors.ENDC+ "%.2f" % (L_exp/L_gal_exp)
  print bcolors.OKGREEN+"A_exp (mean) : "+ bcolors.ENDC+ "%.2f" % (A_exp)  
  
  #n, bins, patches = plt.hist(FIT_HYP[1], 135, facecolor='green', alpha=0.5)
  #plt.show()
  
  chi2_list = [chi2_soph,chi2_hyp,chi2_exp]
  min_chi2 = min(chi2_list)
  if chi2_list.index(min_chi2)==0:
    print "BEST Chi2 is for Sophisticated formula"
    plt.figure(0,figsize=(5, 5))
    plt.xlabel(r' r (pixels) ', fontsize=15)
    plt.ylabel(r' z$_0$ (pixels) ', fontsize=15)
    plt.plot(FIT_SOPH[7],FIT_SOPH[2],'o')
    plt.savefig('z0_r.eps', transparent = False, dpi=300)
    best_chi2 = "Soph"
    A_ext_mean = A_soph
    A_ext_total = L_soph/L_gal_soph
  elif chi2_list.index(min_chi2)==1:
    print "BEST Chi2 is for Sech^2 formula"
    plt.figure(0,figsize=(5, 5))
    plt.xlabel(r' r (pixels) ', fontsize=15)
    plt.ylabel(r' z$_0$ (pixels) ', fontsize=15)
    plt.plot(FIT_HYP[6],FIT_HYP[1],'o')
    plt.savefig('z0_r.eps', transparent = False, dpi=300)
    best_chi2 = "Sech^2"
    A_ext_mean = A_hyp
    A_ext_total = L_hyp/L_gal_hyp    
  else:
    plt.figure(0,figsize=(5, 5))
    plt.xlabel(r' r (pixels) ', fontsize=15)
    plt.ylabel(r' z$_0$ (pixels) ', fontsize=15)
    plt.plot(FIT_EXP[6],FIT_EXP[1],'o')
    plt.savefig('z0_r.eps', transparent = False, dpi=300)      
    print "BEST Chi2 is for Exp formula"
    best_chi2 = "Exp"
    A_ext_mean = A_exp
    A_ext_total = L_exp/L_gal_exp 
 

  plt.figure(1,figsize=(5, 5))
  plt.xlabel(r' r (pixels) ', fontsize=15)
  plt.ylabel(r' n ', fontsize=15)
  plt.ylim(0,110)
  plt.plot(FIT_SOPH[7],FIT_SOPH[1],'o')
  plt.savefig('z0_n.eps', transparent = False, dpi=300)
  plt.clf()
  plt.close()  
  
  PARS_RES = {}
  PARS_RES["INCL"]=[Incl,Incl_err]
  PARS_RES["z0"]=[z0_soph,z0_soph_err,z0_hyp,z0_hyp_err,z0_exp,z0_exp_err]
  PARS_RES["n"]=[n_soph,n_soph_err]
  PARS_RES["BEST_CHI2"]=[best_chi2]
  PARS_RES["Extinction"]=[A_ext_mean,A_ext_total]
  
  
  return y_coord_up,y_coord_down,PARS_RES

def find_nearest(value,array):
	idx = (np.abs(np.array(array)-value)).argmin()
	return array[idx]

def profile(input_image,mask_image,pix2sec,m0,R_min,R_max,sign,x0,y0):
      layer = 1
      hdu = pyfits.open(input_image)
      data = hdu[0].data     
      ySize, xSize = data.shape
      
      hdu_mask = pyfits.open(mask_image)
      mask_data = hdu_mask[0].data  

      if layer==1:
	color1 = 'grey'
	color2 = 'white'
	mecolor = 'black'
	msize=5
      else:
	color1 = 'red'
	color2 = 'red'
	mecolor = 'red'
	msize=3

      # Photometrical cut along the major axis
      r_ini, i_ini = get_slice(input_image,mask_image, x0, y0, 0., layer=0)
      '''
      i_ini = []
      r_ini = []
      for k in range(xSize):
	      if mask_data[int(y0),k]==0.:
		i_ini.append(data[int(y0),k])
		r_ini.append(k-int(x0))
      '''
      r_ini = np.array(r_ini)
      i_ini = np.array(i_ini)      


      r = []
      I = []
      for k in range(len(r_ini)):
	if sign=="+" and r_ini[k]>=0.:
	  r.append(r_ini[k])
	  I.append(i_ini[k])
	if sign=="-" and r_ini[k]<0.:
	  r.append(-r_ini[k])
	  I.append(i_ini[k]) 
 
      r = np.array(r)
      I = np.array(I)
      #r = r*pix2sec
      mag = m0 - 2.5*log10(I)#+ 5.*log10(pix2sec)


      R = []
      Mag = []
      for k in range(len(r)):
	if r[k]>=R_min and mag[k]!=float(inf) and math.isnan(mag[k])==False and r[k]<=R_max-R_max/30.:		#### NOTE: -R_max/30. - go away from the edge of the profile where the noise is very high
	  R.append(r[k])
	  Mag.append(mag[k])
      if sign=='+':
	plot(R, Mag,'o',color=color2,markeredgecolor=mecolor,markersize=msize)
      else:
	plot(-np.array(R), Mag,'o',color=color2,markeredgecolor=mecolor,markersize=msize)
      return R,Mag


def profile_z(input_image,mask_image,pix2sec,m0,R_min,R_max,zmin,Minor_axis,sign,x0,y0):
      layer = 1
      hdu = pyfits.open(input_image)
      data = hdu[0].data     
      ySize, xSize = data.shape
      
      hdu_mask = pyfits.open(mask_image)
      mask_data = hdu_mask[0].data
      
      for k in range(ySize):
	for i in range(xSize):
	  if mask_data[k,i]!=0.:
	    data[k,i]=0.

      if layer==1:
	color1 = 'grey'
	color2 = 'white'
	mecolor = 'black'
	msize=5
      else:
	color1 = 'cyan'
	color2 = 'cyan'
	mecolor = 'cyan'
	msize=4

      # Photometrical cut along the minor axis

      i_ini = []
      r_ini = []
      for y in range(0,ySize,1):
	i_ini.append(sum(data[y,x0-R_max:x0-R_min])+sum(data[y,x0+R_min:x0+R_max]))
	r_ini.append(y-y0)
      r_ini = np.array(r_ini)
      i_ini = np.array(i_ini) 
      #plt.plot(r_ini,i_ini)
      #plt.show()
      r = []
      I = []
      for k in range(len(r_ini)):
	if sign=="+" and r_ini[k]>=0. and r_ini[k]>=zmin and r_ini[k]<Minor_axis:
	  r.append(r_ini[k])
	  I.append(i_ini[k])
	if sign=="-" and r_ini[k]<0. and r_ini[k]<=-zmin and r_ini[k]>-Minor_axis:
	  r.append(-r_ini[k])
	  I.append(i_ini[k])
	if sign=="0" and fabs(r_ini[k])>=zmin:
	  r.append(r_ini[k])
	  I.append(i_ini[k])
	  
      r = np.array(r)
      I = np.array(I)
      #r = r*pix2sec
      mag = m0 - 2.5*log10(I)#+ 5.*log10(pix2sec)


      R = []
      Mag = []
      for k in range(len(r)):
	if mag[k]!=float(inf) and math.isnan(mag[k])==False:
	  R.append(float(r[k]))
	  Mag.append(mag[k])
      if sign=='+':
	plot(R, Mag,'o',color=color2,markeredgecolor=mecolor,markersize=msize)
      if sign=='-':
	plot(-np.array(R), Mag,'o',color=color2,markeredgecolor=mecolor,markersize=msize)

      return R,Mag


def summed_hor_profile(input_image,xc,yc,Radius,zmax,m0,sign='both',min_radius=0.):
	    hdu = pyfits.open(input_image)
	    data = hdu[0].data
	    Radius = int(Radius)
	    min_radius = int(min_radius)
	    zmax = int(zmax)

  
	    #hdu_mask = pyfits.open(mask_image)
	    #data_mask = hdu_mask[layer].data    
    
	    
	    ySize, xSize = data.shape
	    I = []
	    r = []
	    if min_radius==0.:
	      X = range(xc-Radius,xc+Radius,1)
	    else:
	      X = range(xc-Radius,xc-min_radius,1) + range(xc+min_radius,xc+Radius,1)

	    for x in X:
	      if sign=='both':
		SUM = sum(data[yc-zmax:yc+zmax,x])
		I.append(SUM)
		r.append(x-xc)
	      if sign=='+':
		if x-xc>=0:
		  SUM = sum(data[yc-zmax:yc+zmax,x])
		  I.append(SUM)
		  r.append(x-xc)
	      if sign=='-':
		if x-xc<0:
		  SUM = sum(data[yc-zmax:yc+zmax,x])
		  I.append(SUM)
		  r.append(fabs(x-xc))  

	    r = np.array(r)
	    I = np.array(I)

	    # remove nans
	    r_g = []; I_g = []
	    for k in range(len(r)):
	      if np.isnan(m0 - 2.5*log10(I[k]))==False and np.isinf(m0 - 2.5*log10(I[k]))==False:
		r_g.append(r[k])
		I_g.append(m0 - 2.5*log10(I[k]))
	    return np.array(r_g,float),np.array(I_g,float)

def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])
  
def piecewise_kn(x, x0, y0, h1, h2):
    return np.piecewise(x, [x < x0], [lambda x:2.5*log10((fabs(x)/h1))+2.5*log10(sp.special.kn(1,fabs(x)/h1)) - 2.5*log10((fabs(x0)/h1))-2.5*log10(sp.special.kn(1,fabs(x0)/h1)) + y0, lambda x:2.5*log10((fabs(x)/h2))+2.5*log10(sp.special.kn(1,fabs(x)/h2)) - 2.5*log10((fabs(x0)/h2))-2.5*log10(sp.special.kn(1,fabs(x0)/h2)) + y0])  


def disk_edge_r(B,x):
	#*** For edge-on disk SB in mag/arcsec^2 (along r-axis). ***
	I0d = B[0]
	h = B[1]
	x0 = B[2]
	x = np.array(x)
	return I0d * fabs(x-x0)/h * sp.special.kn(1,fabs(x-x0)/h)


def crea_reg(xc,yc,sma_b,smb_b,sma_d,smb_d,R_break):
    fout = open("bulge_disc_boxes.reg", "w")
    fout.write("box(%.1f,%.1f,%.1f,%.1f,%.1f) # color=blue width=2\n" % (xc,yc,2.*sma_b,2.*smb_b,0))
    fout.write("box(%.1f,%.1f,%.1f,%.1f,%.1f) # color=black width=2\n" % (xc,yc,2.*sma_d,2.*smb_d,0))
    fout.write("box(%.1f,%.1f,%.1f,%.1f,%.1f) # color=black width=2\n" % (xc,yc,2.*R_break,2.*yc,0)) #### NOTE: NEW LINE!!!
    fout.close()  

def ds9_show(fileToShow):
    ds9Proc = subprocess.Popen(["ds9", fileToShow,
                                "-scale", "histequ","-cmap","Rainbow"])
    #ds9Proc.wait()

def main(input_image, mask_image, mag_level, Radius, Minor_axis, zmin, m0, pix2sec, xc=0., zc=0., inter_mode=True, min_disc_radius=0., Profile='summed'):  
      print bcolors.OKBLUE+'\n\n************ Bulge and disc domination regions (c) Mosenkov A. 2015 ************' + bcolors.ENDC
      print "Analyzing..."

      hdu = pyfits.open(input_image)
      data = hdu[0].data
      ySize, xSize = data.shape

      if zc==0.:
	zc=ySize/2.
      if xc==0.:
	xc=xSize/2.

      
      hdu_mask = pyfits.open(mask_image)
      mask_data = hdu_mask[0].data
      
      
      PARS = {}
      PARS["BULGE_DOM"] = [Minor_axis/0.7,Minor_axis]
      PARS["DISC_DOM"] = [float(nan),float(nan),Minor_axis/0.5]
      PARS["INNER_DISC"] = [float(nan),float(nan),float(nan),float(nan),float(nan),float(nan)]
      PARS["OUTER_DISC"] = [float(nan),float(nan),float(nan),float(nan),float(nan),float(nan)]
      PARS["BREAK_RADIUS"] = [float(nan),float(nan),float(nan)]
      PARS["BREAK_Z"] = [float(nan),float(nan),float(nan)]
      PARS["THIN_DISC"] = [float(nan),float(nan),float(nan)]
      PARS["THICK_DISC"] = [float(nan),float(nan),float(nan)]
      PARS["DELTA_M0D"] = [float(nan)]
      PARS["THICKNESS"] = [float(nan)]

      I_min = 10**(0.4*(m0-mag_level+5.*log10(pix2sec)))

      repeat=True 
      while repeat!=False:
	    plt.figure(3,figsize=(5, 5))
	    plt.xlabel(r' r (pix) ', fontsize=15)
	    plt.ylabel(r' $\mu$ (mag pix$^{-2}$) ', fontsize=15)
	    plt.gca().invert_yaxis()   
	    
	    
	    
	    #**************************************************************************************
	    #### Bulge region in vertical direction:
	    # Positive z
	    I = []
	    z = []
	    for k in range(ySize):
		      # Select only pixels with SB larger than the limit (to remove noise or faint second disc or halo):
		      if data[k,int(xc)]>=I_min and k-zc>zmin and mask_data[k,int(xc)]==0. and k-zc<Minor_axis:
			I.append(data[k,xc])
			z.append(k-zc)
			
	    I = m0 - 2.5*log10(np.array(I, float))
	    z = np.array(z, float)
	    Res_chi2 = []
	    Res_p = []
	    for zz in z:
	      ind = list(z).index(zz)
	      p , e = optimize.curve_fit(piecewise_linear, z, I,[zz,I[ind],1,1])
	      chi2 = chi2_func(I,piecewise_linear(z, *p))
	      Res_chi2.append(chi2)
	      Res_p.append(p)

	    results_pos = Res_p[Res_chi2.index(min(Res_chi2))]
	    results_pos_bulge = results_pos

	    plt.plot(z,I,'o',color='grey') #,z,piecewise_linear(z, *results_pos),color='blue') 

	    # Negative z
	    I = []
	    z = []
	    for k in range(ySize):
		      # Select only pixels with SB larger than the limit (to remove noise or faint second disc or halo):
		      if data[k,int(xc)]>=I_min and k-zc<-zmin and mask_data[k,int(xc)]==0. and k-zc>-Minor_axis:
			I.append(data[k,xc])
			z.append(k-zc)
			
	    I = m0 - 2.5*log10(np.array(I, float))
	    z = np.array(z, float)

	    Res_chi2 = []
	    Res_p = []
	    for zz in z:
	      ind = list(z).index(zz)
	      p , e = optimize.curve_fit(piecewise_linear, z, I,[zz,I[ind],1,1])
	      chi2 = chi2_func(I,piecewise_linear(z, *p))
	      Res_chi2.append(chi2)
	      Res_p.append(p)

	    results_neg = Res_p[Res_chi2.index(min(Res_chi2))]
	    results_neg_bulge = results_neg

	    plt.plot(z,I,'o',color='grey') #,z,piecewise_linear(z, *results_neg),color='blue') 
	    
	    R_bulge_mi = (fabs(results_neg[0])+fabs(results_pos[0]))/2.
	    #/**************************************************************************************













	    if Profile!='summed':
			#**************************************************************************************
			#### PLOT ALL THE POINTS ON THE PROFILE EXCEPTING MASKED ONES:
			I = []
			r = []
			for k in range(xSize):
				# Select only pixels with SB larger than the limit (to remove noise or faint second disc or halo):
				if mask_data[int(zc),k]==0.:
				  I.append(data[int(zc),k])
				  r.append(k-xc)
			I = m0 - 2.5*log10(np.array(I, float))
			r = np.array(r, float)      
			plt.plot(r,I,'o',color='white') 
			min_mag =  min(I[np.logical_not(np.isnan(I))]) - 0.5
			#/**************************************************************************************
		    
		    
		    
			#**************************************************************************************
			#### Estimation of the disc region:    
			# Positive r
			I = []
			r = []
			for k in range(xSize):
				# Select only pixels with SB larger than the limit (to remove noise or faint second disc or halo):
				if data[int(zc),k]>=I_min and k-xc>0. and k-xc<Radius/2 and mask_data[int(zc),k]==0.:
				  I.append(data[int(zc),k])
				  r.append(k-xc)
			II = I
			I = m0 - 2.5*log10(np.array(I, float))
			r = np.array(r, float)

			Res_chi2 = []
			Res_p = []

			for rr in r:
			  ind = list(r).index(rr)
			  p , e = optimize.curve_fit(piecewise_linear, r, I,[rr,I[ind],1,1])
			  chi2 = chi2_func(I,piecewise_linear(r, *p))
			  Res_chi2.append(chi2)
			  Res_p.append(p)

			results_pos = Res_p[Res_chi2.index(min(Res_chi2))]


			# Negative r
			I = []
			r = []
			for k in range(xSize):
				if data[int(zc),k]>=I_min and k-xc<0. and k-xc>-Radius/2 and mask_data[int(zc),k]==0.:
				  I.append(data[int(zc),k])
				  r.append(fabs(k-xc))
			II = I
			I = m0 - 2.5*log10(np.array(I, float))
			r = np.array(r, float)

			Res_chi2 = []
			Res_p = []
			for rr in r:
			  ind = list(r).index(rr)
			  p , e = optimize.curve_fit(piecewise_linear, r, I,[rr,I[ind],1,1])
			  chi2 = chi2_func(I,piecewise_linear(r, *p))
			  Res_chi2.append(chi2)
			  Res_p.append(p)

			results_neg = Res_p[Res_chi2.index(min(Res_chi2))]
			

			R_bulge_ma = (fabs(results_neg[0])+fabs(results_pos[0]))/2.
			
			m_bulge_ma = (fabs(results_neg[1])+fabs(results_pos[1]))/2.
			
			z_bul = arange(0.,Minor_axis,0.01)
			aver_bul_prof = (piecewise_linear(-z_bul, *results_neg_bulge) + piecewise_linear(z_bul, *results_pos_bulge)) / 2.
			#plt.plot(z_bul,aver_bul_prof,'*')
			#plt.show()
			#print m_bulge_ma
			R_bulge_mi = z_bul[np.where(aver_bul_prof==find_nearest(m_bulge_ma,aver_bul_prof))][0]

			
			#print R_bulge_ma
			PARS["BULGE_DOM"] = [max([R_bulge_ma,R_bulge_mi]),min([R_bulge_ma,R_bulge_mi])]
			
			R_disc_min = 1.2*max([R_bulge_ma,R_bulge_mi])		#### NOTE: Coefficient here should be varified!

			if inter_mode==True and repeat==2:      
			  R_disc_min = max([R_disc_min_left_man,R_disc_min_right_man])
			if min_disc_radius!=0. and repeat!=2:
			  R_disc_min = min_disc_radius
			print '\n\nR_disc_min=',R_disc_min
			h_d = mean([fabs(1.086/results_pos[3] / 1.3),fabs(1.086/results_neg[3] / 1.3)])
			I0_d = mean([10**(0.4*(m0-results_pos[1]+results_pos[3]*results_pos[1])) / 1.9,10**(0.4*(m0-results_neg[1]+results_neg[3]*results_neg[1])) / 1.9])
			m0_d = m0 - 2.5*log10(I0_d)# + 5.*log10(pix2sec)

			# More precise fit to find disc parameters:
			I = []
			r = []
			for k in range(xSize):
				if data[int(zc),k]>=I_min and fabs(k-xc)<=Radius and fabs(k-xc)>=R_disc_min and mask_data[int(zc),k]==0.:
				  I.append(data[int(zc),k])
				  r.append(k-xc)
			h_fin,I0d_fin = fit_eon_radial_prof(r, I, 10**(0.4*(m0-m0_d)), h_d)
			PARS["DISC_DOM"] = [h_fin,m0-2.5*log10(I0d_fin)+5.*log10(pix2sec),R_disc_min]
			
			
			print bcolors.HEADER+"Rough estimation of bulge and disc domination: "+ bcolors.ENDC
			print bcolors.OKGREEN+"Bulge [pix] (Rmax,Rmin): "+ bcolors.ENDC+ "%.2f,\t %.2f" % (max([R_bulge_ma,R_bulge_mi]),min([R_bulge_ma,R_bulge_mi]))
			print bcolors.OKGREEN+"Disc (h [pix], m0d [mag arcsec^-2], Rmin [pix]) : "+ bcolors.ENDC+ "%.2f,\t%.2f,\t%.2f" % (h_fin,m0-2.5*log10(I0d_fin)+5.*log10(pix2sec),R_disc_min)
			axvline(x=R_disc_min,color='red',linewidth=2,ls='--')
			axvline(x=-R_disc_min,color='red',linewidth=2,ls='--')
			
			'''
			if inter_mode==True:
			  ds9_show(input_image)
			  plt.show()
			  
			  R_disc_min_left_man = float(input('Please enter left radius of disc:'))
			  R_disc_min_right_man = float(input('Please enter right radius of disc:'))
			  print R_disc_min_left_man,R_disc_min_right_man
			'''
			#/**************************************************************************************
			
			
	      
			#**************************************************************************************
			#### Disc breaks
			Rmin = R_disc_min

			# Positive radii:
			R,Mag = profile(input_image,mask_image,pix2sec,m0,Rmin,Radius,"+",xc,zc)
			MAG_RIGHT = Mag
			R_RIGHT = R

			if inter_mode==True and repeat==2:
			  x_break_right = R_break_right_man
			  y_break_right = Mag_break_right_man
			  def piecewise_linear_fix1(x, k1, k2):
			      return np.piecewise(x, [x < x_break_right], [lambda x:k1*x + y_break_right-k1* x_break_right, lambda x:k2*x +y_break_right-k2* x_break_right])
			  results_pos1, e = optimize.curve_fit(piecewise_linear_fix1, R, Mag,[1,1]) 
			  results_pos = [R_break_right_man,Mag_break_right_man,results_pos1[0],results_pos1[1]]
			else:
			  Res_chi2 = []
			  Res_p = []
			  for rr in R:
			    ind = list(R).index(rr)
			    p , e = optimize.curve_fit(piecewise_linear, R, Mag,[rr,Mag[ind],1,1])
			    chi2 = chi2_func(Mag,piecewise_linear(R, *p))
			    Res_chi2.append(chi2)
			    Res_p.append(p)
			    
			  
			  results_pos = Res_p[Res_chi2.index(min(Res_chi2))]
			  
			Rlim_pos =  int(max(R))  
			rd_pos = np.linspace(Rmin, int(max(R)), 100)
		  
			
			# Negative radii:   
			R,Mag = profile(input_image,mask_image,pix2sec,m0,Rmin,Radius,"-",xc,zc)
			MAG_LEFT = Mag
			R_LEFT = R

			if inter_mode==True and repeat==2:
			  x_break_left = R_break_left_man
			  y_break_left = Mag_break_left_man
			  def piecewise_linear_fix2(x, k1, k2):
			      return np.piecewise(x, [x < x_break_left], [lambda x:k1*x + y_break_left-k1* x_break_left, lambda x:k2*x +y_break_left-k2* x_break_left])
			  results_neg1, e = optimize.curve_fit(piecewise_linear_fix2, R, Mag,[1,1])
			  results_neg = [R_break_left_man,Mag_break_left_man,results_neg1[0],results_neg1[1]]
			else:
			  Res_chi2 = []
			  Res_p = []
			  for rr in R:
			    ind = list(R).index(rr)
			    try:
			      p , e = optimize.curve_fit(piecewise_linear, R, Mag,[rr,Mag[ind],1,1])
			      chi2 = chi2_func(Mag,piecewise_linear(R, *p))
			      Res_chi2.append(chi2)
			      Res_p.append(p)
			    except:
			      zz = 1
			  
			  results_neg = Res_p[Res_chi2.index(min(Res_chi2))]
			  
			Rlim_neg =  int(max(fabs(R)))
			rd_neg = np.linspace(Rmin, int(max(R)), 100)
			plt.ylim(max(Mag),min_mag)
			plt.xlim(-max([Rlim_pos,int(max(R))]),max([Rlim_pos,int(max(R))]))



			rd_gen = np.linspace(-Radius, Radius, 100)
			plt.plot(rd_gen, m0-2.5*log10(disk_edge_r([I0d_fin,h_fin,0.],rd_gen)),color='black')
		
			
			#********************************************************************
			print Rlim_pos,results_pos[0]
			print Rlim_neg,results_neg[0]
			# RESULTS:
			check_pos = 1
			check_neg = 1
			check_all = 1
			if not check_sign([results_pos[2],results_pos[3]]) or rat(Radius,results_pos[0])<1.1:
			    print bcolors.FAIL+'Positive radii: Inner and outer disc profiles have different signs!'+ bcolors.ENDC
			    plt.plot(rd_pos, piecewise_linear(rd_pos, *results_pos),color='blue')
			    axvline(x=results_pos[0],color='blue',linewidth=2,ls='--')
			    check_pos = 0
			else:
			    plt.plot(rd_pos, piecewise_linear(rd_pos, *results_pos),color='red')
			    axvline(x=results_pos[0],color='black',linewidth=2,ls='--')

			if not check_sign([results_neg[2],results_neg[3]]) or rat(Radius,results_neg[0])<1.1:
			    print bcolors.FAIL+'Negative radii: Inner and outer disc profiles have different signs!'+ bcolors.ENDC
			    plt.plot(-rd_neg, piecewise_linear(rd_neg, *results_neg),color='blue')
			    axvline(x=-results_neg[0],color='blue',linewidth=2,ls='--')
			    check_neg = 0
			else:
			    plt.plot(-rd_neg, piecewise_linear(rd_neg, *results_neg),color='red')
			    axvline(x=-results_neg[0],color='black',linewidth=2,ls='--')

			if not check_sign([results_pos[2],results_pos[3],results_neg[2],results_neg[3]]) and check_pos==1 and check_neg==1:
			    print bcolors.FAIL+'All radii: Inner and outer disc profiles have different signs!'+ bcolors.ENDC
			    check_all = 0

			'''
			plt.plot(rd_pos, piecewise_linear(rd_pos, *results_pos),color='red')
			axvline(x=results_pos[0],color='black',linewidth=2,ls='--')    
			plt.plot(-rd_neg, piecewise_linear(rd_neg, *results_neg),color='red')
			axvline(x=-results_neg[0],color='black',linewidth=2,ls='--')    
			'''
			if inter_mode==True:
			  #ds9_show(input_image)
			  fiG = plt.gcf()
			  plt.show()
			  s = raw_input('Are you happy with this [YES/no]?') or "yes"
			  if s!='yes':
			    input1 = raw_input('Please enter left and right radius of the disc:') or str(R_disc_min)+','+str(R_disc_min)
			    input2 =  raw_input('Please enter left break radius of the disc (x,y):') or str(results_neg[0])+','+str(MAG_LEFT[R_LEFT.index(find_nearest(results_neg[0],R_LEFT))])
			    input3 =  raw_input('Please enter right break radius of the disc (x,y):') or str(results_pos[0])+','+str(MAG_RIGHT[R_RIGHT.index(find_nearest(results_pos[0],R_RIGHT))])
			    
			    R_disc_min_left_man,R_disc_min_right_man = input1.split(',')
			    R_break_left_man,Mag_break_left_man = input2.split(',')
			    R_break_right_man,Mag_break_right_man = input3.split(',')

			    print '\n You entered:'
			    print R_disc_min_left_man,R_disc_min_right_man
			    print R_break_left_man,Mag_break_left_man
			    print R_break_right_man,Mag_break_right_man
			  
			    R_disc_min_left_man = float(R_disc_min_left_man)
			    R_disc_min_right_man = float(R_disc_min_right_man)
			    R_break_left_man = float(R_break_left_man)
			    R_break_right_man = float(R_break_right_man)
			    Mag_break_left_man = float(Mag_break_left_man)
			    Mag_break_right_man = float(Mag_break_right_man)
			    repeat=2
			  else:
			    if check_pos == 0:
			      results_pos[0]=float('nan')
			      results_pos[1]=float('nan')
			      results_pos[2]=float('nan')
			      results_pos[3]=float('nan')
			      
			    if check_neg == 0: 
			      results_neg[0]=float('nan')
			      results_neg[1]=float('nan')
			      results_neg[2]=float('nan')
			      results_neg[3]=float('nan')
			    if check_all == 0:
			      results_pos[0]=float('nan')
			      results_pos[1]=float('nan')
			      results_pos[2]=float('nan')
			      results_pos[3]=float('nan')
			      results_neg[0]=float('nan')
			      results_neg[1]=float('nan')
			      results_neg[2]=float('nan')
			      results_neg[3]=float('nan')
			      
			    R_disc_min_left_man = R_disc_min
			    R_disc_min_right_man = R_disc_min
			    R_break_left_man = results_neg[0]
			    R_break_right_man = results_pos[0]
			    PARS["PARS_MAN"] = [R_disc_min_left_man,R_disc_min_right_man,R_break_left_man,R_break_right_man]
			    repeat = False
			    #plt.clf()
			    #plt.close() 
			else:
			  fiG = plt.gcf()
			  repeat=False

			fiG.savefig('bulge_disc_dom.eps', transparent = False, dpi=300)
			plt.clf()
			plt.close()  
		    
			R_break_pos = results_pos[0]
			R_break_neg = results_neg[0]
			R_break_mean = (R_break_pos+R_break_neg) / 2.
			
			# Inner disc:
			h_in_pos = 1.086/results_pos[2] / 1.3
			h_in_neg = 1.086/results_neg[2] / 1.3
			h_in_mean = (h_in_pos+h_in_neg) / 2.
			m0d_in_pos = results_pos[1]-results_pos[2]*results_pos[1]+5.*log10(pix2sec) + 2.5*log10(1.9)
			m0d_in_neg = results_neg[1]-results_neg[2]*results_neg[1]+5.*log10(pix2sec) + 2.5*log10(1.9)
			m0d_in_mean = m0-2.5*log10((10**(0.4*(m0-m0d_in_pos)) + 10**(0.4*(m0-m0d_in_neg)))/2.)
			
			# Outer disc:
			h_out_pos = 1.086/results_pos[3] / 1.3
			h_out_neg = 1.086/results_neg[3] / 1.3
			h_out_mean = (h_out_pos+h_out_neg) / 2.
			m0d_out_pos = results_pos[1]-results_pos[3]*results_pos[1]+5.*log10(pix2sec) + 2.5*log10(1.9)
			m0d_out_neg = results_neg[1]-results_neg[3]*results_neg[1]+5.*log10(pix2sec) + 2.5*log10(1.9)
			m0d_out_mean = m0-2.5*log10((10**(0.4*(m0-m0d_out_pos)) + 10**(0.4*(m0-m0d_out_neg)))/2.)      

			print bcolors.HEADER+"\n\nINNER AND OUTER DISCS: Results of the fitting (Rmin=%.1f, Rmax=%.1f):" % (Rmin,Radius) + bcolors.ENDC
			print bcolors.OKGREEN+"Break radius [pix]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)" % (R_break_mean,R_break_neg,R_break_pos)
			print bcolors.OKGREEN+"h_in [pix]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)" % (h_in_mean,h_in_neg,h_in_pos)
			print bcolors.OKGREEN+"m0d_in [mag arcsec^-2]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)"  % (m0d_in_mean,m0d_in_neg,m0d_in_pos)
			print bcolors.OKGREEN+"h_out [pix]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)" % (h_out_mean,h_out_neg,h_out_pos)
			print bcolors.OKGREEN+"m0d_out [mag arcsec^-2]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)\n"  % (m0d_out_mean,m0d_out_neg,m0d_out_pos)       
			#/**************************************************************************************

	    else:
			#**************************************************************************************
			#### PLOT ALL THE POINTS ON THE PROFILE EXCEPTING MASKED ONES:
			r,I = summed_hor_profile(input_image,xc,zc,Radius,Minor_axis,m0,sign='both')    
			plt.plot(r,I,'o',color='white') 
			min_mag =  min(I[np.logical_not(np.isnan(I))]) - 0.5
			#/**************************************************************************************
		    
		    
		    
			#**************************************************************************************
			#### Estimation of the disc region:    
			# Positive r
			r,I = summed_hor_profile(input_image,xc,zc,Radius,Minor_axis,m0,sign='+')    

			Res_chi2 = []
			Res_p = []

			for rr in r:
			  ind = list(r).index(rr)
			  p , e = optimize.curve_fit(piecewise_linear, r, I,[rr,I[ind],1,1])
			  chi2 = chi2_func(I,piecewise_linear(r, *p))
			  Res_chi2.append(chi2)
			  Res_p.append(p)

			results_pos = Res_p[Res_chi2.index(min(Res_chi2))]
			#plt.plot(r,piecewise_linear(r, *results_pos),'*',color='red')
			#plt.show()
			#exit()


			# Negative r
			r,I = summed_hor_profile(input_image,xc,zc,Radius,Minor_axis,m0,sign='-')    

			Res_chi2 = []
			Res_p = []
			for rr in r:
			  ind = list(r).index(rr)
			  p , e = optimize.curve_fit(piecewise_linear, r, I,[rr,I[ind],1,1])
			  chi2 = chi2_func(I,piecewise_linear(r, *p))
			  Res_chi2.append(chi2)
			  Res_p.append(p)

			results_neg = Res_p[Res_chi2.index(min(Res_chi2))]
			#plt.plot(r,I,'*',color='red')
			#plt.show()
			#exit()

			

			R_bulge_ma = (fabs(results_neg[0])+fabs(results_pos[0]))/2.
			
			m_bulge_ma = (fabs(results_neg[1])+fabs(results_pos[1]))/2.
			
			z_bul = arange(0.,Minor_axis,0.01)
			aver_bul_prof = (piecewise_linear(-z_bul, *results_neg_bulge) + piecewise_linear(z_bul, *results_pos_bulge)) / 2.
			#plt.plot(z_bul,aver_bul_prof,'*')
			#plt.show()
			#print m_bulge_ma
			R_bulge_mi = z_bul[np.where(aver_bul_prof==find_nearest(m_bulge_ma,aver_bul_prof))][0]

			
			#print R_bulge_ma
			PARS["BULGE_DOM"] = [max([R_bulge_ma,R_bulge_mi]),min([R_bulge_ma,R_bulge_mi])]
			
			R_disc_min = 1.2*max([R_bulge_ma,R_bulge_mi])		#### NOTE: Coefficient here should be varified!

			if inter_mode==True and repeat==2:      
			  R_disc_min = max([R_disc_min_left_man,R_disc_min_right_man])
			if min_disc_radius!=0. and repeat!=2:
			  R_disc_min = min_disc_radius
			print '\n\nR_disc_min=',R_disc_min
			h_d = mean([fabs(1.086/results_pos[3] / 1.3),fabs(1.086/results_neg[3] / 1.3)])
			I0_d = mean([10**(0.4*(m0-results_pos[1]+results_pos[3]*results_pos[1])) / 1.9,10**(0.4*(m0-results_neg[1]+results_neg[3]*results_neg[1])) / 1.9])
			m0_d = m0 - 2.5*log10(I0_d)

			# More precise fit to find disc parameters:
			r,I = summed_hor_profile(input_image,xc,zc,Radius,Minor_axis,m0,sign='both',min_radius=R_disc_min)
			#plt.plot(r,I,'*',color='red')
			#plt.show()
			#exit()
			h_fin,I0d_fin = fit_eon_radial_prof(r, 10**(0.4*(m0-I)), 10**(0.4*(m0-m0_d)), h_d)
			PARS["DISC_DOM"] = [h_fin,m0-2.5*log10(I0d_fin)+5.*log10(pix2sec),R_disc_min]
			#plt.show()
			#exit()
			
			print bcolors.HEADER+"Rough estimation of bulge and disc domination: "+ bcolors.ENDC
			print bcolors.OKGREEN+"Bulge [pix] (Rmax,Rmin): "+ bcolors.ENDC+ "%.2f,\t %.2f" % (max([R_bulge_ma,R_bulge_mi]),min([R_bulge_ma,R_bulge_mi]))
			print bcolors.OKGREEN+"Disc (h [pix], m0d [mag arcsec^-2], Rmin [pix]) : "+ bcolors.ENDC+ "%.2f,\t%.2f,\t%.2f" % (h_fin,m0-2.5*log10(I0d_fin)+5.*log10(pix2sec),R_disc_min)
			axvline(x=R_disc_min,color='red',linewidth=2,ls='--')
			axvline(x=-R_disc_min,color='red',linewidth=2,ls='--')
			#/**************************************************************************************
			
			
	      
			#**************************************************************************************
			#### Disc breaks
			Rmin = R_disc_min

			# Positive radii:   
			R,Mag = summed_hor_profile(input_image,xc,zc,Radius,Minor_axis,m0,sign='+',min_radius=Rmin)
			MAG_RIGHT = Mag
			R_RIGHT = R

			if inter_mode==True and repeat==2:
			  x_break_right = R_break_right_man
			  y_break_right = Mag_break_right_man
			  def piecewise_linear_fix1(x, k1, k2):
			      return np.piecewise(x, [x < x_break_right], [lambda x:k1*x + y_break_right-k1* x_break_right, lambda x:k2*x +y_break_right-k2* x_break_right])
			  results_pos1, e = optimize.curve_fit(piecewise_linear_fix1, R, Mag,[1,1]) 
			  results_pos = [R_break_right_man,Mag_break_right_man,results_pos1[0],results_pos1[1]]
			else:
			  Res_chi2 = []
			  Res_p = []
			  for rr in R:
			    ind = list(R).index(rr)
			    p , e = optimize.curve_fit(piecewise_linear, R, Mag,[rr,Mag[ind],1,1])
			    chi2 = chi2_func(Mag,piecewise_linear(R, *p))
			    Res_chi2.append(chi2)
			    Res_p.append(p)
			    
			  
			  results_pos = Res_p[Res_chi2.index(min(Res_chi2))]
			  
			Rlim_pos =  int(max(R))  
			rd_pos = np.linspace(Rmin, int(max(R)), 100)
		  
			
			# Negative radii:   
			R,Mag = summed_hor_profile(input_image,xc,zc,Radius,Minor_axis,m0,sign='-',min_radius=Rmin)
			MAG_LEFT = Mag
			R_LEFT = R

			if inter_mode==True and repeat==2:
			  x_break_left = R_break_left_man
			  y_break_left = Mag_break_left_man
			  def piecewise_linear_fix2(x, k1, k2):
			      return np.piecewise(x, [x < x_break_left], [lambda x:k1*x + y_break_left-k1* x_break_left, lambda x:k2*x +y_break_left-k2* x_break_left])
			  results_neg1, e = optimize.curve_fit(piecewise_linear_fix2, R, Mag,[1,1])
			  results_neg = [R_break_left_man,Mag_break_left_man,results_neg1[0],results_neg1[1]]
			else:
			  Res_chi2 = []
			  Res_p = []
			  for rr in R:
			    ind = list(R).index(rr)
			    try:
			      p , e = optimize.curve_fit(piecewise_linear, R, Mag,[rr,Mag[ind],1,1])
			      chi2 = chi2_func(Mag,piecewise_linear(R, *p))
			      Res_chi2.append(chi2)
			      Res_p.append(p)
			    except:
			      zz = 1
			  
			  results_neg = Res_p[Res_chi2.index(min(Res_chi2))]
			  
			Rlim_neg =  int(max(fabs(R)))
			rd_neg = np.linspace(Rmin, int(max(R)), 100)
			plt.ylim(max(Mag),min_mag)
			plt.xlim(-max([Rlim_pos,int(max(R))]),max([Rlim_pos,int(max(R))]))



			rd_gen = np.linspace(-Radius, Radius, 100)
			plt.plot(rd_gen, m0-2.5*log10(disk_edge_r([I0d_fin,h_fin,0.],rd_gen)),color='black')
		
			
			#********************************************************************
			print Rlim_pos,results_pos[0]
			print Rlim_neg,results_neg[0]
			# RESULTS:
			check_pos = 1
			check_neg = 1
			check_all = 1
			if not check_sign([results_pos[2],results_pos[3]]) or rat(Radius,results_pos[0])<1.1:
			    print bcolors.FAIL+'Positive radii: Inner and outer disc profiles have different signs!'+ bcolors.ENDC
			    plt.plot(rd_pos, piecewise_linear(rd_pos, *results_pos),color='blue')
			    axvline(x=results_pos[0],color='blue',linewidth=2,ls='--')
			    check_pos = 0
			else:
			    plt.plot(rd_pos, piecewise_linear(rd_pos, *results_pos),color='red')
			    axvline(x=results_pos[0],color='black',linewidth=2,ls='--')

			if not check_sign([results_neg[2],results_neg[3]]) or rat(Radius,results_neg[0])<1.1:
			    print bcolors.FAIL+'Negative radii: Inner and outer disc profiles have different signs!'+ bcolors.ENDC
			    plt.plot(-rd_neg, piecewise_linear(rd_neg, *results_neg),color='blue')
			    axvline(x=-results_neg[0],color='blue',linewidth=2,ls='--')
			    check_neg = 0
			else:
			    plt.plot(-rd_neg, piecewise_linear(rd_neg, *results_neg),color='red')
			    axvline(x=-results_neg[0],color='black',linewidth=2,ls='--')

			if not check_sign([results_pos[2],results_pos[3],results_neg[2],results_neg[3]]) and check_pos==1 and check_neg==1:
			    print bcolors.FAIL+'All radii: Inner and outer disc profiles have different signs!'+ bcolors.ENDC
			    check_all = 0

			'''
			plt.plot(rd_pos, piecewise_linear(rd_pos, *results_pos),color='red')
			axvline(x=results_pos[0],color='black',linewidth=2,ls='--')    
			plt.plot(-rd_neg, piecewise_linear(rd_neg, *results_neg),color='red')
			axvline(x=-results_neg[0],color='black',linewidth=2,ls='--')    
			'''
			if inter_mode==True:
			  #ds9_show(input_image)
			  fiG = plt.gcf()
			  plt.show()
			  s = raw_input('Are you happy with this [YES/no]?') or "yes"
			  if s!='yes':
			    input1 = raw_input('Please enter left and right radius of the disc:') or str(R_disc_min)+','+str(R_disc_min)
			    input2 =  raw_input('Please enter left break radius of the disc (x,y):') or str(results_neg[0])+','+str(MAG_LEFT[list(R_LEFT).index(find_nearest(results_neg[0],R_LEFT))])
			    input3 =  raw_input('Please enter right break radius of the disc (x,y):') or str(results_pos[0])+','+str(MAG_RIGHT[list(R_RIGHT).index(find_nearest(results_pos[0],R_RIGHT))])
			    
			    R_disc_min_left_man,R_disc_min_right_man = input1.split(',')
			    R_break_left_man,Mag_break_left_man = input2.split(',')
			    R_break_right_man,Mag_break_right_man = input3.split(',')

			    print '\n You entered:'
			    print R_disc_min_left_man,R_disc_min_right_man
			    print R_break_left_man,Mag_break_left_man
			    print R_break_right_man,Mag_break_right_man
			  
			    R_disc_min_left_man = float(R_disc_min_left_man)
			    R_disc_min_right_man = float(R_disc_min_right_man)
			    R_break_left_man = float(R_break_left_man)
			    R_break_right_man = float(R_break_right_man)
			    Mag_break_left_man = float(Mag_break_left_man)
			    Mag_break_right_man = float(Mag_break_right_man)
			    repeat=2
			  else:
			    if check_pos == 0:
			      results_pos[0]=float('nan')
			      results_pos[1]=float('nan')
			      results_pos[2]=float('nan')
			      results_pos[3]=float('nan')
			      
			    if check_neg == 0: 
			      results_neg[0]=float('nan')
			      results_neg[1]=float('nan')
			      results_neg[2]=float('nan')
			      results_neg[3]=float('nan')
			    if check_all == 0:
			      results_pos[0]=float('nan')
			      results_pos[1]=float('nan')
			      results_pos[2]=float('nan')
			      results_pos[3]=float('nan')
			      results_neg[0]=float('nan')
			      results_neg[1]=float('nan')
			      results_neg[2]=float('nan')
			      results_neg[3]=float('nan')
			      
			    R_disc_min_left_man = R_disc_min
			    R_disc_min_right_man = R_disc_min
			    R_break_left_man = results_neg[0]
			    R_break_right_man = results_pos[0]
			    PARS["PARS_MAN"] = [R_disc_min_left_man,R_disc_min_right_man,R_break_left_man,R_break_right_man]
			    repeat = False
			    #plt.clf()
			    #plt.close() 
			else:
			  fiG = plt.gcf()
			  repeat=False

			fiG.savefig('bulge_disc_dom.eps', transparent = False, dpi=300)
			plt.clf()
			plt.close()  
		    
			R_break_pos = results_pos[0]
			R_break_neg = results_neg[0]
			R_break_mean = (R_break_pos+R_break_neg) / 2.
			
			# Inner disc:
			h_in_pos = 1.086/results_pos[2] / 1.3
			h_in_neg = 1.086/results_neg[2] / 1.3
			h_in_mean = (h_in_pos+h_in_neg) / 2.
			m0d_in_pos = results_pos[1]-results_pos[2]*results_pos[1]+5.*log10(pix2sec) + 2.5*log10(1.9)
			m0d_in_neg = results_neg[1]-results_neg[2]*results_neg[1]+5.*log10(pix2sec) + 2.5*log10(1.9)
			m0d_in_mean = m0-2.5*log10((10**(0.4*(m0-m0d_in_pos)) + 10**(0.4*(m0-m0d_in_neg)))/2.)
			
			# Outer disc:
			h_out_pos = 1.086/results_pos[3] / 1.3
			h_out_neg = 1.086/results_neg[3] / 1.3
			h_out_mean = (h_out_pos+h_out_neg) / 2.
			m0d_out_pos = results_pos[1]-results_pos[3]*results_pos[1]+5.*log10(pix2sec) + 2.5*log10(1.9)
			m0d_out_neg = results_neg[1]-results_neg[3]*results_neg[1]+5.*log10(pix2sec) + 2.5*log10(1.9)
			m0d_out_mean = m0-2.5*log10((10**(0.4*(m0-m0d_out_pos)) + 10**(0.4*(m0-m0d_out_neg)))/2.)      

			print bcolors.HEADER+"\n\nINNER AND OUTER DISCS: Results of the fitting (Rmin=%.1f, Rmax=%.1f):" % (Rmin,Radius) + bcolors.ENDC
			print bcolors.OKGREEN+"Break radius [pix]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)" % (R_break_mean,R_break_neg,R_break_pos)
			print bcolors.OKGREEN+"h_in [pix]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)" % (h_in_mean,h_in_neg,h_in_pos)
			print bcolors.OKGREEN+"m0d_in [mag arcsec^-2]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)"  % (m0d_in_mean,m0d_in_neg,m0d_in_pos)
			print bcolors.OKGREEN+"h_out [pix]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)" % (h_out_mean,h_out_neg,h_out_pos)
			print bcolors.OKGREEN+"m0d_out [mag arcsec^-2]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)\n"  % (m0d_out_mean,m0d_out_neg,m0d_out_pos)       
			#/**************************************************************************************
    






























      #### Disc thickness:
      plt.figure(7,figsize=(6, 5))

      plt.xlabel(r' z (pix) ', fontsize=15)
      plt.ylabel(r' $\mu$ (mag pix$^{-2}$) ', fontsize=15)
      plt.gca().invert_yaxis()

      Rmin = R_disc_min

 
      R_gen,Mag_gen = profile_z(input_image,mask_image,pix2sec,m0,Rmin,Radius,zmin,Minor_axis,"0",xc,zc)
      Mag_gen = 10**(0.4*(m0-np.array(Mag_gen)))
      res_gen, e_gen = curve_fit(disk_edge_exp,R_gen,Mag_gen,p0=(max(Mag_gen),max(R_gen)/3.,0.))
      
      #plt.show()
      #exit()
 
      #********************************************************************
      # RESULTS:
     
      z0_gen = res_gen[1]


      print bcolors.HEADER+"\n\nGENERAL THICKNESS OF THE DISC (Rmin=%.1f, Rmax=%.1f):" % (Rmin,Radius) + bcolors.ENDC
      print bcolors.OKGREEN+"z0=2*hz, (exp-law) [pix]:" + bcolors.ENDC+ "  %.3f " % (2.*z0_gen)
      #********************************************************************






































      #### Thin and thick discs:
      #plt.figure(7,figsize=(6, 5))

      #plt.xlabel(r' z (pix) ', fontsize=15)
      #plt.ylabel(r' $\mu$ (mag pix$^{-2}$) ', fontsize=15)
      #plt.gca().invert_yaxis()

      Rmin = R_disc_min
      #zmin = Minor_axis


      # Positive radii:   
      R,Mag = profile_z(input_image,mask_image,pix2sec,m0,Rmin,Radius,zmin,Minor_axis,"+",xc,zc)

      Res_chi2 = []
      Res_p = []
      for rr in R:
	ind = list(R).index(rr)
	p , e = optimize.curve_fit(piecewise_linear, R, Mag,[rr,Mag[ind],1,1])
	chi2 = chi2_func(Mag,piecewise_linear(R, *p))
	Res_chi2.append(chi2)
	Res_p.append(p)

      results_pos_z = Res_p[Res_chi2.index(min(Res_chi2))]
      Rlim_pos =  int(max(R))
      #rd = np.linspace(0, int(max(R)), 100)
      #plt.plot(rd, piecewise_linear(rd, *results_pos_z),color='red')


      # Negative radii:   
      R,Mag = profile_z(input_image,mask_image,pix2sec,m0,Rmin,Radius,zmin,Minor_axis,"-",xc,zc)
   
      Res_chi2 = []
      Res_p = []
      for rr in R:
	ind = list(R).index(rr)
	p , e = optimize.curve_fit(piecewise_linear, R, Mag,[rr,Mag[ind],1,1])
	chi2 = chi2_func(Mag,piecewise_linear(R, *p))
	Res_chi2.append(chi2)
	Res_p.append(p)
 
      results_neg_z = Res_p[Res_chi2.index(min(Res_chi2))]

      #rd = np.linspace(0, int(max(R)), 100)
      #plt.plot(-rd, piecewise_linear(rd, *results_neg_z),color='red')
      
      rd = np.linspace(-int(max(R_gen)), int(max(R_gen)), 100)
      plt.plot(R_gen,m0-2.5*log10(disk_edge_exp(R_gen,*res_gen)),color='black')
      
      plt.xlim(-max([Rlim_pos,int(max(R))]),max([Rlim_pos,int(max(R))]))



 
      #********************************************************************
      # RESULTS:

      #print rat(2.*z0_gen, 1.086/results_pos_z[2]),rat(1.086/results_pos_z[2],1.086/results_pos_z[3])
      #exit()
      #print check_sign([results_pos_z[2],results_pos_z[3]]),rat(1.086/results_pos_z[2],1.086/results_pos_z[3])
      #exit()
      all_bad = 0
      if not check_sign([results_pos_z[2],results_pos_z[3],results_neg_z[2],results_neg_z[3]]) or 1.086/results_pos_z[2]+1.086/results_neg_z[2]>1.086/results_pos_z[3]+1.086/results_neg_z[3]:
	  #print bcolors.FAIL+'All radii: Thin and thick disc profiles have different signs!'+ bcolors.ENDC
	  results_neg_z[0]=float('nan')
	  results_neg_z[1]=float('nan')
	  results_neg_z[2]=float('nan')
	  results_neg_z[3]=float('nan')
	  results_pos_z[0]=float('nan')
	  results_pos_z[1]=float('nan')
	  results_pos_z[2]=float('nan')
	  results_pos_z[3]=float('nan')
	  all_bad = 1

      if not check_sign([results_pos_z[2],results_pos_z[3]]) or rat(1.086/results_pos_z[2],1.086/results_pos_z[3])>10. or all_bad == 1: #or rat(2.*z0_gen, 1.086/results_pos_z[2])<1.2
	  print bcolors.FAIL+'Positive radii: Thin and thick disc profiles have different signs!'+ bcolors.ENDC
	  results_pos_z[0]=float('nan')
	  results_pos_z[1]=float('nan')
	  results_pos_z[2]=float('nan')
	  results_pos_z[3]=float('nan')
      else:
	  rd = np.linspace(0, int(max(R_gen)), 100)
	  plt.plot(rd, piecewise_linear(rd, *results_pos_z),color='red')
 
      if not check_sign([results_neg_z[2],results_neg_z[3]]) or rat(1.086/results_neg_z[2],1.086/results_neg_z[3])>10. or all_bad == 1: #or rat(2.*z0_gen, 1.086/results_neg_z[2])<1.2
	  print bcolors.FAIL+'Negative radii: Thin and thick disc profiles have different signs!'+ bcolors.ENDC
	  results_neg_z[0]=float('nan')
	  results_neg_z[1]=float('nan')
	  results_neg_z[2]=float('nan')
	  results_neg_z[3]=float('nan')
      else:
	  rd = np.linspace(0, int(max(R_gen)), 100)
	  plt.plot(-rd, piecewise_linear(rd, *results_neg_z),color='red')





	  
      z_break_pos = results_pos_z[0]
      z_break_neg = results_neg_z[0]
      z_break_mean = (z_break_pos+z_break_neg) / 2.
      
      z0_thin_pos = 1.086/results_pos_z[2]
      z0_thin_neg = 1.086/results_neg_z[2]
      z0_thin_mean = (z0_thin_pos+z0_thin_neg) / 2.
      
      z0_thick_pos = 1.086/results_pos_z[3]
      z0_thick_neg = 1.086/results_neg_z[3]
      z0_thick_mean = (z0_thick_pos+z0_thick_neg) / 2.
      
      delta_m0d_thin_thick = fabs(mean([results_pos_z[1]-results_pos_z[3]*results_pos_z[1],results_neg_z[1]-results_neg_z[3]*results_neg_z[1]])+5.*log10(pix2sec) - mean([results_pos_z[1]-results_pos_z[2]*results_pos_z[1],results_neg_z[1]-results_neg_z[2]*results_neg_z[1]])+5.*log10(pix2sec))

      print bcolors.HEADER+"\n\nTHIN AND THICK DISCS: Results of the fitting (Rmin=%.1f, Rmax=%.1f):" % (Rmin,Radius) + bcolors.ENDC
      print bcolors.OKGREEN+"Z break radius [pix]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)" % (z_break_mean,z_break_neg,z_break_pos)
      print bcolors.OKGREEN+"z0=2*hz, thin disc (exp-law) [pix]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)" % (z0_thin_mean,z0_thin_neg,z0_thin_pos)
      print bcolors.OKGREEN+"z0=2*hz, thick disc (exp-law) [pix]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)" % (z0_thick_mean,z0_thick_neg,z0_thick_pos)
      print bcolors.OKGREEN+"Delta Z-m0d (mag arcsec^-2):" + bcolors.ENDC+ "  %.3f " % (delta_m0d_thin_thick)


      plt.savefig('thin_thick_discs.eps', transparent = False, dpi=300)
      plt.clf()
      plt.close()  

      #********************************************************************


    
      PARS["INNER_DISC"] = [h_in_mean,h_in_neg,h_in_pos,m0d_in_mean,m0d_in_neg,m0d_in_pos]
      PARS["OUTER_DISC"] = [h_out_mean,h_out_neg,h_out_pos,m0d_out_mean,m0d_out_neg,m0d_out_pos]
      PARS["BREAK_RADIUS"] = [R_break_mean,R_break_neg,R_break_pos]
      PARS["THICKNESS"] = [z0_gen]
      PARS["BREAK_Z"] = [z_break_mean,z_break_neg,z_break_mean]
      PARS["THIN_DISC"] = [z0_thin_mean,z0_thin_neg,z0_thin_pos]
      PARS["THICK_DISC"] = [z0_thick_mean,z0_thick_neg,z0_thick_pos]

      # Differnece between central SB of thin and thick discs:
      PARS["DELTA_M0D"] = [delta_m0d_thin_thick]
      
      crea_reg(xc,zc,R_bulge_ma,R_bulge_mi,R_disc_min,ySize/2.,mean([results_pos[0],results_neg[0]]))
      print 'Done!'
      if inter_mode==False:
	PARS["PARS_MAN"] = [float(nan),float(nan),float(nan),float(nan)]
      return PARS,'thin_thick_discs.eps','bulge_disc_dom.eps',"bulge_disc_boxes.reg"



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Bulge and disc domination analysis + break analysis of the discs")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("maskImage", help="Input mask image (with masked objects which have intensities >0")
    parser.add_argument("magLevel", help="Input the outer isophote level in [mag/arcsec^2]",type=float)
    parser.add_argument("Radius", help="Input the outer isophote major axis in [pix]",type=float)
    parser.add_argument("MinorAxis", help="Input the outer isophote minor axis in [pix]",type=float)
    parser.add_argument("zmin", help="Input the half width of the dust lane in [pix]",type=float)
    parser.add_argument("ZeroPoint", help="Input Zero Point in [mag/arcsec^2]",type=float) 
    parser.add_argument("Scale", help="Input scale in [arcsec/pix]",type=float)
    parser.add_argument("xc", nargs='?', const=0.,help="Input x-coordinate of the galaxy center [pix]",type=float,default=0.)
    parser.add_argument("zc", nargs='?', const=0.,help="Input z-coordinate of the galaxy center [pix]",type=float,default=0.)


    args = parser.parse_args()


    input_image = args.inputImage
    mask_image = args.maskImage
    mag_level = args.magLevel
    Radius = args.Radius
    MinorAxis = args.MinorAxis
    zmin = args.zmin
    m0 = args.ZeroPoint
    pix2sec = args.Scale
    xc = args.xc
    zc = args.zc

    input_image = args.inputImage
    main(input_image,mask_image,mag_level,Radius,MinorAxis,zmin,m0,pix2sec,xc=xc,zc=zc)
















#main('cropped_r.fits',11,120.,5.,28.,0.4,zc=0.,xc=0.)



#main('WHT_g_a_r_norm.fits',3.e-06,21.097,0.25)
#pars = main('Spitzer_1_a_r_norm.fits',3.e-06,21.097,0.25)
#print pars
