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
import warnings
from scipy import special
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy import optimize
import heapq
from astropy.io import fits as pyfits
from scipy import interpolate
import argparse

#PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
#PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
#sys.path.append(PATH_TO_PACKAGE+'/FindFluxes')
import polygon_reg
import crea_skeleton

#warp_method = 'whole'

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


def chi2_func(ref,mod):
  return sum( (ref-mod)**2)

def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])


def piecewise_linear_double(x, xl, yl, xr, yr, kl, km, kr):
    return np.piecewise(x, [x < xl, (x>=xl) & (x<=xr), x>xr], [lambda x:kl*x + yl-kl*xl, lambda x:km*x + yr-km*xr, lambda x:kr*x + yr-kr*xr])

def disk_edge_exp(z,I0d,z0,z_c):
          #*** For edge-on disk SB in mag/arcsec^2 (along z-axis). n - index of this law. ***
          # B[0] = I0d
          # B[2] = z0 
          # B[3] = z_c
          
          return I0d * exp(-fabs(z-z_c)/(2.*z0))
        
def takeClosest(num,collection):
   return min(collection,key=lambda x:abs(x-num))

def warp_geometry_part(x,y,x0,y0,Radius,x_break_left,x_break_right,y_break_left,y_break_right):
      # For the left part of the curve (to the left from the center):
      if x_break_left==0. and x_break_right==0.:
        x_left = []; y_left = []
        for k in range(len(x)):
          if x[k]<x0:
                x_left.append(x[k])
                y_left.append(y[k])
          
        Res_chi2 = []
        Res_p = []
        Res_k = []
        Res_e = []
        for xx in x_left:
          ind = x_left.index(xx)
          p , e = optimize.curve_fit(piecewise_linear, x_left, y_left,[xx,y_left[ind],1,1])
          chi2 = chi2_func(y_left,piecewise_linear(x_left, *p))
          Res_chi2.append(chi2)
          Res_p.append(p)
          Res_k.append(fabs(p[3]))
          Res_e.append(e)
          
          
        #results_neg = Res_p[Res_chi2.index(min(Res_chi2))]        ### TODO: CHECK
        
        results_neg = Res_p[Res_k.index(min(Res_k))]        ### TODO: CHECK
        errors = Res_e[Res_k.index(min(Res_k))]
        errors = np.sqrt(np.diag(errors))
        
        x_break_left = results_neg[0]
        y_break_left = results_neg[1]
        k_1_left = results_neg[2]
        k_2_left = results_neg[3]
        

        x_break_left_er = errors[0]
        y_break_left_er = errors[1]
        k_1_left_er = errors[2]
        k_2_left_er = errors[3]
  
        
        r_warp_left_er = x_break_left_er
        angle_warp_left_er = degrees(fabs(k_1_left_er/(k_1_left**2 + 1.) - k_2_left_er/(k_2_left**2 + 1.)))


    
    
            
            #plt.plot(x,y,'o')
            #xx = np.linspace(min(x_left),max(x_left),100)
            #yy = piecewise_linear(xx, x_break_left, y_break_left, k_1_left, k_2_left)
            #plt.plot(xx,yy,'o')
            #plt.show()


        # For the right part of the curve (to the right from the center):
        x_right = []; y_right = []
        for k in range(len(x)):
          if x[k]>x0:
                x_right.append(x[k])
                y_right.append(y[k])
          
        Res_chi2 = []
        Res_p = []
        Res_k = []
        Res_e = []
        for xx in x_right:
          ind = x_right.index(xx)
          p , e = optimize.curve_fit(piecewise_linear, x_right, y_right,[xx,y_right[ind],1,1])
          chi2 = chi2_func(y_right,piecewise_linear(x_right, *p))
          Res_chi2.append(chi2)
          Res_p.append(p)
          Res_e.append(e)
          Res_k.append(fabs(p[2]))

        results_pos = Res_p[Res_chi2.index(min(Res_chi2))]        ### TODO: CHECK
        errors = Res_e[Res_chi2.index(min(Res_chi2))]        ### TODO: CHECK
        errors = np.sqrt(np.diag(errors))
        #results_pos = Res_p[Res_k.index(min(Res_k))]        ### TODO: CHECK
        
        x_break_right = results_pos[0]
        y_break_right = results_pos[1]
        k_1_right = results_pos[2]
        k_2_right = results_pos[3]

        x_break_right_er = errors[0]
        y_break_right_er = errors[1]
        k_1_right_er = errors[2]
        k_2_right_er = errors[3]
  
        
        r_warp_right_er = x_break_right_er
        angle_warp_right_er = degrees(fabs(k_1_right_er/(k_1_right**2 + 1.) - k_2_right_er/(k_2_right**2 + 1.)))
    

            #xx = np.linspace(min(x_right),max(x_right),100)
            #yy = piecewise_linear(xx, x_break_right, y_break_right, k_1_right, k_2_right)
            #plt.plot(xx,yy,'o')
            #plt.show()
        
      else:
        x_left = []; y_left = []
        for k in range(len(x)):
          if x[k]<x0:
                x_left.append(x[k])
                y_left.append(y[k])
        
        if y_break_left==0. and y_break_right==0.:
          y_break_left = y_left[x_left.index(takeClosest(x_break_left,x_left))]

        
        def piecewise_linear_fix1(x, k1, k2):
                return np.piecewise(x, [x < x_break_left], [lambda x:k1*x + y_break_left-k1* x_break_left, lambda x:k2*x +y_break_left-k2* x_break_left])
        p , e = optimize.curve_fit(piecewise_linear_fix1, x_left, y_left,[1,1])
        results_neg = p
        errors = np.sqrt(np.diag(e))
        k_1_left = results_neg[0]
        k_2_left = results_neg[1]
        
        k_1_left_er = errors[0]
        k_2_left_er = errors[1]    
    
        angle_warp_left_er = degrees(fabs(k_1_left_er/(k_1_left**2 + 1.) - k_2_left_er/(k_2_left**2 + 1.)))
        

        # For the right part of the curve (to the right from the center):
        x_right = []; y_right = []
        for k in range(len(x)):
          if x[k]>x0:
                x_right.append(x[k])
                y_right.append(y[k])
        
        if y_break_left==0. and y_break_right==0.:
          y_break_right = y_right[x_right.index(takeClosest(x_break_right,x_right))]
          
        def piecewise_linear_fix2(x, k1, k2):
                return np.piecewise(x, [x < x_break_right], [lambda x:k1*x + y_break_right-k1*x_break_right, lambda x:k2*x + y_break_right-k2*x_break_right])
        p , e = optimize.curve_fit(piecewise_linear_fix2, x_right, y_right,[1,1])
        results_pos = p
        errors = np.sqrt(np.diag(e))
        
        k_1_right = results_pos[0]
        k_2_right = results_pos[1]  

        k_1_right_er = errors[0]
        k_2_right_er = errors[1]    
    
        angle_warp_right_er = degrees(fabs(k_1_right_er/(k_1_right**2 + 1.) - k_2_right_er/(k_2_right**2 + 1.)))      

        r_warp_left_er = float(nan)
        r_warp_right_er = float(nan)
      
      # Left warp
      k_left = k_1_left
      r_warp_left = sqrt((x0 - x_break_left)**2 + (y0 - y_break_left)**2)        # Radius at which the warp begins
      x_end_left = min(x_left)                                                        # Left edge of the warp
      y_end_left = k_left*x_end_left + y_break_left-k_left*x_break_left                # Left edge of the warp
      d_warp_left = sqrt((x_break_left-x_end_left)**2 + (y_break_left-y_end_left)**2)        # Length of the warp
      k_middle = k_2_left
      teta_warp_left = fabs(degrees(math.atan2(k_left-k_middle,1.+k_middle*k_left)))        # True warp angle
      kk = (y_end_left-y0) / (x_end_left-x0)
      psi_warp_left = fabs(degrees(math.atan2(kk-k_middle,1.+k_middle*kk)))                # Reshetnokov's angle


      if y_end_left<k_middle*x_end_left+y_break_left-k_middle*x_break_left:
         psi_warp_left = - psi_warp_left        # down bending
         teta_warp_left = - teta_warp_left        # down bending


      # Right warp
      k_right = k_2_right
      x_end_right = max(x_right)
      y_end_right = k_right*x_end_right + y_break_right-k_right*x_break_right
      r_warp_right = sqrt((x0 - x_break_right)**2 + (y0 - y_break_right)**2)        # Radius at which the warp begins

      d_warp_right = sqrt((x_break_right-x_end_right)**2 + (y_break_right-y_end_right)**2)        # Length of the warp
      k_middle = k_1_right
      teta_warp_right = fabs(degrees(math.atan2(k_right-k_middle,1.+k_middle*k_right)))        # True warp angle
      kk = (y_end_right-y0) / (x_end_right-x0)
      psi_warp_right = fabs(degrees(math.atan2(kk-k_middle,1.+k_middle*kk)))                # Reshetnokov's angle


      if y_end_right<k_middle*x_end_right+y_break_right-k_middle*x_break_right:
         psi_warp_right = - psi_warp_right        # down bending
         teta_warp_right = - teta_warp_right        # down bending




    
      left_warp=[x_break_left,y_break_left,k_left,r_warp_left,d_warp_left,psi_warp_left,teta_warp_left]
      right_warp=[x_break_right,y_break_right,k_right,r_warp_right,d_warp_right,psi_warp_right,teta_warp_right]
    
      line_left = [x_end_left,y_end_left,x_break_left,y_break_left]
      line_left_medium = [x_break_left,y_break_left,x0,x0*k_2_left-k_2_left*x_break_left+y_break_left]
      line_right = [x_break_right,y_break_right,x_end_right,y_end_right]
      line_right_medium = [x0,x0*k_1_right-k_1_right*x_break_right+y_break_right,x_break_right,y_break_right]
    
      ERRORS = [r_warp_left_er,angle_warp_left_er,r_warp_right_er,angle_warp_right_er]

      return left_warp,right_warp,line_left,line_left_medium,line_right_medium,line_right,ERRORS
      



def z0_est(INT,x,y,left_warp,right_warp):
        step = 1.
        x_break_left,y_break_left,k_left,r_warp_left,d_warp_left,psi_warp_left,teta_warp_left = left_warp
        x_break_right,y_break_right,k_right,r_warp_right,d_warp_right,psi_warp_right,teta_warp_right = right_warp
      
        # Left:
        x_end_left = min(x)
        x_range_left = arange(x_end_left,x_break_left,step)
        y_range_left = x_range_left*k_left - k_left*x_break_left + y_break_left
        angle = 90. - teta_warp_left
        
        z0_left = []
        x_true_left = []
        y_true_left = []
        for k in range(len(x_range_left)):
            r,I = get_slice(INT, x_range_left[k], y_range_left[k], angle)
        try:
          popt_exp, pcov_exp = curve_fit(disk_edge_exp,r,I,p0=(max(I),mean(fabs(r)),0.))
          z0_left.append(popt_exp[1])
          z_c = popt_exp[2]
          x_true_left.append(x_cen[k]+z_c*cos(radians(angle)))
          y_true_left.append(y_cen[k]+z_c*sin(radians(angle)))
        except:
          zz = 1      
      
      
        # Right:
        x_end_right = max(x)
        x_range_right = arange(x_break_right,x_end_right,step)
        y_end_right = k_right*x_end_right + y_break_right-k_right*x_break_right
        y_range_right = x_range_right*k_right - k_right*x_break_right + y_break_right
        angle = 90. - teta_warp_right

        z0_right = []
        x_true_right = []
        y_true_right = []
        for k in range(len(x_range_right)):
            r,I = get_slice(INT, x_range_right[k], y_range_right[k], angle)
        try:
          popt_exp, pcov_exp = curve_fit(disk_edge_exp,r,I,p0=(max(I),mean(fabs(r)),0.))
          z0_right.append(popt_exp[1])
          z_c = popt_exp[2]
          x_true_right.append(x_cen[k]+z_c*cos(radians(angle)))
          y_true_right.append(y_cen[k]+z_c*sin(radians(angle)))
        except:
          zz = 1   
          
        return z0_left,z0_right




def warp_geometry_whole(x,y,x0,y0,Radius,x_break_left,x_break_right,y_break_left,y_break_right):
      # For the left part of the curve (to the left from the center):
      x=np.array(x)
      y=np.array(y)
      if x_break_left==0. and x_break_right==0.:
        Res_chi2 = []
        Res_p = []
        Res_e = []
        Res_k = []
        X = list(x)
        if y[X.index(min(x))]>y0:
          kll = -1.
        else:
          kll = 1.     

        if y[X.index(max(x))]>y0:
          krr = 1.
        else:
          krr = -1.         

        for r in arange(0,Radius,1.):
          #idx = np.abs(x - x0 + r).argmin()      
          yb_left = y[X.index(takeClosest(x0-r,x))]
          #idx = np.abs(x - x0 - r).argmin()
          yb_right = y[X.index(takeClosest(x0+r,x))]      
          p , e = optimize.curve_fit(piecewise_linear_double, x, y,[x0-r,yb_left,x0+r,yb_right,kll,0.,krr])
          chi2 = chi2_func(y,piecewise_linear_double(x, *p))
          Res_chi2.append(chi2)
          Res_p.append(p)
          Res_e.append(e)
          Res_k.append(fabs(p[3]))
          
          
          
        results = Res_p[Res_chi2.index(min(Res_chi2))]        ### TODO: CHECK   
        errors = Res_e[Res_chi2.index(min(Res_chi2))]
        errors = np.sqrt(np.diag(errors))
        #results = Res_p[Res_k.index(min(Res_k))]                ### TODO: CHECK
        
        x_break_left = results[0]
        y_break_left = results[1]
        x_break_right = results[2]
        y_break_right = results[3]  
        k_left = results[4]
        k_middle = results[5]
        k_right = results[6]
    

        x_break_left_er = errors[0]
        y_break_left_er = errors[1]
        x_break_right_er = errors[2]
        y_break_right_er = errors[3]  
        k_left_er = errors[4]
        k_middle_er = errors[5]
        k_right_er = errors[6]    
        
        r_warp_left_er = x_break_left_er
        r_warp_right_er = x_break_right_er
        angle_warp_left_er = degrees(fabs(k_middle_er/(k_middle**2 + 1.) - k_left_er/(k_left**2 + 1.)))
        angle_warp_right_er = degrees(fabs(k_middle_er/(k_middle**2 + 1.) - k_right_er/(k_right**2 + 1.)))
        
        
        
        
        #print k_left,k_middle,k_right
        #plt.plot(x,y,'o')
        #plt.plot(x,piecewise_linear_double(x, *results),'o')
        #plt.show()

      else: 
        if y_break_left==0. and y_break_right==0.:
                X = list(x)
                y_break_left = y[X.index(takeClosest(x_break_left,x))]
                y_break_right = y[X.index(takeClosest(x_break_right,x))]
          
        def piecewise_linear_double_fix1(x, kl, km, kr):
                return np.piecewise(x, [x < x_break_left, (x>=x_break_left) & (x<=x_break_right), x>x_break_right], [lambda x:kl*x + y_break_left-kl*x_break_left, lambda x:km*x + y_break_right-km*x_break_right, lambda x:kr*x +y_break_right-kr* x_break_right])
          
        p , e = optimize.curve_fit(piecewise_linear_double_fix1, x, y,[1,0,1])
        results = p
        errors = np.sqrt(np.diag(e))
        
        k_left = results[0]
        k_middle = results[1]  
        k_right = results[2]
        
        k_left_er = errors[0]
        k_middle_er = errors[1]
        k_right_er = errors[2]

        r_warp_left_er = float(nan)
        r_warp_right_er = float(nan)
        angle_warp_left_er = degrees(fabs(k_middle_er/(k_middle**2 + 1.) - k_left_er/(k_left**2 + 1.)))
        angle_warp_right_er = degrees(fabs(k_middle_er/(k_middle**2 + 1.) - k_right_er/(k_right**2 + 1.)))
      
      # Left warp
      r_warp_left = sqrt((x0 - x_break_left)**2 + (y0 - y_break_left)**2)        # Radius at which the warp begins
      x_end_left = min(x)                                                        # Left edge of the warp
      y_end_left = k_left*x_end_left + y_break_left-k_left*x_break_left                # Left edge of the warp
      d_warp_left = sqrt((x_break_left-x_end_left)**2 + (y_break_left-y_end_left)**2)        # Length of the warp
      teta_warp_left = fabs(degrees(math.atan2(k_left-k_middle,1.+k_middle*k_left)))        # True warp angle
      kk = (y_end_left-y0) / (x_end_left-x0)
      psi_warp_left = fabs(degrees(math.atan2(kk-k_middle,1.+k_middle*kk)))                # Reshetnokov's angle


      if y_end_left<k_middle*x_end_left+y_break_left-k_middle*x_break_left:
         psi_warp_left = - psi_warp_left        # down bending
         teta_warp_left = - teta_warp_left        # down bending


      # Right warp
      r_warp_right = sqrt((x0 - x_break_right)**2 + (y0 - y_break_right)**2)        # Radius at which the warp begins
      x_end_right = max(x)                                                        # Right edge of the warp
      y_end_right = k_right*x_end_right + y_break_right-k_right*x_break_right                # Right edge of the warp
      d_warp_right = sqrt((x_break_right-x_end_right)**2 + (y_break_right-y_end_right)**2)        # Length of the warp
      teta_warp_right = fabs(degrees(math.atan2(k_right-k_middle,1.+k_middle*k_right)))        # True warp angle
      kk = (y_end_right-y0) / (x_end_right-x0)
      psi_warp_right = fabs(degrees(math.atan2(kk-k_middle,1.+k_middle*kk)))                # Reshetnokov's angle


      if y_end_right<k_middle*x_end_right+y_break_right-k_middle*x_break_right:
         psi_warp_right = - psi_warp_right        # down bending
         teta_warp_right = - teta_warp_right        # down bending




    
      left_warp=[x_break_left,y_break_left,k_left,r_warp_left,d_warp_left,psi_warp_left,teta_warp_left]
      right_warp=[x_break_right,y_break_right,k_right,r_warp_right,d_warp_right,psi_warp_right,teta_warp_right]
    
      line_left = [x_end_left,y_end_left,x_break_left,y_break_left]
      line_left_medium = [x_break_left,y_break_left,x_break_right,y_break_right]
      line_right = [x_break_right,y_break_right,x_end_right,y_end_right]
      line_right_medium = [x_break_left,y_break_left,x_break_right,y_break_right]
    
      ERRORS = [r_warp_left_er,angle_warp_left_er,r_warp_right_er,angle_warp_right_er]

      return left_warp,right_warp,line_left,line_left_medium,line_right_medium,line_right,ERRORS



      
def find_perp(x,y,a):
  # Find perpendicular line (angle) to the tangent in the given point
  # glowingpython.blogspot.be/2013/02/visualizing-tangent.html
  spl = interpolate.splrep(x,y)
  small_x = arange(a-5,a+5)
  fa = interpolate.splev(a,spl,der=0)
  fprime = interpolate.splev(a,spl,der=1)
  #tan = fa+fprime*(small_x-a)
  #plot(a,fa,'om',small_x,tan,'--r')
  fprime_perp = -1./fprime
  angle = degrees(math.atan(fprime_perp)) + 90.
  
  return angle
      

def get_slice(data, xOrig, yOrig, posang, layer=0):
    """
    Function gets slice in fits file along specified line.
    Parameters:
        fitsName -- name of fits file.
        xOrig, yOrig -- coordinates of reference point
        posang [degrees] -- position angle of line. posang=0 -> vertical line (North-South slice).
                            positive posang values are for counterclockwise rotation, i.e. slice
                            with posang=45 is from upper-left corner till buttom right
        layer -- number of hdu in multilayered images
    """
    xOrig, yOrig = yOrig, xOrig

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
        if data[ix, iy]!=0.:
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
        if data[ix, iy]!=0.:    
          rArray.append(y-yOrig)
          iArray.append(I)
    rArray = np.array(rArray,float)
    iArray = np.array(iArray,float)
    return rArray, iArray
      

def main(input_image,iso_reg,x0,y0,warp_method='whole',centerline='poly',x_break_left=0.,x_break_right=0.,y_break_left=0.,y_break_right=0.):
        print(bcolors.OKBLUE+'\n\n************ Warp analysis 2015 ************' + bcolors.ENDC)
        print("Analyzing...")
        PARAMS = {}
        PARAMS["RIGHT_WARP"] = [float(nan),float(nan),float(nan),float(nan),float(nan),float(nan),float(nan),float(nan),float(nan),float(nan)]
        PARAMS["LEFT_WARP"] = [float(nan),float(nan),float(nan),float(nan),float(nan),float(nan),float(nan),float(nan),float(nan),float(nan)]
        PARAMS["WARPS_ERRORS"] = [float(nan),float(nan),float(nan),float(nan)]
        
        hdulist = pyfits.open(input_image)
        inframe = hdulist[0].data
        prihdr = hdulist[0].header
        nx, ny =inframe.shape[1], inframe.shape[0]
        
        # 1. Read region file to download the isophotes
        reg = open(iso_reg,'r')

        X_cen = []; Y_cen = []
        std_delta_y = 0.
        STD_DELTA_Y = []
        for Line in reg:
            x = []
            y = []
            if 'polygon' in Line:
                coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]
                for k in range(0,len(coords)-1,2):
                    x.append(float(coords[k]))
                    y.append(float(coords[k+1]))
            else:
                continue
            min_x_ind = x.index(min(x))
            max_x_ind = x.index(max(x))

            x = np.array(x)
            y = np.array(y)
            
            Radius = max(np.sqrt((x-x0)**2+(y-y0)**2))
            if True:
                if len(STD_DELTA_Y)>2:
                    std_delta_y = mean(STD_DELTA_Y)
                print(Line)
                x_cen,y_cen,Std_delta_y = crea_skeleton.main(input_image,Line,std_delta_y)
                STD_DELTA_Y.append(Std_delta_y)
                print('Skeleton')
                # Picture
                points_reg = input_image.split('.fits')[0]+'_warp_points.reg'
                
                points_ds9 = open(points_reg,'w')
                for k in range(len(x_cen)):
                    #print >>points_ds9, "point(%.1f,%.1f) # point=x color=red width=1" % (x_cen[k],y_cen[k])
                    points_ds9.write("circle(%.1f,%.1f,%.1f) # color=red width=1\n" % (x_cen[k],y_cen[k],1.0))
                points_ds9.close()
                exit()
            else:
                print('Poly')
                ZZ = np.polyfit(x, y, 5)
                ff = np.poly1d(ZZ)
                y_cen = ff(x)
                x_cen = x

            '''
            except:
              # Upper curve
              x_up = []
              y_up = []
              if max_x_ind>min_x_ind:
                  step = -1
              else:
                  step = 1

              for k in range(max_x_ind,min_x_ind,step):
                  x_up.append(x[k])
                  y_up.append(y[k])
              
              points = zip(x_up,y_up)
              sorted_points = sorted(points)
              new_x_up = [point[0] for point in sorted_points]
              new_y_up = [point[1] for point in sorted_points]  
              
              new_x_up = np.array(new_x_up)
              new_y_up = np.array(new_y_up)
              
              f = sp.interpolate.interp1d(new_x_up,new_y_up)
              x_new_up = np.arange(ceil(min(new_x_up)), floor(max(new_x_up)))
              y_new_up = f(x_new_up)



              # Bottom curve
              x_bot = []
              y_bot = []
              
              if max_x_ind>min_x_ind:
                  bot_range = range(0,min_x_ind,1) + range(max_x_ind,len(x),1) 
              else:
                  bot_range = range(0,max_x_ind,1) + range(min_x_ind,len(x),1)
                  
              for k in bot_range:
                  x_bot.append(x[k])
                  y_bot.append(y[k])
              
              points = zip(x_bot,y_bot)
              sorted_points = sorted(points)
              new_x_bot = [point[0] for point in sorted_points]
              new_y_bot = [point[1] for point in sorted_points]  
              
              new_x_bot = np.array(new_x_bot)
              new_y_bot = np.array(new_y_bot)

              
              f = sp.interpolate.interp1d(new_x_bot,new_y_bot)
              x_new_bot = np.arange(ceil(min(new_x_bot)), floor(max(new_x_bot)))
              y_new_bot = f(x_new_bot)

              x_cen = []
              y_cen = []
              
              x_new_bot = list(x_new_bot)
              
              
              for k in range(len(x_new_up)):
                  try:
                      y_cen.append((y_new_bot[x_new_bot.index(x_new_up[k])] + y_new_up[k]) / 2.)
                      x_cen.append(x_new_up[k])
                  except:
                      zz = 1
            '''        
            X_cen.append(x_cen)
            Y_cen.append(y_cen)

        

        min_x = 0.
        max_x = 0.
        X_CEN = []; Y_CEN = []
        for k in range(len(X_cen)):
          for i in range(len(X_cen[k])):
              if X_cen[k][i]<min_x or X_cen[k][i]>max_x:
                X_CEN.append(X_cen[k][i])
                Y_CEN.append(Y_cen[k][i])
          min_x = min(X_cen[k])
          max_x = max(X_cen[k])
        
        X_CEN = np.array(X_CEN)
        Y_CEN = np.array(Y_CEN)
        


        # Another approach: pixels with maximum intensities inside the polygon (isophote)
        INT,coords = polygon_reg.Pixels_inside(input_image,iso_reg,mask_image=None)
        
        x_cen = X_CEN
        y_cen = Y_CEN
        xx_cen = X_CEN
        yy_cen = Y_CEN
        

        if centerline=='poly':
          ZZ = np.polyfit(x_cen, y_cen, 3.)
          ff = np.poly1d(ZZ)
          y_new = ff(x_cen)
        else:
          y_new = y_cen


        XX_CEN = x_cen
        YY_CEN = y_cen
        x_cen = list(x_cen)
        y_new = list(y_new)

        new_x_cen = []
        new_y_new = []
        kkkk = -1
        for i in x_cen:
                kkkk = kkkk + 1
                if i not in new_x_cen:
                    new_x_cen.append(i)
                    new_y_new.append(y_new[kkkk])
        
        
        x_cen = np.array(new_x_cen,float)
        y_new = np.array(new_y_new)

        if warp_method=='whole':
            left_warp,right_warp,line_left,line_left_medium,line_right_medium,line_right,errors = warp_geometry_whole(x_cen,y_new,x0,y0,Radius,x_break_left,x_break_right,y_break_left,y_break_right)
        else:
            left_warp,right_warp,line_left,line_left_medium,line_right_medium,line_right,errors = warp_geometry_part(x_cen,y_new,x0,y0,Radius,x_break_left,x_break_right,y_break_left,y_break_right)
        
        
        
        # Find ratio of luminosities: lum(x>Rwarp) / lum(total)
        lum_ratio_left = sum(INT[:,min(xx_cen):int(left_warp[0])]) / sum(INT)
        lum_ratio_right = sum(INT[:,int(right_warp[0]):max(xx_cen)]) / sum(INT)
  
      
        z0_left,z0_right = z0_est(INT,x_cen,y_new,left_warp,right_warp)


        PARAMS["LEFT_WARP"]=[left_warp[0],left_warp[1],left_warp[2],left_warp[3],left_warp[4],left_warp[5],left_warp[6],lum_ratio_left,mean(z0_left)]

        PARAMS["RIGHT_WARP"]=[right_warp[0],right_warp[1],right_warp[2],right_warp[3],right_warp[4],right_warp[5],right_warp[6],lum_ratio_right,mean(z0_right)]
        
        PARAMS["WARPS_ERRORS"] = errors

        #final results:
        print(bcolors.HEADER+"\n\nResults of the analysis:")
        print(bcolors.OKGREEN+"Warp radius - mean (left,right) [pix]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)" % (mean([left_warp[3],right_warp[3]]),left_warp[3],right_warp[3]))
        print(bcolors.OKGREEN+"Warp length - mean (left,right) [pix]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)" % (mean([left_warp[4],right_warp[4]]),left_warp[4],right_warp[4]))
        print(bcolors.OKGREEN+"Psi warp angle (left,right) [deg]:" + bcolors.ENDC+ "  %.3f,%.3f" % (left_warp[5],right_warp[5]))
        print(bcolors.OKGREEN+"Teta warp angle (left,right) [deg]:" + bcolors.ENDC+ "  %.3f,%.3f" % (left_warp[6],right_warp[6]))
        print(bcolors.OKGREEN+"Warp fraction to the total flux - mean (left,right) [pix]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)" % (mean([lum_ratio_left,lum_ratio_right]),lum_ratio_left,lum_ratio_right))
        print(bcolors.OKGREEN+"z0 of the warp - mean (left,right) [pix]:" + bcolors.ENDC+ "  %.3f (%.3f,%.3f)" % (mean([mean(z0_left),mean(z0_right)]),mean(z0_left),mean(z0_right) ))
        print(bcolors.OKGREEN+"Warps erros [r_warp_left_er,angle_warp_left_er,r_warp_right_er,angle_warp_right_er]:" + bcolors.ENDC+ "  %.2f %.2f\t%.2f %.2f" % (errors[0],errors[1],errors[2],errors[3]))
        
        # Picture
        points_reg = input_image.split('.fits')[0]+'_warp_points.reg'
        lines_reg = input_image.split('.fits')[0]+'_warp_lines.reg'
        warps_png = input_image.split('.fits')[0]+'_warps.png'
        
        points_ds9 = open(points_reg,'w')
        for k in range(len(x_cen)):
            #print >>points_ds9, "point(%.1f,%.1f) # point=x color=red width=1" % (x_cen[k],y_cen[k])
            points_ds9.write("circle(%.1f,%.1f,%.1f) # color=red width=1\n" % (x_cen[k],y_new[k],1.0))
        points_ds9.close()

        
        lines_ds9 = open(lines_reg,'w')
        lines_ds9.write("line(%.1f,%.1f,%.1f,%.1f) # line=0 0 color=blue width=3\n" % (line_left[0],line_left[1],line_left[2],line_left[3]))
        lines_ds9.write("line(%.1f,%.1f,%.1f,%.1f) # line=0 0 color=black width=3\n" % (line_left_medium[0],line_left_medium[1],line_left_medium[2],line_left_medium[3]))
        lines_ds9.write("line(%.1f,%.1f,%.1f,%.1f) # line=0 0 color=black width=3\n" % (line_right_medium[0],line_right_medium[1],line_right_medium[2],line_right_medium[3]))
        lines_ds9.write("line(%.1f,%.1f,%.1f,%.1f) # line=0 0 color=blue width=3\n" % (line_right[0],line_right[1],line_right[2],line_right[3]))

        lines_ds9.close()        


        zoom = nx/(2.*np.max(fabs(x_cen-x0)))
        command1 = "ds9 -geometry %ix%i -zoom %.2f " % (1.01*nx,260+ny,zoom)
        if not np.isnan(float(np.max(inframe)/8.)):
            max_int = np.max(inframe)/8.
        else:
            max_int = inframe[int(y0),int(x0)]
        #command1 = command1 + "%s -scale log -cmap grey -invert -contour yes -contour limits %.1f %.1f -contour smooth 5 -contour nlevels 15 -contour scale log -regions load %s -regions load %s " % (input_image,int_level,max_int,lines_reg,points_reg)
        
        command1 = command1 + "%s -scale pow -scale mode minmax -cmap grey -invert -regions load %s -regions load %s -regions load %s " % (input_image,iso_reg,lines_reg,points_reg)
        
        command1 = command1 + "-colorbar no -saveimage png %s -exit" % (warps_png)
        subprocess.call(command1, shell=True)


        return PARAMS
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Warp analysis")
    parser.add_argument("input_image", help="Input image")
    parser.add_argument("iso_reg", help="Input region file with the outer isophote")
    parser.add_argument("x0", help="Input x-coordinate of the center")
    parser.add_argument("y0", help="Input y-coordinate of the center")
    parser.add_argument("method", nargs='?', const=0., help="Optional: Input method to analyze warps (default 'whole' - analyze at once, otherwise - by left and right parts)",type=str,default='whole')
    parser.add_argument("centerline", nargs='?', const=0., help="Optional: Input method how to find centerline for isophotes (default 'poly' - fast but more precise in case of strong warps; 'skel' - using skeleton approach from skimage (very slow!), useful for stromg smoothed warps.)",type=str,default='whole')
    parser.add_argument("x_break_left", nargs='?', const=0., help="Optional: Input x-coordinate of the left warp (start point)",type=float,default=0.)
    parser.add_argument("x_break_right", nargs='?', const=0., help="Optional: Input x-coordinate of the right warp (start point)",type=float,default=0.)    
    parser.add_argument("y_break_left", nargs='?', const=0., help="Optional: Input y-coordinate of the left warp (start point)",type=float,default=0.)
    parser.add_argument("y_break_right", nargs='?', const=0., help="Optional: Input y-coordinate of the right warp (start point)",type=float,default=0.)      
    
    args = parser.parse_args()
    
    input_image = args.input_image
    iso_reg = args.iso_reg
    x0 = float(args.x0)
    y0 = float(args.y0)
    method = args.method
    centerline = args.centerline
    x_break_left = float(args.x_break_left)
    x_break_right = float(args.x_break_right)
    y_break_left = float(args.y_break_left)
    y_break_right = float(args.y_break_right)
  
    main(input_image,iso_reg,x0,y0,method,centerline,x_break_left=x_break_left,x_break_right=x_break_right,y_break_left=y_break_left,y_break_right=y_break_right)

'''
xmin = 35.5
ymin = 15.5
xmax = 170
ymax = 57.5


inner_level=21.5; outer_level=25.5
# 1,g
input_image= '1_g_i.fits'
m0 = 28.55366

import crea_iso
import plot_iso_map_with_warps
x0 = 102.66
y0 = 35.23
pix2sec = 0.396
name = 'SDSS J140639.64+272242.4, g'
lines_reg = '1_g_i_warp_lines.reg'
points_reg = '1_g_i_warp_points.reg'
iso_reg = crea_iso.main(input_image,m0,pix2sec,x0,y0,inner_level,outer_level)
main(input_image,iso_reg,x0,y0)#,method,centerline,x_break_left=x_break_left,x_break_right=x_break_right,y_break_left=y_break_left,y_break_right=y_break_right)
plot_iso_map_with_warps.main(name,input_image,iso_reg,lines_reg,points_reg,'N',x0=xmin,x1=xmax,y0=ymin,y1=ymax):
'''

'''
inner_level=19.5; outer_level=24.5
# 1,r
input_image= '1_r_i.fits'
m0 = 27.32076

import crea_iso
import plot_iso_map_with_warps
x0 = 102.66
y0 = 35.23
pix2sec = 0.396
name = 'r'
lines_reg = '1_r_i_warp_lines.reg'
points_reg = '1_r_i_warp_points.reg'
iso_reg = crea_iso.main(input_image,m0,pix2sec,x0,y0,inner_level,outer_level)
main(input_image,iso_reg,x0,y0)#,method,centerline,x_break_left=x_break_left,x_break_right=x_break_right,y_break_left=y_break_left,y_break_right=y_break_right)
plot_iso_map_with_warps.main(name,input_image,iso_reg,lines_reg,points_reg,'N',x0=xmin,x1=xmax,y0=ymin,y1=ymax):
'''
'''
inner_level=19.5; outer_level=24.5
# 1,i
input_image= '1_I_i.fits'
m0 = 27.320877

import crea_iso
import plot_iso_map_with_warps
x0 = 102.66
y0 = 35.23
pix2sec = 0.396
name = 'i'
lines_reg = '1_I_i_warp_lines.reg'
points_reg = '1_I_i_warp_points.reg'
iso_reg = crea_iso.main(input_image,m0,pix2sec,x0,y0,inner_level,outer_level)
main(input_image,iso_reg,x0,y0)#,method,centerline,x_break_left=x_break_left,x_break_right=x_break_right,y_break_left=y_break_left,y_break_right=y_break_right)
plot_iso_map_with_warps.main(name,input_image,iso_reg,lines_reg,points_reg,'N',x0=xmin,x1=xmax,y0=ymin,y1=ymax):
'''

'''
xmin = 57.7
ymin = 15.5
xmax = 162.5
ymax = 54.5


inner_level=21.5; outer_level=25.0
# 2,g
input_image= '1_g_i.fits'
m0 = 28.453256

import crea_iso
import plot_iso_map_with_warps
x0 = 107.2
y0 = 35.6
pix2sec = 0.396
name = 'SDSS J153538.63+464229.5, g'
lines_reg = '1_g_i_warp_lines.reg'
points_reg = '1_g_i_warp_points.reg'
iso_reg = crea_iso.main(input_image,m0,pix2sec,x0,y0,inner_level,outer_level)
main(input_image,iso_reg,x0,y0,x_break_left=94.5,x_break_right=127.3,y_break_left=36.5,y_break_right=34.3)#,method,centerline,x_break_left=x_break_left,x_break_right=x_break_right,y_break_left=y_break_left,y_break_right=y_break_right)
plot_iso_map_with_warps.main(name,input_image,iso_reg,lines_reg,points_reg,'N',x0=xmin,x1=xmax,y0=ymin,y1=ymax):

'''


'''
inner_level=21.5; outer_level=25.0
# 2,r
input_image= '1_r_i.fits'
m0 = 28.2335

import crea_iso
import plot_iso_map_with_warps
x0 = 107.2
y0 = 35.6
pix2sec = 0.396
name = 'r'
lines_reg = '1_r_i_warp_lines.reg'
points_reg = '1_r_i_warp_points.reg'
iso_reg = crea_iso.main(input_image,m0,pix2sec,x0,y0,inner_level,outer_level)
main(input_image,iso_reg,x0,y0,x_break_left=94.5,x_break_right=127.3,y_break_left=36.5,y_break_right=34.3)#,method,centerline,x_break_left=x_break_left,x_break_right=x_break_right,y_break_left=y_break_left,y_break_right=y_break_right)
plot_iso_map_with_warps.main(name,input_image,iso_reg,lines_reg,points_reg,'N',x0=xmin,x1=xmax,y0=ymin,y1=ymax):
'''
'''
inner_level=21.5; outer_level=25.0
# 2,i
input_image= '1_I_i.fits'
m0 = 27.926155

import crea_iso
import plot_iso_map_with_warps
x0 = 107.2
y0 = 35.6
pix2sec = 0.396
name = 'i'
lines_reg = '1_I_i_warp_lines.reg'
points_reg = '1_I_i_warp_points.reg'
iso_reg = crea_iso.main(input_image,m0,pix2sec,x0,y0,inner_level,outer_level)
main(input_image,iso_reg,x0,y0,x_break_left=94.5,x_break_right=127.3,y_break_left=36.5,y_break_right=34.3)#,method,centerline,x_break_left=x_break_left,x_break_right=x_break_right,y_break_left=y_break_left,y_break_right=y_break_right)
plot_iso_map_with_warps.main(name,input_image,iso_reg,lines_reg,points_reg,'N',x0=xmin,x1=xmax,y0=ymin,y1=ymax):
'''







#x,y = loadtxt('/home/mosenkov/CurrentWork/BTA_warped/interpolated_images/1/res_i/coords.txt', usecols=[0,1],dtype=float, unpack=True,skiprows=0,delimiter='\t')
#warp_geometry(x,y,300,100,167.219,x_break_left=0.,x_break_right=0.)
