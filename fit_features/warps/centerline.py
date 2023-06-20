#!/usr/bin/python
# -*- coding:  cp1251 -*-
# python ~/CurrentWork/ImaPrep/IMAN/Warps/warps_new.py cropped_i.fits isophotes.reg 303 96
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
import pyfits
from scipy import interpolate
import argparse


def find_bissect(x,k1,b1,k2,b2):
  A1 = -k1
  B1 = 1.
  C1 = -b1
  
  A2 = -k2
  B2 = 1.
  C2 = -b2
  sqrt1 = sqrt(A1**2+B1**2)
  sqrt2 = sqrt(A2**2+B2**2)
  
  Sol1 = 1. / (B1/sqrt1-B2/sqrt2) * (-x*(A1/sqrt1-A2/sqrt2)-C1/sqrt1+C2/sqrt2) 
  Sol2 = 1. / (B1/sqrt1+B2/sqrt2) * (-x*(A1/sqrt1+A2/sqrt2)-C1/sqrt1-C2/sqrt2) 
  return Sol1,Sol2
  




def Ellipse(a,b,an,x0,y0): #an is the rotational angle
    points=100 #Number of points whicnh needs to construct the elipse
    cos_a,sin_a=cos(an*pi/180),sin(an*pi/180)
    the=linspace(0,2*pi,points)
    #Here goes the general ellpse, x0, y0 is the origin of the ellipse in xy plane
    X=a*cos(the)*cos_a-sin_a*b*sin(the)+x0
    Y=a*cos(the)*sin_a+cos_a*b*sin(the)+y0
    return X,Y

def centerline(x,y):
	      x = list(x)
	      y = list(y)
	      points = zip(x,y)
	      sorted_points = sorted(points)
	      new_x = [point[0] for point in sorted_points]
	      new_y = [point[1] for point in sorted_points]  

	      min_x_ind = x.index(min(new_x))
	      max_x_ind = x.index(max(new_x))

	      x = np.array(new_x)
	      y = np.array(new_y)
	      
	      '''
	      step = 3.

	      for X in x:
		if X

	      ZZ = np.polyfit(x_cen, y_cen, 3.)
	      ff = np.poly1d(ZZ)
	      y_new = ff(x_cen)
	      '''
	      
	      ZZ = np.polyfit(x, y, 1.)
	      ff = np.poly1d(ZZ)
	      y_cen = ff(x)
	      x_cen = x
	      x_cen = list(x_cen)
	      y_cen = list(y_cen)
	      
	      return x_cen,y_cen,ff

x,y=Ellipse(5,8,45.,0,0)
x = list(x)
y = list(y)
x_cen,y_cen,ff = centerline(x,y)
y_cen1,x_cen1,ff1 = centerline(y,x)

x_cen = list(x_cen)
x_cen1 = list(x_cen1)
y_cen = list(y_cen)
y_cen1 = list(y_cen1)
XX = x_cen + x_cen1
YY = y_cen + y_cen1


for k in range(4):
  XXX1,YYY1,F1 = centerline(XX,YY)
  YYY2,XXX2,F2 = centerline(YY,XX)
  YY = YYY1 + YYY2
  XX = XXX1 + XXX2




'''
from scipy.optimize import fsolve
f = ff - ff1
xcen = fsolve(f, 0.)
ycen = f(xcen)
print xcen,ycen

# Find tangents:
dff = sp.misc.derivative(ff, xcen, dx=1e-6)
dff1 = sp.misc.derivative(ff1, xcen, dx=1e-6)

x_cen = np.array(x_cen)
ybis1,ybis2 = find_bissect(x_cen,dff,ycen,dff1,ycen)
plt.plot(x_cen,ybis1,'s',x_cen,ybis2,'s')
'''

plt.plot(x,y,'o',XX,YY,'o')



x = np.array(x)
Y = -x
xlim(-20,20)
ylim(-20,20)
plt.plot(x,Y)

plt.show()
  