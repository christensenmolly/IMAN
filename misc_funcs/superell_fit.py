#!/usr/bin/python
# -*- coding:  cp1251 -*-

#*** Common modules ***
import random as random_number
import sys
import math
import numpy as np
from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.gridspec as gridspec
from matplotlib.ticker import NullFormatter
from numpy import *
from pylab import *
import os
import shutil
import subprocess
from os.path import exists
import re
import random
import pylab
from scipy.odr.odrpack import *
from scipy import special
from scipy.optimize import fsolve
from scipy.optimize import fmin
from scipy import signal
#import peakdetect
import pyfits
from numpy.linalg import eig, inv
tmp_out = sys.stdout

def fitEllipse(x,y):
    #http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a

def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])


def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))


def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])

def ellipse_angle_of_rotation2( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    if b == 0:
        if a < c:
            return 0
        else:
            return np.pi/2
    else:
        if a < c:
            return np.arctan(2*b/(a-c))/2
        else:
            return np.pi/2 + np.arctan(2*b/(a-c))/2  



def superell1(a,b,h,k,theta,C0,delta,rand_coord):
  C0 = C0 + 2.
  #http://users.cs.cf.ac.uk/Paul.Rosin/resources/papers/superellipse2.pdf
  #http://en.wikipedia.org/wiki/Superellipse
  #http://dl.acm.org/citation.cfm?id=1157195
  #http://ac.els-cdn.com/S1474034613000281/1-s2.0-S1474034613000281-main.pdf?_tid=963e0fd8-2795-11e4-b19a-00000aab0f6c&acdnat=1408448595_bd3343407e99448a7af09ca46f6aa489
  x1 = arange(h-a,h+a,delta)
  y1 = k + b*(1.-(fabs((x1-h)/a))**C0)**(1./C0)
  x2 = arange(h+a-delta,h-a,-delta)
  y2 = k - b*(1.-(fabs((x2-h)/a))**C0)**(1./C0)
  x = np.concatenate((x1,x2))
  y = np.concatenate((y1,y2))   
  
  xx = []
  yy = []
  for i in range(len(x)):      
    xx.append(((x[i]-h)*cos(radians(theta))+(y[i]-k)*sin(radians(theta))+h)*(1.+random.uniform(-rand_coord,rand_coord)))
    yy.append(((y[i]-k)*cos(radians(theta))-(x[i]-h)*sin(radians(theta))+k)*(1.+random.uniform(-rand_coord,rand_coord)))

  return xx,yy


def fitting_iso(x,y,out):
      f = open(r"iso_in.txt", "w")			# Coordinates of the pixels of a given isophote
      sys.stdout = f
      print "pixel"
      for i in range(len(x)):
	      print "%i %i" % (int(round(x[i])), int(round(y[i]))) #(int(round(x[i]+0.5)), int(round(y[i]+0.5)))
      print "-1 -1"
      sys.stdout = tmp_out
      f.close()
      
      # Running Rosin code:
      subprocess.call("se_epsilon3 iso_in.txt rosin_out.dat", shell=True)     
      xc_se,yc_se,a_se,b_se,c0_se,theta_se = loadtxt('rosin_out.dat', usecols=[0,1,2,3,4,5], unpack=True, skiprows=0)
      theta_se = 2.*math.pi-theta_se
      if max([a_se,b_se])==b_se:
	theta_se = math.pi/2. - theta_se
	BB = min([a_se,b_se])
	a_se = max([a_se,b_se])
	b_se = BB
      if out==1:  
	print "=================================="
	print "Results of SuperEllipse (Rosin):"
	print "xc,yc: %i %i" % (xc_se,yc_se)
	print "a,b: %.2f %.2f" % (a_se,b_se)
	print "c0: %.2f" % (c0_se)
	print "theta: %.2f" % (degrees(theta_se))
	print "==================================\n"  
      
 
      return a_se,b_se,c0_se,xc_se,yc_se
 

def main(x_iso,y_iso,xc,yc,out):
	  #print "\n==============================================="  
	  #print "================ISOPHOTE FITTING================"

	      
	  x_iso1 = []; x_iso2 = []; y_iso1 = []; y_iso2 = []
	  x_iso3 = []; x_iso3 = []; y_iso3 = []; y_iso4 = []	  
	  for k in range(len(x_iso)):
	    if y_iso[k]>yc:
	      x_iso1.append(x_iso[k])
	      y_iso1.append(y_iso[k])
	  x_iso2 = sorted(x_iso1)
	  for k in range(len(x_iso2)):
	    y_iso2.append(y_iso1[x_iso1.index(x_iso2[k])])
	    
	  for k in range(len(x_iso)):
	    if y_iso[k]<=yc:
	      x_iso3.append(x_iso[k])
	      y_iso3.append(y_iso[k])
	  x_iso4 = sorted(x_iso3, reverse=True)
	  for k in range(len(x_iso4)):
	    y_iso4.append(y_iso3[x_iso3.index(x_iso4[k])])

	  xiso = np.concatenate((x_iso2,x_iso4))
	  yiso = np.concatenate((y_iso2,y_iso4))   
  
	  a_se,b_se,c0_se,xc_se,yc_se = fitting_iso(xiso,yiso,out)
	  '''
	  a_se1,b_se1,c0_se1,xc_se1,yc_se1 = fitting_iso(x_iso2,y_iso2,out)
	  plt.plot(x_iso2,y_iso2,'x')
	  plt.show()  
	  
	  
	  a_se2,b_se2,c0_se2,xc_se2,yc_se2 = fitting_iso(x_iso4,y_iso4,out)
	  a_se = mean([a_se1,a_se2])
	  b_se = mean([b_se1,b_se2])
	  c0_se = mean([c0_se1,c0_se2])
	  xc_se = mean([xc_se1,xc_se2])
	  yc_se = mean([yc_se1,yc_se2])
	  '''
	  
	  a = fitEllipse(xiso,yiso)
	  center = ellipse_center(a)
	  phi = ellipse_angle_of_rotation(a)
	  axes = ellipse_axis_length(a)

	  print("center = ",  center)
	  print("angle of rotation = ",  phi)
	  print("axes = ", axes)
	  
	  os.remove('iso_in.txt')
	  os.remove('rosin_out.dat')
	  #a, b = axes
	  #R = np.arange(0,2.*np.pi, 0.01)
	  #xx = center[0] + a*np.cos(R)*np.cos(phi) - b*np.sin(R)*np.sin(phi)
	  #yy = center[1] + a*np.cos(R)*np.sin(phi) + b*np.sin(R)*np.cos(phi)
	  #plt.plot(xiso,yiso)
	  #plt.plot(xx,yy, color = 'red')
	  #plt.show()
	  #return a_se,b_se,c0_se,xc_se,yc_se
	  return axes[0],axes[1],c0_se,center[0],center[1]
      

  


  

