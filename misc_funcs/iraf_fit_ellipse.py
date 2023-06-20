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
import argparse

import pyfits

import iraf_ell_plots
tmp_out = sys.stdout

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
        
        

def ring_median_filter(galaxy_image,sma,ellip,x0,y0,B4):
  
  hdulist = pyfits.open(galaxy_image, do_not_scale_image_data=True)
  img = hdulist[0].data


  (dimy,dimx) = img.shape
  shutil.copy(galaxy_image,'galaxy_med.fits') 
  hdulist1 = pyfits.open('galaxy_med.fits', do_not_scale_image_data=True, mode='update')
  img1 = hdulist1[0].data  

  for k in range(dimx):
    for i in range(dimy):
      img1[i,k]=0.  
  
  q = 1. - ellip
  smb = sma*q
  xc = median(x0)
  yc = median(y0)
  C0 = -183.1*B4**3 + 66.26*B4**2 - 12.67*B4 + 0.02504
  #print c0
  
  for p in range(0,len(sma)):
    if p<len(sma)-1:
      #print p,len(sma)
      c0 = (C0[p+1]+C0[p])/2.
      Sma_out = sma[p+1]
      Smb_out = smb[p+1]
      Sma_in = sma[p]
      Smb_in = smb[p]
      Int = []
      II = []
      KK = []
      for x in arange(xc-Sma_out,xc+Sma_out,1):
	for y in arange(yc-Smb_out,yc+Smb_out,1):
	  if (fabs(x-xc))**(c0+2.)/Sma_out**(c0+2.) + (fabs(y-yc))**(c0+2.)/Smb_out**(c0+2.)<=1 and (fabs(x-xc))**(c0+2.)/Sma_in**(c0+2.) + (fabs(y-yc))**(c0+2.)/Smb_in**(c0+2.)>1:
	    K = ceil(x)-1
	    I = ceil(y)-1
	    Int.append(img[I,K])
	    II.append(I)
	    KK.append(K)
	    
      Median_Int = mean(Int)
      print Median_Int
      for k in range(len(II)):
	X = KK[k]
	Y = II[k]
	img1[Y,X] = fabs(img[Y,X] - Median_Int)
    

      

      
  hdulist1.flush()








def main(input_image,xc,yc,m0,pix2sec,minsma=0.,maxsma=0.,step=1.,outp_format='eps'):
    print bcolors.OKBLUE+'\n\n************ IRAF ellipse/fitting analysis (c) Mosenkov A. 2015 ************' + bcolors.ENDC
    print "Analyzing..."

    hdulist = pyfits.open(input_image, do_not_scale_image_data=True)
    scidata = hdulist[0].data
    ny,nx = np.shape(scidata)

    if maxsma==0.:
      maxsma = max([xc,nx-xc,yc,ny-yc])

    output_image = input_image.split('/')[-1].split('.fits')[0]+'_iraf_ell.'+outp_format

    iraf_ell_plots.main_ell(input_image,output_image,xc,yc,step,minsma,maxsma,m0,pix2sec,outp_format)
    try:
        os.remove('ell.cl')
    except:
        z=1

    sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = loadtxt('ellipse.txt', usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 6, dtype='str')
    shutil.move('ellipse.txt',input_image.split('/')[-1].split('.fits')[0]+'_iraf_ell.txt')

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
    
    #return sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4
    return output_image,input_image.split('/')[-1].split('.fits')[0]+'_iraf_ell.txt'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="IRAF/ellipse analysis")
    parser.add_argument("input_image", help="Input image")
    parser.add_argument("x0", help="Input x-coordinate of the center")
    parser.add_argument("y0", help="Input y-coordinate of the center")        
    parser.add_argument("m0", help="Input Zero Point in [mag/arcsec^2]")         
    parser.add_argument("Scale", help="Input scale in [arcsec/pix]")
    parser.add_argument("minsma", nargs='?', const=0., help="Input minimum radius where to start fitting (in [pix])",type=float,default=0.)
    parser.add_argument("maxsma", nargs='?', const=0., help="Input maximum radius where to end fitting (in [pix])",type=float,default=0.)
    parser.add_argument("step", nargs='?', const=0., help="Input step in [pix]",type=float,default=1.)
    parser.add_argument("--o", nargs='?', const=0., help="Format of the output files: eps, png, jpg, pdf",type=str,default='eps')
    
    
    args = parser.parse_args()
    
    input_image = args.input_image
    x0 = float(args.x0)
    y0 = float(args.y0)
    m0 = float(args.m0)
    pix2sec = float(args.Scale)
    minsma = args.minsma
    maxsma = args.maxsma
    step = args.step
    outp_format = args.o
    
    main(input_image,x0,y0,m0,pix2sec,minsma=minsma,maxsma=maxsma,step=step,outp_format=outp_format)

#main('cropped_z.fits',225,70,28.,0.4,3.25,2.)