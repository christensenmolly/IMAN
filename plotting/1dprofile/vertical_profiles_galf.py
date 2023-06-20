#! /usr/bin/env python

import os
import sys

from pylab import *
import pyfits

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import glob
from matplotlib.pyplot import cm
import re
from scipy.odr.odrpack import *

def get_slice(fitsName, xOrig, yOrig, posang, layer=0):
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
    hdu = pyfits.open(fitsName)
    data = hdu[layer].data
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
            rArray.append(y-yOrig)
            iArray.append(I)
    rArray = np.array(rArray,float)
    iArray = np.array(iArray,float)
    return rArray, iArray

def line(B, x):
    return (1.0857/B[0])*fabs(x-B[2]) + B[1]


def main(input_image,m0,pix2sec): 
      plt.figure(0,figsize=(5, 5))



      plt.xlabel(r' z (arcsec) ', fontsize=15)
      plt.ylabel(r' $\mu$ (mag arcsec$^{-2}$) ', fontsize=15)
      plt.gca().invert_yaxis()

	
      layers = [1,2]
      for layer in layers:
	    hdu = pyfits.open(input_image)
	    data = hdu[layer].data
	    ySize, xSize = data.shape
	    I = []
	    r = []
	    for y in range(0,ySize,1):
	      #I.append(sum(data[y,:]))
	      II = []
	      for x in range(0,xSize,1):
		if hdu[1].data[y,x]>0.:
		  II.append(data[y,x])
	      I.append(sum(II))
	      r.append(y-ySize/2.)
	    r = np.array(r)
	    I = np.array(I)

	    
	    if layer==1:
	      color1 = 'white'
	      color2 = 'white'
	      mecolor = 'black'
	      msize=5
	      plot(r*pix2sec, m0 - 2.5*log10(I)+ 5.*log10(pix2sec),'o',color=color1,markeredgecolor=mecolor,markersize=msize)
	    else:
	      color1 = 'red'
	      color2 = 'red'
	      mecolor = 'red'
	      msize=3
	      plot(r*pix2sec, m0 - 2.5*log10(I)+ 5.*log10(pix2sec),'r-',color=color1,lw=2)
	    
	    
	    if layer==1:
	      ylim(np.max(m0 - 2.5*log10(I[~np.isnan(I)])+ 5.*log10(pix2sec)),np.min(m0 - 2.5*log10(I[~np.isnan(I)])+ 5.*log10(pix2sec))-0.5)

      xlim(-max(r*pix2sec),max(r*pix2sec))
      plt.savefig(input_image.split('.fits')[0] + '_prof_ver.eps', transparent = False, dpi=300)     
      plt.clf()
      plt.close() 

if __name__ == '__main__':
    input_image = str(sys.argv[1])
    m0 = float(sys.argv[2])
    pix2sec = float(sys.argv[3])
    main(input_image,m0,pix2sec)
