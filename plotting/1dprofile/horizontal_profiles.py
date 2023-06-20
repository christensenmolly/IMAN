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


def main(input_image,m0,pix2sec,noise):
      plt.figure(0,figsize=(5, 5))



      plt.xlabel(r' r (arcsec) ', fontsize=15)
      plt.ylabel(r' $\mu$ (mag arcsec$^{-2}$) ', fontsize=15)
      plt.gca().invert_yaxis()

	
      layers = [0]
      for layer in layers:
	    #print image
	    hdu = pyfits.open(input_image)
	    data = hdu[layer].data
	    ySize, xSize = data.shape
	    I = []
	    r = []
	    I_noise = []
	    r_noise = []
	    for x in range(0,xSize,1):
	      SUM = sum(data[:,x])
	      if SUM>3.*noise:
		I.append(SUM)
		r.append(x-xSize/2.)
	      else:
		I_noise.append(SUM)
		r_noise.append(x-xSize/2.)
	    r = np.array(r)
	    I = np.array(I)
	    r_noise = np.array(r_noise)
	    I_noise = np.array(I_noise)
    
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
	    
	    plot(r*pix2sec, m0 - 2.5*log10(I)+ 5.*log10(pix2sec),'o',color=color1,markeredgecolor=mecolor,markersize=msize)
	    plot(r_noise*pix2sec, m0 - 2.5*log10(I_noise)+ 5.*log10(pix2sec),'o',color='white',markeredgecolor='black',markersize=msize)
	    '''
	    r1, i1 = get_slice(input_image, xSize/2., ySize/2., 0., layer=layer)
	    r1 = np.array(r1)
	    i1 = np.array(i1)      
	    plot(r1*pix2sec, m0 - 2.5*log10(i1)+ 5.*log10(pix2sec),'o',color=color2,markeredgecolor=mecolor,markersize=msize)      
	    
	    rr1 = r1*pix2sec
	    mag1 = m0 - 2.5*log10(i1)+ 5.*log10(pix2sec)
	    rr = []
	    mag = []
	    for k in range(len(rr1)):
	      if mag1[k]!=float(inf) and math.isnan(mag1[k])==False:
		#print 'here'
		mag.append(mag1[k])
		rr.append(rr1[k])
	    rr = np.array(rr)
	    mag = np.array(mag)
	    '''

      plt.savefig(input_image.split('.fits')[0] + '_prof_hor.eps', transparent = False, dpi=300)     
      #plt.clf()
      #plt.close()
      plt.show()
      return input_image.split('.fits')[0] + '_prof_hor.eps'
      
if __name__ == '__main__':
    input_image = str(sys.argv[1])
    m0 = float(sys.argv[2])
    pix2sec = float(sys.argv[3])
    main(input_image,m0,pix2sec)