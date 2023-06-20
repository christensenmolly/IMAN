#!/usr/bin/env python

from math import sin
from math import cos
from math import pi
from math import modf
from math import sqrt
import numpy as np
import argparse
from types import SimpleNamespace

from astropy.io import fits
from photutils.isophote import IsophoteList
from photutils.isophote import build_ellipse_model
from astropy.stats import sigma_clipped_stats
import math

def ds9_to_np(x):
    '''
    Function to convert ds9 pixel coordinates to numpy format, i.e
    0.5 pix is 0, 1.5 pix is 1
    '''
    return int(round(x)) - 1


def build_model(sma, inten, ell, PA, output_image, input_image=None, rmax=None, center=True, fix_ellipses=True):
    print('Building IRAF/ellipse model image...')
    if input_image is None:
        if rmax is None:
            rmax = max(sma)
    else:
        hdulist = fits.open(input_image)
        data = hdulist[0].data        
        ny,nx = np.shape(data)
            
    if center and input_image is None:
        nx = int(math.ceil(2.*(rmax))) +5
        ny = nx      

        # Determine the center:
        if nx%2==False:
            xc = nx/2. + 0.5
        else:
            xc = nx/2. + 1.

        if ny%2==False:
            yc = ny/2. + 0.5
        else:
            yc = ny/2. + 1.
        
        xc = ds9_to_np(xc)-1
        yc = ds9_to_np(yc)-1
        
        x0 = np.array([xc] * len(sma), float)
        y0 = np.array([yc] * len(sma), float)
    

    if fix_ellipses:
        ell = np.array([0.0] * len(ell), float)
        PA = np.array([0.] * len(PA), float)  
    
    ellipses = IsophoteList([])

    #ellipses.append(SimpleNamespace(x0=10., y0=10., sma=0., intens=inten[0], int_err=inten_err[0],
    #                              eps=ell[0], ellip_err=errell[0], pa=np.radians(PA[0]), pa_err=np.radians(errPA[0]), grad=0., a3=0., b3=0., #a4=0., b4=0., grad_r_error=0., ndata=1, nflag=1, niter=1, stop_code=0, fit_ok=False))

    for k in range(len(sma)-1):
      if sma[k]<=rmax:
        ellipse = SimpleNamespace(x0=x0[k], y0=y0[k], sma=sma[k], intens=inten[k],
                                  eps=0., pa=0., grad=0, a3=0, b3=0, a4=0, b4=0, fit_ok=False)
        ellipses.append(ellipse)    

    model = build_ellipse_model((ny,nx), ellipses)

    fits.PrimaryHDU(data=model).writeto(output_image, overwrite=True)
    print('Done!')

def main(input_image=None, iraf_ell_file='ellipse.txt', output_model='ellipse.fits', center=True, rmax=None, fix_ellipses=True):

    
    
    sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = np.loadtxt(iraf_ell_file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 5, dtype='str')

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
    
    
    if input_image is None:
        if rmax is None:
            rmax = max(sma)
    else:
        hdulist = fits.open(input_image)
        data = hdulist[0].data        
        ny,nx = np.shape(data)
            
    if center and input_image is None:
        nx = int(math.ceil(2.*(rmax))) +5
        ny = nx      

        # Determine the center:
        if nx%2==False:
            xc = nx/2. + 0.5
        else:
            xc = nx/2. + 1.

        if ny%2==False:
            yc = ny/2. + 0.5
        else:
            yc = ny/2. + 1.
        
        xc = ds9_to_np(xc)-1
        yc = ds9_to_np(yc)-1
        
        x0 = np.array([xc] * len(x0), float)
        y0 = np.array([yc] * len(y0), float)
    

    if fix_ellipses:
        ell = np.array([0.0] * len(ell), float)
        PA = np.array([0.] * len(PA), float)  
    
    ellipses = IsophoteList([])

    #ellipses.append(SimpleNamespace(x0=10., y0=10., sma=0., intens=inten[0], int_err=inten_err[0],
    #                              eps=ell[0], ellip_err=errell[0], pa=np.radians(PA[0]), pa_err=np.radians(errPA[0]), grad=0., a3=0., b3=0., #a4=0., b4=0., grad_r_error=0., ndata=1, nflag=1, niter=1, stop_code=0, fit_ok=False))

    for k in range(len(sma)-1):
      if sma[k]<=rmax:
        ellipse = SimpleNamespace(x0=x0[k], y0=y0[k], sma=sma[k], intens=inten[k],
                                  eps=0., pa=0., grad=0, a3=0, b3=0, a4=0, b4=0, fit_ok=False)
        #print(ellipse)
        ellipses.append(ellipse)    
    #print(ny,nx)

    model = build_ellipse_model((ny,nx), ellipses)

    fits.PrimaryHDU(data=model).writeto(output_model, overwrite=True)
    print('Done!')