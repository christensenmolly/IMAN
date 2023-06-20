#!/usr/bin/python

import sys
import math
import numpy as np
from scipy import stats
import scipy as sp
#import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
#import matplotlib.patches as patches
#import matplotlib.path as path
#from matplotlib.ticker import NullFormatter
from numpy import *
#from pylab import *
import os
import shutil
import subprocess
import random
#import re
#import glob

from astropy.io import fits as pyfits
#from astropy.modeling import models, fitting
#import itertools
from scipy import ndimage
from astropy.stats import sigma_clipped_stats
from photutils import centroid_com, centroid_1dg, centroid_2dg
import argparse
import shapely
from shapely.geometry import Point


LOCAL_DIR = "/imp/psf"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))
sys.path.append(os.path.join(IMAN_DIR, 'Ellipse_photometry'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/rotate'))
sys.path.append(os.path.join(IMAN_DIR, 'detect_objects'))
sys.path.append(os.path.join(IMAN_DIR, 'sextractor'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/misc'))

import crop_image
import rebin_image
import run_SExtractor
import rotate_image
import azimProfile
import norm_image



def find_borders(data, xc, yc, std_d=None):
    k = 3
    N = 0
    Npix = 10
    go_further=True
    while go_further==True:
        try:
            I1 = data[yc-k,xc-k:xc+k+1]
            I2 = data[yc-k:yc+k+1,xc-k]
            I3 = data[yc+k,xc-k:xc+k+1]
            I4 = data[yc-k:yc+k+1,xc+k]
            I = I1+I2+I3+I4
            mean, median, std = sigma_clipped_stats(I, sigma=3.0)
            if std_d==None:
              std_d = std
            k = k + 1
            if median<=std_d:
              N = N + 1
            if N>=Npix:
              go_further=False        
        except:
            go_further=False

    return k-1,std_d




def center_psf(input_file, R_find=True, R=None, output_file=None, x_l=None, y_l=None, x_r=None, y_r=None, Std=None, verbosity=True):
        hdulist1 = pyfits.open(input_file)
        data = hdulist1[0].data
        ny,nx = np.shape(data)

        if x_l!=None and y_l!=None and x_r!=None and y_r!=None:
            xc,yc = centroid_2dg(data[y_l:y_r+1,x_l:x_r+1])
            xc = int(x_l+xc+0.5)
            yc = int(y_l+yc+0.5)
        else:
            xc,yc = centroid_2dg(data)  # Analyzing the whole image, slow!
            xc = int(xc+0.5)
            yc = int(yc+0.5)

        if verbosity: print('Coordinates of the center are:')
        if verbosity: print('\t DS9 (count from 1): %i, %i' % (xc+1,yc+1))
        if verbosity: print('\t Numpy (count from 0): %i, %i' % (xc,yc))

        

        if x_l!=None and y_l!=None and x_r!=None and y_r!=None:
            R = min([xc-x_l,x_r-xc,yc-y_l,y_r-yc])
        elif x_l==None and y_l==None and x_r==None and y_r==None and R==None:
            R = min([xc,nx-1-xc,yc,ny-1-yc])
        else:
            R = int(R)

        if R_find==True:
            R,std = find_borders(data,xc,yc)
        else:
          std = 0.

        if Std!=None:
          std = Std
        
        while (yc-R<0 or yc+R+1>ny-1 or xc-R<0 or xc+R+1>nx-1):
            if verbosity: print('WARNING: Radius %i is larger than the borders of the input image!' % (R))
            #R = int(ceil(min([yc,ny-yc,xc,nx-xc])))
            R = R-1

        new_data = data[yc-R:yc+R+1,xc-R:xc+R+1]
        ny1,nx1 = np.shape(new_data)

          
        
        hdu = pyfits.PrimaryHDU(new_data)
        if output_file==None:
            output_file = image_file.split('/')[-1].split('.fits')[0]+'_crop.fits'
        hdu.writeto(output_file,clobber=True)

        
        # Mask contaminents:
        run_SExtractor.call_SE(output_file, snr=3., min_pix=5)
        
        hdulist_segm = pyfits.open('segm.fits')
        data_segm = hdulist_segm[0].data

        hdulist_data = pyfits.open(output_file)
        data = hdulist_data[0].data
        ny,nx = np.shape(data)
        I_c = data_segm[int(ny/2.),int(nx/2.)]
        for k in range(ny):
            for i in range(nx):
                if (data_segm[k,i]>0 and data_segm[k,i]!=I_c) or data[k,i]==0.:  # to replace 0 to noise
                    data[k,i] = random.gauss(0.0, std)

        hdu = pyfits.PrimaryHDU(data)
        hdu.writeto(output_file,clobber=True)
        os.remove('segm.fits')

        if verbosity: print('Done!')


    



def backgr_around_star(scidata, xc, yc, Rin, Rout=None):
  nx, ny =scidata.shape[1], scidata.shape[0]
  if Rout==None:
    Rin = int(ceil(1.1*Rin)) # NOTE: 1.1
    Rout = Rin + max([5,int(ceil(1.3*Rin))]) # NOTE: 5 or 1.3
  else:
    Rin = int(ceil(1.1*Rin)) # NOTE: 1.1
    Rout = int(ceil(Rout))    
  backgr = []
  for y in range(int(yc)-Rout,int(yc)+Rout,1):
    for x in range(int(xc)-Rout,int(xc)+Rout,1):
      if (x-int(xc))**2+(y-int(yc))**2>Rin**2 and (x-int(xc))**2+(y-int(yc))**2<=Rout**2 and x>=0 and y>=0 and x<nx and y<ny:
        backgr.append(scidata[y,x])
  mean, median, std = sigma_clipped_stats(backgr, sigma=3.0)
  return median,std



def sky_subtration(input_image, output_image, sky_level):
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    header = hdulist[0].header
    outData = data - sky_level
    outHDU = pyfits.PrimaryHDU(outData, header=header)
    outHDU.writeto(output_image, clobber=True)


def crop_star(input_image, region_file, output_star_image=None, star_number=None, factor=1, PA=0., verbosity=True, azim_model=True, norm=True):
        if star_number==-1:
            p = subprocess.Popen(["ds9",input_image,"-scale","histequ","-invert","-regions","load",region_file])
            star_number = str(1)
            star_number = str(raw_input('Please choose the star to crop (%s): ' % (star_number))) or star_number
            star_number = int(star_number)

        ## Read in the input image
        hdulist = pyfits.open(input_image)
        header1 = hdulist[0].header
        scidata = hdulist[0].data

        ff = open(region_file, "r")
        lines = ff.readlines()
        
        N = 0
        stars = []; stars_circles = []; star_numbers = []
        for k in range(len(lines)):
           if 'circle(' in lines[k] or 'ellipse(' in lines[k] or 'annulus(' in lines[k]:
                N = N + 1
                star = lines[k].split(',')
                xc_star = float(star[0].split('(')[1])
                yc_star = float(star[1])
                if 'circle(' in lines[k]:
                    a_star_in = float(star[2].split(')')[0])
                    a_star_out = 1.25 * a_star_in
                else:
                    a_star_in = float(star[2])
                    a_star_out = float(lines[k].split(',')[4])

                if 'text={' in lines[k]:
                    N = int(lines[k].split('text={')[1].split('}')[0])
                
                if star_number is not None:
                    if N==star_number:
                        XC_star = xc_star
                        YC_star = yc_star
                        A_star_in = a_star_in
                        A_star_out = a_star_out
                        break
                else:
                    stars.append([xc_star,yc_star,a_star_in,a_star_out,N])
                    star_numbers.append(N)
                    stars_circles.append(Point(xc_star,yc_star).buffer(a_star_out))

        if star_number is None:
            good_stars = []
            for k in range(len(stars_circles)):
                intersect = False
                for i in range(len(stars_circles)):
                    if i!=k and stars_circles[k].intersects(stars_circles[i]):
                        if verbosity: print('Star #%i intersects others!' % (star_numbers[k]))
                        intersect = True
                        break
                if intersect==False:
                    good_stars.append(k)
            
            if good_stars == []:
                if verbosity: print('All selected stars intersect others! Please choose the best one by hand!')
                p = subprocess.Popen(["ds9",input_image,"-scale","histequ","-invert","-regions","load",region_file])
                star_number = str(1)
                star_number = str(input('Please choose the star to crop (%s): ' % (star_number))) or star_number
                star_number = int(star_number)
                N = star_number
                try:
                    [XC_star,YC_star,A_star_in,A_star_out,N] = stars[star_numbers.index(star_number)]
                except:
                    if verbosity: print('Star with number %s is not found in the region file. Exiting.')
                    exit()
            else:
                I_max = 0; A_max = 0
                for k in good_stars: ### TODO:!!!!!!!
                    [xc_star,yc_star,a_star_in,a_star_out,N] = stars[k]
                    '''
                    i_max = np.max(scidata[int(yc_star-a_star_out):int(yc_star+a_star_out),int(xc_star-a_star_out):int(xc_star+a_star_out)])
                    if i_max>I_max:
                        I_max = i_max
                        best_star = k
                    '''
                    a_max = a_star_out
                    if a_max>A_max:
                        A_max = a_max
                        best_star = k                        
                        
                [XC_star,YC_star,A_star_in,A_star_out,N] = stars[best_star]

        if verbosity: print('The star #%i at %i,%i was chosen as the best PSF star.' % (N, int(XC_star), int(YC_star)) )
        
        ## We now work with the best PSF star
        ## Cut out the image with the PSF star
        xc_best = int(XC_star+0.5) #http://wcs2kml.googlecode.com/svn-history/r25/trunk/python/wcslib.py
        yc_best = int(YC_star+0.5)

        R_best = int(ceil(A_star_out)) * 1.3

        '''
        if PA!=0.:
            while fabs(PA)>90.:
                PA = fabs(PA) - 90.
            Coeff = sqrt(2.)*fabs(cos(radians(45.-PA)))      
        else:
            Coeff = sqrt(2.)
        '''
                  
        xmin = xc_best - R_best + 1
        xmax = xc_best + R_best
        ymin = yc_best - R_best + 1
        ymax = yc_best + R_best
        
        '''
        # TODO: Check below:
        else:
            xmin = xc_best - int(ceil(1.25*Coeff*R_best)) + 1	# larger by 25% of the cut-out radius * sqrt(2) <-- if the square will be rotated
            xmax = xc_best + int(ceil(1.25*Coeff*R_best))
            ymin = yc_best - int(ceil(1.25*Coeff*R_best)) + 1
            ymax = yc_best + int(ceil(1.25*Coeff*R_best))
        '''
        
        ## Measure background around the star
        backgr_level,backgr_std = backgr_around_star(scidata, xc_best, yc_best, R_best)
                
        ## Subtract this background
        sky_subtration(input_image, 'sky_subtr_psf.fits', backgr_level)
        if azim_model:
            azimProfile.main('sky_subtr_psf.fits', output_model=output_star_image, azim_tab='azim_model.txt', mask_image=None, xcen=xc_best, ycen=yc_best, ell=0., posang=0., sma_min=-1, sma_max=A_star_in, step=1., sigma_sky=None, sigma_cal=None, outside_frame=False, center_model=True, stop_where_negative=True) # TODO Sometimes it freezes!!!!!
            os.remove('azim_model.txt')

        else:    
            ## Center PSF image
            crop_image.main('sky_subtr_psf.fits', xmin, ymin, xmax, ymax, 'psf_star.fits', hdu=0)
            os.remove('sky_subtr_psf.fits')


                    
            ## Rotate the PSF image to match rotated frames (if needed)
            if PA!=0.:
                rotate_image.main('psf_star.fits', 'psf_star_rot.fits', PA, hdu_inp=0)
                shutil.move('psf_star_rot.fits', 'psf_star.fits')

            ## Center the image
            if output_star_image is None:
                output_star_image = input_image.split('/')[-1].split('.fits')[0]+'_psf.fits'
            center_psf('psf_star.fits', R_find=True, R=None, output_file=output_star_image, x_l=None, y_l=None, x_r=None, y_r=None, Std=backgr_std, verbosity=verbosity)

            ## Rebin the image if needed
            if factor!=1.:
                rebin_image.downsample(output_star_image, factor, 'psf_star_tmp.fits', set_wcs=False)
                os.rename('psf_star_tmp.fits', output_star_image)
            else:
                if verbosity: print('There is no need to rebin star psf!')
            
            os.remove('psf_star.fits')
        if norm:
            norm_image.main(output_star_image, 'tmp.fits')
            os.rename('tmp.fits', output_star_image)
        
        return output_star_image



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Crop psf star from image using region file with psf stars")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("regionFile", help="Input region file with psf stars", type=str) 
    parser.add_argument("--output_star_image", help="Optional: Output image with region", type=str, default= None)
    parser.add_argument("--posang", nargs='?', const=1, help="Optional: Position angle of the star", type=float, default=0.)
    parser.add_argument("--star_number", nargs='?', const=1, help="Optional: Number of PSF star from the region file", type=int, default=1)
    parser.add_argument("--factor", help="Optional: Sampling factor of the PSF image.", type=float, default=1.) 

    
    args = parser.parse_args()

    input_image = args.inputImage
    region_file = args.regionFile
    output_star_image = args.output_star_image
    posang = args.posang
    factor = args.factor
    star_number = args.star_number

    crop_star(input_image, region_file, output_star_image=output_star_image, star_number=star_number, factor=factor, PA=posang)