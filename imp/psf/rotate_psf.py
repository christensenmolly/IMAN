#!/usr/bin/python
# -*- coding:  cp1251 -*-
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

from astropy.io import fits as pyfits
from astropy.modeling import models, fitting
import itertools
from scipy import ndimage

tmp_out = sys.stdout


LOCAL_DIR = "/imp/psf"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/rotate'))

#from Sergey_pipelinelib.rotima import crotima
import rotate_image
from astropy.stats import sigma_clipped_stats
from photutils import centroid_com, centroid_1dg, centroid_2dg, centroids

# Colors to highlight the output text
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

def find_nearest(array,value):
    arr = np.abs(array-value)
    cen = np.where(arr==np.min(arr))
    xcc = cen[1][0]
    ycc = cen[0][0]
    return xcc,ycc



def rotate_psf(input_file, xc, yc, PosAng, find_center='max'):
  #find_center = ['max','com','old']
  if PosAng!=0.:
            hdulist = pyfits.open(input_file)
            data = hdulist[0].data
            ny,nx = np.shape(data)
            
            #xCenRot, yCenRot = crotima(input_file, 'rot_psf.fits',xc,yc,PosAng, set_wcs = False)
            xCenRot, yCenRot = rotate_image.main(input_file, PosAng, xc=xc, yc=yc, output_image='rot_psf.fits', hdu_inp=0, cval=float('nan'), cropping=False, verbosity=True)
            xCenRot = int(xCenRot) - 1
            yCenRot = int(yCenRot) - 1

            hdulist1 = pyfits.open('rot_psf.fits')
            data1 = hdulist1[0].data

            R = int(min([xc,yc,nx-xc,ny-yc]))

            if find_center=='max':
                cen =  np.where(data1==np.max(data1))
                xcc = cen[1][0]
                ycc = cen[0][0]
                #print xcc,ycc 
                #exit()
            elif find_center=='com':
                xcc, ycc = centroid_com(data1)
                xcc = int(floor(xcc))
                ycc = int(floor(ycc))
                #print xcc,ycc
                #exit()
            elif find_center=='gauss':
                table = centroids.fit_2dgaussian(data1)
                #print arr['x_mean']
                #print arr['y_mean']
                xcc = int(floor(table.x_mean[0]))
                ycc = int(floor(table.y_mean[0]))
                #exit()
                #print arr
            else:
                norm_data = data/np.sum(data)
                if nx%2==0:
                    xc_psf = int(nx/2. + 1)
                else:
                    xc_psf = int(nx/2. + 0.5)
                if ny%2==0:
                    yc_psf = int(ny/2. + 1)
                else:
                    yc_psf = int(ny/2. + 0.5)
                #print 'here',yc_psf,xc_psf
                I_norm = norm_data[yc_psf,xc_psf]
                
                norm_rot_data = data1/np.sum(data1)
                
                xcc,ycc = find_nearest(norm_rot_data,I_norm)
                #print xcc,ycc
                #exit()
        
            new_data = data1[ycc-R:ycc+R+1,xcc-R:xcc+R+1]
            
            hdu = pyfits.PrimaryHDU(new_data/np.sum(new_data))
            hdu.writeto('rot_psf.fits',clobber=True)
            #exit()
  else:
            shutil.copy(input_file,'rot_psf.fits') 

def center_psf(input_image, output_image, xc, yc, R=None):
        hdulist = pyfits.open(input_image)
        data = hdulist[0].data
        ny,nx = np.shape(data)
        if R==None:
            R = int(min([xc,yc,nx-xc,ny-yc]))
        
        new_data = data[yc-R:yc+R+1,xc-R:xc+R+1]
        
        hdu = pyfits.PrimaryHDU(new_data/np.sum(new_data))
        hdu.writeto(output_image,clobber=True)

