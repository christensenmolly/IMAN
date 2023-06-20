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
import collections
import warnings
warnings.filterwarnings("ignore")
tmp_out = sys.stdout
FNULL = open(os.devnull, 'w')


PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0].split('IMP_NEW')[0]
sys.path.append(PATH_TO_PACKAGE + '/DECA/NEW')

import rotate_psf
import galfit_input


def show_complete(i,N):
  percent = 100. * float(i) / float(N)
  sys.stdout.write("\r%2d%%" % percent)
  sys.stdout.flush()


def estim_PA(segm_image, xc, yc):
      from skimage.morphology import skeletonize
      hdulist = pyfits.open(segm_image)
      data = hdulist[0].data
      nx, ny =data.shape[1], data.shape[0]
      
      star_int = data[int(yc),int(xc)]
      ind =  np.where(data==star_int)
      star_data = data[min(ind[0]):max(ind[0]),min(ind[1]):max(ind[1])]
      nx, ny =star_data.shape[1], star_data.shape[0]
      
      final = np.zeros(shape=(ny,nx))
      np.putmask(final, star_data==star_int, 1.)
      
      skeleton = skeletonize(final)
      #print skeleton
      #plt.imshow(skeleton, cmap=plt.cm.gray)
      #plt.show()


  

def main(file_Image, input_psf_file, xc_psf, yc_psf, xc, yc, sky_level, xmin, ymin, xmax, ymax):
      PA = 0.
      mag = 10.
      file_Sigma = 'none'
      file_Mask = 'none'
      file_constraints = 'none'
      file_galfit_outimage = 'model.fits'
      magZP = 18.536
      generalScaleValue = 1.375
      sampling = 1.
      Chi2 = []
      PAS = []
      
      while PA<360.:
            #print PA
            # Rotate PSF
            if PA!=0.:
                rotate_psf.rotate_psf(input_psf_file, xc_psf, yc_psf, PA)
            else:
                shutil.copy(input_psf_file, 'rot_psf.fits')
            
            # Create input galfit file:
            f = open('galfit.inp', "w") 
            sys.stdout = f  
            galfit_input.header(file_Image,file_Sigma,'rot_psf.fits','mask_psf.fits',file_constraints,file_galfit_outimage,xmin,xmax,ymin,ymax,magZP,generalScaleValue,sampling,write_to=None)
            galfit_input.Agn(1,xc,yc,mag)
            galfit_input.Sky(2,sky_level)
            sys.stdout = tmp_out
            f.close()
            os.chmod(r"galfit.inp",777) 
            subprocess.call("galfit %s" % ("galfit.inp"), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
            try:
                ff = open('galfit.01', 'r')
                for line in ff:
                    if '#  Chi^2/nu = ' in line:
                        Chi2.append(float(line.split('#  Chi^2/nu = ')[-1].split(',')[0]))
                        #print PA,'\t',float(line.split('#  Chi^2/nu = ')[-1].split(',')[0])
                    if '#  Integrated magnitude' in line:
                        mag = float(line.split(' 3) ')[-1].split()[0])
                        #print mag
                    if '#  Position x, y' in line:
                        xc = float(line.split(' 1) ')[-1].split()[0])
                        yc = float(line.split(' 1) ')[-1].split()[1])
                        #print xc,yc
                ff.close()
                PAS.append(PA)
                os.remove("galfit.01")
            except:
                z=1
            PA = PA + 1.
            show_complete(abs(PA),360.)
	
      #print PAS
      #print Chi2
      #print PAS
      #print Chi2
      #plt.plot(PAS,Chi2,'o')
      #plt.show()
      best_PA = PAS[Chi2.index(min(Chi2))]
      #print best_PA
      os.remove('fit.log')
      os.remove('model.fits')
      os.remove('rot_psf.fits')
      os.remove('mask_psf.fits')
      os.remove('galfit.inp')
      return best_PA


def main_fast(file_Image, input_psf_file, xc_psf, yc_psf, xc, yc, sky_level, xmin, ymin, xmax, ymax):
      mag = 10.
      file_Sigma = 'none'
      file_Mask = 'none'
      file_constraints = 'none'
      file_galfit_outimage = 'model.fits'
      magZP = 18.536
      generalScaleValue = 1.375
      sampling = 1.
      Chi2 = []
      PAS = []

      PA = -2.      
      # 1. Find best PA in [-2:92] 
      while PA<=92:
            # Rotate PSF
            #print PA
            if PA!=0.:
                rotate_psf.rotate_psf(input_psf_file, xc_psf, yc_psf, PA)
            else:
                shutil.copy(input_psf_file, 'rot_psf.fits')
            
            # Create input galfit file:
            f = open('galfit.inp', "w") 
            sys.stdout = f  
            galfit_input.header(file_Image,file_Sigma,'rot_psf.fits','mask_psf.fits',file_constraints,file_galfit_outimage,xmin,xmax,ymin,ymax,magZP,generalScaleValue,sampling,write_to=None)
            galfit_input.Agn(1,xc,yc,mag)
            galfit_input.Sky(2,sky_level)
            sys.stdout = tmp_out
            f.close()
            os.chmod(r"galfit.inp",777) 
            subprocess.call("galfit %s" % ("galfit.inp"), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
            try:
                ff = open('galfit.01', 'r')
                for line in ff:
                    if '#  Chi^2/nu = ' in line:
                        Chi2.append(float(line.split('#  Chi^2/nu = ')[-1].split(',')[0]))
                    if '#  Integrated magnitude' in line:
                        mag = float(line.split(' 3) ')[-1].split()[0])
                    if '#  Position x, y' in line:
                        xc = float(line.split(' 1) ')[-1].split()[0])
                        yc = float(line.split(' 1) ')[-1].split()[1])
                ff.close()
                PAS.append(PA)
                os.remove("galfit.01")
            except:
                z=1
            PA = PA + 1.
            #show_complete(abs(PA),94.)
      
      # Find best PA:
      best_PA = PAS[Chi2.index(min(Chi2))]
      #print best_PA
      #plt.xlabel('PA (deg)',fontsize=15)
      #plt.ylabel('Chi2', fontsize=15)      
      #plt.plot(PAS,Chi2,'o')
      #plt.show()


      #print '\n'
      # Consider the three other quadrants:
      for PA in [best_PA+90.-2.,best_PA+180.-2.,best_PA+270.-2.]:
            #print PA
            PA_max = PA + 2.
            while PA<=PA_max:
                # Rotate PSF
                if PA!=0.:
                    rotate_psf.rotate_psf(input_psf_file, xc_psf, yc_psf, PA)
                else:
                    shutil.copy(input_psf_file, 'rot_psf.fits')
                
                # Create input galfit file:
                f = open('galfit.inp', "w") 
                sys.stdout = f  
                galfit_input.header(file_Image,file_Sigma,'rot_psf.fits','mask_psf.fits',file_constraints,file_galfit_outimage,xmin,xmax,ymin,ymax,magZP,generalScaleValue,sampling,write_to=None)
                galfit_input.Agn(1,xc,yc,mag)
                galfit_input.Sky(2,sky_level)
                sys.stdout = tmp_out
                f.close()
                os.chmod(r"galfit.inp",777) 
                subprocess.call("galfit %s" % ("galfit.inp"), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                try:
                    ff = open('galfit.01', 'r')
                    for line in ff:
                            if '#  Chi^2/nu = ' in line:
                                Chi2.append(float(line.split('#  Chi^2/nu = ')[-1].split(',')[0]))
                            if '#  Integrated magnitude' in line:
                                mag = float(line.split(' 3) ')[-1].split()[0])
                            if '#  Position x, y' in line:
                                xc = float(line.split(' 1) ')[-1].split()[0])
                                yc = float(line.split(' 1) ')[-1].split()[1])
                    ff.close()
                    PAS.append(PA)
                    os.remove("galfit.01")
                except:
                    z=1
                PA = PA + 1.

      
      #print PAS
      #print Chi2
      #print PAS
      #print Chi2
      best_PA = PAS[Chi2.index(min(Chi2))]
      #print best_PA
      #plt.plot(PAS,Chi2,'o')
      #plt.show()

      os.remove('fit.log')
      os.remove('model.fits')
      os.remove('rot_psf.fits')
      os.remove('mask_psf.fits')
      return best_PA


'''
file_Image = '3_ESO149-001_WISE_3.4.fits'
input_psf_file = 'PSF_W1.V4_rebin.fits'
xc_psf = 120.
yc_psf = 120.
xc = 813.1076
yc = 564.2584
sky_level = 2.1378
xmin = 715
xmax = 911
ymin = 473
ymax = 642

main(file_Image, input_psf_file, xc_psf, yc_psf, xc, yc, sky_level, xmin, ymin, xmax, ymax)
'''
