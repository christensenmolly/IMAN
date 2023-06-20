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
import pyfits
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

def main(input_psf_image, output_psf_image, Radius_in, Radius_out):
    # Read in the input PSF 
    hdulist = pyfits.open(input_psf_image)
    data = hdulist[0].data
    nx, ny =data.shape[1], data.shape[0]
    
    # Determine the center:
    if nx%2==False:
      xc = nx/2. + 0.5
    else:
      xc = nx/2. + 1.

    if ny%2==False:
      yc = ny/2. + 0.5
    else:
      yc = ny/2. + 1.
    
    m0 = 20.
    # Create the mask: fit only the inner part of the kernel within Radius_in
    data_mask = np.zeros(shape=(ny,nx))
    y,x = np.ogrid[-int(yc):data.shape[0]-int(yc), -int(xc):data.shape[1]-int(xc)]
    mask = x*x + y*y <= Radius_in*Radius_in
    data_mask[mask] = 1

    y,x = np.ogrid[-int(yc):data.shape[0]-int(yc), -int(xc):data.shape[1]-int(xc)]
    mask = x*x + y*y > max([xc,yc])*max([xc,yc])
    data_mask[mask] = 1
    
    for k in range(ny):
      for i in range(nx):
	if data[k,i]<0.:
	  data_mask[k,i]=1
    
    hdu = pyfits.PrimaryHDU(data_mask)
    hdu.writeto('mask_psf.fits',clobber=True)
    #exit()
    # Create input galfit file:
    f = open('galfit.inp', "w") 
    sys.stdout = f  
    galfit_input.header(input_psf_image,'none','none','mask_psf.fits','none','model.fits',1,nx,1,ny,10.,1.,1.,write_to=None)
    galfit_input.Sersic(component=1,xc=xc,yc=yc,meb=float('nan'),
		      Mb = m0 - 2.5*log10(np.sum(data)),reb=nx/10.,
		      n=2.,q=1.,PA=0.)
    sys.stdout = tmp_out
    f.close()
    os.chmod(r"galfit.inp",0777) 
    subprocess.call("galfit %s" % ("galfit.inp"), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    
    # Read in the ouput galfit file:
    f = open('galfit.01', 'r')
    for line in f:
      if '#  Integrated magnitude' in line:
	Mb = float(line.split()[1])
      if '#  R_e (effective radius)   [pix]' in line:
	reb = float(line.split()[1])
      if '#  Sersic index n (de Vaucouleurs n=4)' in line:
	n = float(line.split()[1])    
      if '#  Axis ratio (b/a)' in line:
	q = float(line.split()[1])
      if '#  Position angle (PA) [deg: Up=0, Left=90]' in line:
	PA = float(line.split()[1])
    f.close()
    # Create the model of larger size:
    nx_new = int(xc + Radius_out)
    ny_new = int(yc + Radius_out)
    if nx_new%2:
      nx_new = nx_new + 1

    if ny_new%2:
      ny_new = ny_new + 1

    if nx_new%2==False:
      xc_new = nx_new/2. + 0.5
    else:
      xc_new = nx_new/2. + 1.

    if ny_new%2==False:
      yc_new = ny_new/2. + 0.5
    else:
      yc_new = ny_new/2. + 1.
    
    
    f = open('galfit.inp', "w") 
    sys.stdout = f  
    galfit_input.header('none','none','none','none','none','model.fits',1,ny_new,1,ny_new,10.,1.,1.,write_to=None)
    galfit_input.Sersic(component=1,xc=xc_new,yc=yc_new,meb=float('nan'),
		      Mb=Mb, reb=reb,
		      n=n, q=q, PA=PA)
    sys.stdout = tmp_out
    f.close()
    os.chmod(r"galfit.inp",0777) 
    subprocess.call("galfit -o3 %s" % ("galfit.inp"), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    
    # Combine model and real PSF:
    hdulist_model = pyfits.open('model.fits')
    model = hdulist_model[0].data
    nx_mod, ny_mod =model.shape[1], model.shape[0]
    
    final_data = np.copy(model)
    
    x_l = int(xc_new)-int(math.floor(xc))
    x_r = x_l + nx
    y_l = int(yc_new)-int(math.floor(yc))
    y_r = y_l + ny
    
    final_data[y_l:y_r, x_l:x_r] = data
    hdu = pyfits.PrimaryHDU(final_data)
    hdu.writeto(output_psf_image, clobber=True)    
    print 'Done!'
   
main('PSF_W1.V4_rebin.fits', 'psf.fits', 70, 170)    
    
  