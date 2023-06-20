#!/usr/bin/python
# -*- coding:  cp1251 -*-
import sys
import math
import numpy as np
import re
from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LogNorm
import matplotlib.image as mpimg
from numpy import *
from pylab import *
import os
import shutil

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
sys.path.append(PATH_TO_PACKAGE+'/Warps')
sys.path.append(PATH_TO_PACKAGE+'/Plot')
import crea_iso


import imp_rebin
import test_for_deepness
import crea_rgbimages
import imp_center
import aplpy
import crea_rgbimages_ds9
import pyfits
import mask_indiv

box_size_arcsec = 10.


import plot_image
import scipy.ndimage as ndimage
import backEst

def remove_func(files):
    for file in files:
        if os.path.exists(file):
            os.remove(file)

def roundup(x, round_number=10):
    return int(math.ceil(x / float(round_number))) * round_number

def convert_to_ZP(input_image, output_image, m0_input, m0_reference):
  if m0_input!=m0_reference:
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    header = hdulist[0].header
    data = data*10**(0.4*(m0_reference-m0_input))
    
    hdu = pyfits.PrimaryHDU(data,header)
    hdu.writeto(output_image,clobber=True)    
    #print 'here'

def main(input_image, output_file,m0,pix2sec, outer_level=25., sigma=None,text=None, x_l=None, y_l=None, x_r=None, y_r=None, same_fov=True, stretch='arcsinh',Distance=None,plot_scale_size=True,L_bar=None):

    # 0. Cut out the deep image if needed:
    '''
    if sigma==None:
      sky,sigma = backEst.estimate_sky(input_image, cleanName = None, method = 'random', degree=0, sky_subtr=False, x_obj = None, y_obj = None, coords = 'world', manual=False, reg_mask=None, ret_only_sky=False, annulus_width=None, inner_ellipse=None, Annulus_coeff=None)
    '''

    if x_l!=None and y_l!=None and x_r!=None and y_r!=None:
        imp_center.crop(input_image, 'cropp.fits',x_l,y_l,x_r,y_r)
        input_image = 'cropp.fits'


    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    ny,nx = np.shape(data)

    # Plot the final image
    fig =  plt.figure(0)
    ax=plt.subplot()
    fsize = 15
    
    min_int = 10**(0.4*(m0-outer_level)) * pix2sec
    vmax = np.max(data)
    color = 'terrain_r'#'nipy_spectral_r'#
    mg2 = ax.imshow(data, cmap=color, norm=LogNorm(vmin=min_int, vmax=vmax),interpolation='none', origin="lower", aspect='equal')



    if text!=None:
            #plt.text(0.05, 0.03, text, fontsize=fsize, color='blue',transform=ax.transAxes)
            plt.figtext(0.5, 0.88, text,color='blue', size='xx-large', horizontalalignment='center',verticalalignment='top')

    ax.set_frame_on(False)
    ax.set_xticks([]); ax.set_yticks([])
    plt.axis('off')
    if Distance!=None:
      kpc_per_arc = 4.84*float(Distance)/1000.
      #print kpc_per_arc
      if L_bar==None:
        L_bar = roundup((nx/5.)*pix2sec*kpc_per_arc, round_number=1)

      plot_image.crea_scale_bar(ax,nx,ny,pix2secs[0],l_bar=L_bar,kpc_per_arc=kpc_per_arc, plot_scale_size=plot_scale_size)
    if L_bar==None:
        L_bar = roundup((nx/5.)*pix2sec)
    if L_bar>=60.:
      L_bar = roundup(L_bar/60., round_number=1)
      angle_units='arcmin'
    else:
      angle_units='arcsec'
    
    #cmap = matplotlib.cm.jet
    #cmap.set_bad('w',1.)
    
    plot_image.crea_scale_bar(ax,nx,ny,pix2sec,l_bar=L_bar,kpc_per_arc=None,angle_units=angle_units, plot_scale_size=plot_scale_size)
    plt.savefig(output_file, bbox_inches='tight', pad_inches=0)
    plt.clf()
    plt.close()

    plot_image.axes_direct(input_image, output_file,'tmp.'+ output_file.split('.')[-1])

    
    # Remove files:
    remove_func(['cropp.fits'])

    print 'Done!'
 