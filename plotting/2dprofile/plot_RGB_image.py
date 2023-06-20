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
import glob

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
sys.path.append(PATH_TO_PACKAGE+'/Warps')
sys.path.append(PATH_TO_PACKAGE+'/Plot')
sys.path.append(PATH_TO_PACKAGE+'/IMP_NEW')
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

def main(input_images, output_file, reproject = True, text=None, x_l=None, y_l=None, x_r=None, y_r=None, Distance=None, pix2sec=None, hor_pos=0.03, vert_pos=0.8):
    
    if x_l!=None and y_l!=None and x_r!=None and y_r!=None:
        Input_images = []
        for input_image in input_images:
            imp_center.crop(input_image, input_image.split('.fits')[0]+'_crop.fits',x_l,y_l,x_r,y_r)
            Input_images.append(input_image.split('.fits')[0]+'_crop.fits')
    else:
        Input_images = input_images

    if reproject==True:
        aplpy.make_rgb_cube([Input_images[0],Input_images[1],Input_images[2]], 'cube_tmp.fits')
        aplpy.make_rgb_image('cube_tmp.fits', 'tmp.png')
        os.remove('cube_tmp.fits')
    else:
        aplpy.make_rgb_image([Input_images[0],Input_images[1],Input_images[2]], 'tmp.png')

    # Get image resolution
    if pix2sec == None:
        pix2sec,comments = imp_rebin.resolution(input_images[0])
        pix2sec = float(pix2sec)



    hdulist = pyfits.open(Input_images[0])
    data = hdulist[0].data
    ny,nx = np.shape(data)

    # Plot final image
    fig =  plt.figure(0)
    ax = plt.subplot()
    fsize = 13
    
    #min_int = 10**(0.4*(m0-outer_level)) * pix2sec
    #vmax = np.max(data)
    #color = 'terrain_r'#'nipy_spectral_r'#
    #mg2 = ax.imshow(data, cmap=color, norm=LogNorm(vmin=min_int, vmax=vmax),interpolation='none', origin="lower", aspect='equal')

    img = mpimg.imread('tmp.png')
    mg2 = ax.imshow(img)

    
    if text!=None:
            ax.text(hor_pos, vert_pos, text, fontsize=fsize, color='white',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline')
            #plt.figtext(0.5, 0.88, text,color='blue', horizontalalignment='center',verticalalignment='center')
    length_arc = 30.
    length = length_arc/pix2sec
    yc = ny-ny/7.
    xc = nx - length/2.-nx/20.
    

    ax.set_frame_on(False)
    ax.set_xticks([]); ax.set_yticks([])
    plt.axis('off')
    
    if Distance!=None:
      kpc_per_arc = 4.84*float(Distance)/1000.
      L_bar = roundup((nx/5.)*pix2sec*kpc_per_arc, round_number=1)

      plot_image.crea_scale_bar(ax,nx,ny,pix2secs[0],l_bar=L_bar,kpc_per_arc=kpc_per_arc)
    L_bar = roundup((nx/5.)*pix2sec)
    if L_bar>=60.:
      L_bar = roundup(L_bar/60., round_number=1)
      angle_units='arcmin'
    else:
      angle_units='arcsec'
    
    #plot_image.crea_scale_bar(ax,nx,ny,pix2sec,l_bar=L_bar,kpc_per_arc=None,angle_units=angle_units)
    ax.errorbar(xc, yc, xerr=length/2., color='lime', capsize=2, c='red')
    plt.savefig(output_file, bbox_inches='tight', pad_inches=0)
    plt.clf()
    plt.close()

    #plot_image.axes_direct(input_image, output_file,'tmp.'+ output_file.split('.')[-1])

    
    # Remove tmp files:
    files = glob.glob('*_crop.fits')
    files.append('tmp.png')
    for file in files:
        os.remove(file)

    print 'Done!'
 