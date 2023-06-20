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
import matplotlib.image as mpimg
from numpy import *
from pylab import *
import os
import shutil

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

def main(deep_image, mask_image, blue_image, green_image, red_image, output_file, m0_deep, m0_blue, m0_green, m0_red, outer_level=25., sigma=None,text=None, x_l=None, y_l=None, x_r=None, y_r=None, same_fov=True, stretch='arcsinh',Distance=None):

    # 0. Cut out the deep image if needed:
    if sigma==None:
      sky,sigma = backEst.estimate_sky(deep_image, cleanName = None, method = 'random', degree=0, sky_subtr=False, x_obj = None, y_obj = None, coords = 'world', manual=False, reg_mask=None, ret_only_sky=False, annulus_width=None, inner_ellipse=None, Annulus_coeff=None)


    if x_l!=None and y_l!=None and x_r!=None and y_r!=None:
        imp_center.crop(deep_image, 'deep_crop.fits',x_l,y_l,x_r,y_r)
        deep_image = 'deep_crop.fits'
    
    #crea_rgbimages_ds9.crea_hist_image('deep_crop.fits')
    #crea_rgbimages_ds9.crea_hist_image('deep_crop.fits',sigma,m0_deep,0.396)
    #exit()
    #crea_rgbimages_ds9.crea_hist_image('deep_crop.fits',"ds9.png",sigma,m0_deep,0.396)
    
    # Create DS9 histogram image:
    crea_rgbimages_ds9.crea_histogram_image(deep_image,"ds9.png",sigma,cmap="grey",scale="histequ", invert=True)
    #exit()
    

    # 1. Rebin all images to deep_image
    rebin_images = []
    factors = []
    pix2secs = []
    # Determine resolution of the deep_image
    deep_res = imp_rebin.resolution(deep_image)
    deep_res = float(deep_res[0])
    
    for image in [deep_image,blue_image, green_image, red_image]:
        if image != deep_image:
            output_image = image.split('.fits')[0] + '_rebin.fits'
            if same_fov==True and x_l!=None and y_l!=None and x_r!=None and y_r!=None:
                imp_center.crop(image, output_image,x_l,y_l,x_r,y_r)
            elif same_fov==False:
                imp_rebin.rebin(deep_image, image, output_image)
            else:
               shutil.copy(image, output_image) 
            #convert_to_ZP(output_image,'tmp.fits')
            if image==blue_image:
	      m0_input = m0_blue
	    elif image==green_image:
	      m0_input = m0_green
	    elif image==red_image:
	      m0_input = m0_red
	    else:
	      m0_input = m0_deep
            convert_to_ZP(output_image, 'tmp.fits', m0_input, m0_deep)
            if m0_input!=m0_deep:
	      shutil.move('tmp.fits', output_image) 
            rebin_images.append(output_image)
        
        pix2sec = float(imp_rebin.resolution(image)[0])
        factors.append(deep_res/pix2sec)
        pix2secs.append(pix2sec)


    hdulist = pyfits.open(deep_image)
    data = hdulist[0].data
    ny,nx = np.shape(data)
    
    # fill zero values if exist:
    for image in rebin_images:
	hdulist_tmp = pyfits.open(image)
	data_tmp = hdulist_tmp[0].data
	header_tmp = hdulist_tmp[0].header
	for k in range(ny):
	  for i in range(nx):
	    if data_tmp[k,i]==0.:
	      data_tmp[k,i] = data[k,i]
	hdu = pyfits.PrimaryHDU(data_tmp, header_tmp)
	hdu.writeto(image,clobber=True)          
    


    # Plot isophotes around objects
    #import crea_iso
    #exit()
    crea_iso.main(deep_image,m0_deep,pix2secs[0],None,None,outer_level,outer_level,smooth=5,nLevels=1)
    #exit()
    
    # Create RGB image
    aplpy.make_rgb_image([rebin_images[2],rebin_images[1],rebin_images[0]], "rgb_aplpy.png",stretch_r=stretch, stretch_g=stretch, stretch_b=stretch)
    image = mpimg.imread("rgb_aplpy.png")
    ny,nx,n_rgb =  np.shape(image)

    
    # Create mask for the outer ispohotes

    data_fill = np.zeros(shape=(ny,nx))
    new_data,mask = mask_indiv.mask_array(data, data_fill, 'isophotes.reg')



    image_hist = mpimg.imread("ds9.png")
    #ny,nx,n_rgb =  np.shape(image_hist)
    #print ny,nx,n_rgb
    #exit()

    # Create overlayed image
    over_data = np.copy(image)
    for k in range(ny):
        for i in range(nx):
            if mask[ny-1-k,i]==0:
                over_data[k,i] = image_hist[k,i]

    # Apply filter to smooth the image
    img = ndimage.gaussian_filter(over_data, sigma=(2, 2, 0), order=0)


    # Plot the final image
    fig =  plt.figure(0)
    ax=plt.subplot()
    fsize = 15

    mg2 = ax.imshow(img,origin="upper", aspect='equal')

    # Overlay contours
    f = open('con1.reg','r')
    for Line in f:
          x = []
          y = []
	  if 'polygon(' in Line:
	    coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]
	    for k in range(0,len(coords)-1,2):
	      x.append(float(coords[k]))
	      y.append(ny-float(coords[k+1]))

	  x = np.array(x)
	  y = np.array(y)
	  plt.plot(x,y,'r-',color='white',lw=0.1,alpha=0.1)#markersize=1.,alpha=0.1)
    f.close()

    if text!=None:
            #plt.text(0.05, 0.03, text, fontsize=fsize, color='blue',transform=ax.transAxes)
            plt.figtext(0.5, 0.88, text,color='blue', horizontalalignment='center',verticalalignment='center')

    ax.set_frame_on(False)
    ax.set_xticks([]); ax.set_yticks([])
    plt.axis('off')
    if Distance!=None:
      kpc_per_arc = 4.84*float(Distance)/1000.
      #print kpc_per_arc
      L_bar = roundup((nx/5.)*pix2sec*kpc_per_arc, round_number=1)

      plot_image.crea_scale_bar(ax,nx,ny,pix2secs[0],l_bar=L_bar,kpc_per_arc=kpc_per_arc)
    L_bar = roundup((nx/5.)*pix2sec)
    if L_bar>=60.:
      L_bar = roundup(L_bar/60., round_number=1)
      angle_units='arcmin'
    else:
      angle_units='arcsec'
    plot_image.crea_scale_bar(ax,nx,ny,pix2secs[0],l_bar=L_bar,kpc_per_arc=None,angle_units=angle_units)
    #plot_image.axes_direct(ax,deep_image,nx,ny,pix2secs[0],reverse=True)
    plt.savefig(output_file, bbox_inches='tight', pad_inches=0)
    plt.clf()
    plt.close()

    plot_image.axes_direct(deep_image, output_file,'tmp.'+ output_file.split('.')[-1])
    shutil.move('tmp.'+ output_file.split('.')[-1],output_file)
    
    # Remove files:
    remove_func(["ds9.png",'isophotes.reg',"rgb_aplpy.png","deep_crop.fits","con1.reg"]+rebin_images)

    print 'Done!'
    '''       
    # 2. Define deepest isophotes for deep_image and others
    limits = []
    images = [deep_image,blue_image, green_image, red_image]
    m0 = [m0_deep, m0_blue, m0_green, m0_red]
    int_min = [] # intensity minima for rebinned images
    int_max = [] # intensity maxima for rebinned images
    
    for k in range(len(images)):
        limit = test_for_deepness.sky_in_boxes(images[k], mask_image, m0[k], pix2secs[k], box_size=int(math.ceil(box_size_arcsec/pix2secs[k])), Nboxes=1000,n_sigma=3)
        limit_int,max_int = test_for_deepness.sky_in_boxes(images[k], mask_image, m0[k], pix2secs[k], box_size=int(math.ceil(box_size_arcsec/pix2secs[k])), Nboxes=1000,n_sigma=3,units='DN',upper=True)
        limits.append(limit)
        
        int_min.append( limit_int ) 
        int_max.append( max_int ) 
        
    print 'Depth of the deep image is %.2f mag/arcsec^2' % (limits[0])
    print 'Depth of the blue image is %.2f mag/arcsec^2' % (limits[1])
    print 'Depth of the green image is %.2f mag/arcsec^2' % (limits[2])
    print 'Depth of the red image is %.2f mag/arcsec^2' % (limits[3])
    min_limit = max([limits[1],limits[2],limits[3]])
    '''     
    
'''
m00 = 2.5 * log10(53.907456) + 24. 
deep_image = 'f1761_rdeep.rec.fits'
mask_image = None
blue_image = 'f1761_g.rec.fits'
green_image = 'f1761_r.rec.fits'
red_image = 'f1761_i.rec.fits'
output_image = 'image.png'

main(deep_image, mask_image, blue_image, green_image, red_image, output_image, m00, m00, m00, m00, outer_level=26., text=None, x_l=2690, y_l=680, x_r=3892, y_r=1788, same_fov=True)
'''