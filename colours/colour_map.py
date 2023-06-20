#!/usr/bin/python
# DESCRIPTION:
# Script to create a colour map as a residual between two images of the same object in different wavebands. It takes into account PSF for each image (to reduce to the same resolution) and Galactic extinction in each band. Note that both images should have the same scale and orientation.
# MINIMAL USAGE: python colour_map.py [input_image_1] [input_image_2] [psf_image_1] [psf_image_2] [pixel_scale] [A1] [A2] [ZP1] [ZP2] 


import os
import sys
from pylab import *
import astropy.io.fits as pyfits
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import glob
from matplotlib.pyplot import cm
import re
from scipy.odr.odrpack import *
import pyparsing
import matplotlib.pyplot as pyplot
import matplotlib as mpl
import numpy as np
import subprocess
import argparse
import shutil
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.ndimage as ndimage

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
sys.path.append(PATH_TO_PACKAGE+'/misc_funcs')
sys.path.append(PATH_TO_PACKAGE+'/plotting/2dprofile')
import convolve
import plot_2d_profile



def discrete_cmap(N=8):
    # define individual colors as hex values
    cpool = [ '#000000', '#00EE00', '#0000EE', '#00EEEE', '#EE0000','#FFFF00', '#EE00EE', '#FFFFFF']
    cmap_i8 = mpl.colors.ListedColormap(cpool[0:N], 'i8')
    mpl.cm.register_cmap(cmap=cmap_i8)
    return cmap_i8

def color_ranges(color):
  if color == 'g-i':
    return 0,1.5  
  if color == 'g-r':
    return -0.5,1.5  
  if color == 'u-g':
    return -0.5,2.5 
  if color == 'u-r':
    return 1.0,3.0 
  else:
    return -3.0,3.0 

def line_reg(header1,ima_pix2sec):
  nx = int(header1['NAXIS1'])
  ny = int(header1['NAXIS2'])
  scale = int(round(nx/8.*ima_pix2sec,-1))
  x2 = nx*9.8/10.
  x1 = x2 - scale/ima_pix2sec
  y1 = ny/7.
  y2 = y1
  
  return x1,y1,x2,y2,scale

def main(color,input_image_1,input_image_2,psf_image_1=None,psf_image_2=None,pix2sec=1.,A1=0.,A2=0.,m01=0.,m02=0.,zoom=None,rms1=0.,rms2=0.,text=None,scale_bar=None,show_bar_label=False, sigma_smooth=5., output_image='colour_map.fits', SB_min1=26.):
  if psf_image_1 is not None and psf_image_2 is not None:
    # convolve images:
    convolve.convolution(input_image_1, psf_image_2, input_image_1.split('.fits')[0]+'_conv.fits')
    convolve.convolution(input_image_2, psf_image_1, input_image_2.split('.fits')[0]+'_conv.fits')
  else:
      shutil.copy(input_image_1, input_image_1.split('.fits')[0]+'_conv.fits')
      shutil.copy(input_image_2, input_image_2.split('.fits')[0]+'_conv.fits')
      
  input_image_1 = input_image_1.split('.fits')[0]+'_conv.fits'
  input_image_2 = input_image_2.split('.fits')[0]+'_conv.fits'
  
  color_min,color_max = color_ranges(color)


  hdu1 = pyfits.open(input_image_1)
  data1 = hdu1[0].data
  header1 = hdu1[0].header
  ySize, xSize = data1.shape
  xc = xSize/2.
  yc = ySize/2.
  try:
    m01 = header['M0']
  except:
    zz =1
  hdu2 = pyfits.open(input_image_2)
  data2 = hdu2[0].data
  header2 = hdu2[0].header
  try:
    m02 = header['M0']
  except:
    zz =1
  
  data1 = ndimage.gaussian_filter(data1, sigma=sigma_smooth, order=0)
  data2 = ndimage.gaussian_filter(data2, sigma=sigma_smooth, order=0)
  
  div12 = data1/data2#ndimage.gaussian_filter(data1/data2, sigma=5., order=0)
  
  data_col = m01 - 2.5*np.log10(div12) - A1 - m02 + A2
  


  #for k in range(ySize):
  #  for i in range(xSize):
  #    if data1[k,i]>1.*rms1 and data2[k,i]>1.*rms2:
  #               mag1 = m01 - 2.5*log10(data1[k,i]) - A1
  #               mag2 = m02 - 2.5*log10(data2[k,i]) - A2
  #               data_col[k,i]= mag1 - mag2
  #    else:
  #               data_col[k,i]=float(nan)
                 
  if zoom is not None:
      yl = int(yc-yc/zoom)
      yr = int(yc+yc/zoom)
      xl = int(xc-xc/zoom)
      xr = int(xc+xc/zoom)
      data_col = data_col[yl:yr,xl:xr]
  
  ySize, xSize = data_col.shape
  
  #data_col = ndimage.gaussian_filter(data_col, sigma=3., order=0)
  
  data_col[data1<(10**(0.4*(m01-SB_min1)))*(pix2sec**2)]=np.nan

  hdu = pyfits.PrimaryHDU(data_col)
  hdu.writeto(output_image, overwrite=True)
  return 0

  
  exit()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  fig = plt.figure(6,figsize=(5, 5))
  ax = plt.gca()
  plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
  #im = ax.imshow(data_col,cmap=discrete_cmap(), vmin=color_min, vmax=color_max,interpolation='none', origin="lower", aspect='equal')
  im = ax.imshow(data_col,cmap='rainbow', vmin=color_min, vmax=color_max,interpolation='none', origin="lower", aspect='equal')
  
  
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)

  ax.set_axis_off()
  ax.set_xlim(0,xSize)
  ax.set_ylim(0,ySize)
  cbar = plt.colorbar(im, cax=cax)
  cbar.ax.tick_params(labelsize=10)
  x0 = 0
  y0 = 0
  x1 = xSize
  y1 = ySize  

  if scale_bar!=None:     
    plot_2d_profile.crea_scale_bar(ax,x0,x1,y0,y1,pix2sec,float(scale_bar),show_bar_label)

  if text!=None:
            plt.figtext(0.5, 0.88, text,color='black', size='x-large', horizontalalignment='center',verticalalignment='top')
  plt.draw()
  plt.savefig('color.eps', transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.)
  plt.clf()
  plt.close()   




  # Color along minor axis
  plt.figure(4,figsize=(5, 5))



  plt.xlabel(r' z (pix) ', fontsize=15)
  plt.ylabel(r' %s ' % (color), fontsize=15)
  
  z = np.arange(0,ySize,1) - ySize/2.
  plt.plot(z,data_col[:,int(xSize/2.)],'*',color='white')
  
  

  I = []
  z = []
  for y in range(0,ySize,1):
                 #I.append(sum(data[y,:]))
                 II = []
                 for x in range(0,xSize,1):
                   if data_col[y,x]>0.:
                     II.append(data_col[y,x])
                 I.append(sum(II))
                 z.append(y-ySize/2.)
  z = np.array(z)
  I = np.array(I)
  
  plt.plot(z,I*100./sum(I),'o',color='black')
  ylim(0,2)
  
  plt.savefig('color_vertical.eps', transparent = False, dpi=300)
  plt.clf()
  plt.close()     
  

  # Color along major axis
  plt.figure(5,figsize=(5, 5))



  plt.xlabel(r' r (pix) ', fontsize=15)
  plt.ylabel(r' %s ' % (color), fontsize=15)
  
  r = np.arange(0,xSize,1) - xSize/2.
  #plt.plot(r,data_col[int(ySize/2.),:],'*',color='white')
  
  

  I = []
  r = []
  for x in range(0,xSize,1):
                 #I.append(sum(data[y,:]))
                 II = []
                 for y in range(0,ySize,1):
                   if data_col[y,x]>0.:
                     II.append(data_col[y,x])
                 I.append(sum(II))
                 r.append(x-xSize/2.)
  r = np.array(r)
  I = np.array(I)
  
  plt.plot(r,I*100./sum(I),'o',color='black')
  ylim(0,2)
  
  plt.savefig('color_horizontal.eps', transparent = False, dpi=300)
  plt.clf()
  plt.close()
  os.remove(input_image_1)
  os.remove(input_image_2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="1D and 2D colour maps")
    parser.add_argument("Color", help="Input colour, e.g. 'g-i")
    parser.add_argument("input_image_1", help="Input first image")
    parser.add_argument("input_image_2", help="Input second image")
    parser.add_argument("--psf_image_1", help="Input first psf image",type=str,default=None)
    parser.add_argument("--psf_image_2", help="Input second psf image",type=str,default=None)
    parser.add_argument("--Scale", help="Input scale in [arcsec/pix]",type=float,default=1.)
    parser.add_argument("--A1", nargs='?', const=0., help="Input Galactic extinction in first band",type=float,default=0.)
    parser.add_argument("--A2", nargs='?', const=0., help="Input Galactic extinction in second band",type=float,default=0.)
    parser.add_argument("--m01", nargs='?', const=0., help="Input Zero Point in first band in [mag/arcsec^2]",type=float,default=0.)
    parser.add_argument("--m02", nargs='?', const=0., help="Input Zero Point in second band in [mag/arcsec^2]",type=float,default=0.)     
    parser.add_argument("--scale_bar", nargs='?', const=0., help="Create the scale bar: length [arcsec]",type=str,default=None)
    parser.add_argument("--text", nargs='?', const=0., help="Create title of the image",type=str,default=None)
    parser.add_argument("--show_bar_label", nargs='?', const=0., help="Show the scale bar label: yes",type=str,default='yes')

    parser.add_argument("--rms1", nargs='?', const=0., help="rms1",type=float,default=0.)
    parser.add_argument("--rms2", nargs='?', const=0., help="rms2",type=float,default=0.)
    parser.add_argument("--o", help="Output_image",type=str,default='colour_map.fits')
    args = parser.parse_args()
    
    color = args.Color
    input_image_1 = args.input_image_1
    input_image_2 = args.input_image_2
    psf_image_1 = args.psf_image_1
    psf_image_2 = args.psf_image_2    
    pix2sec = args.Scale
    A1 = args.A1
    A2 = args.A2
    m01 = args.m01
    m02 = args.m02
    scale_bar = args.scale_bar
    show_bar_label = args.show_bar_label    
    text = args.text
    
    rms1 = args.rms1
    rms2 = args.rms2
    output_image = args.o
    
    main(color,input_image_1,input_image_2,psf_image_1,psf_image_2,pix2sec,A1=A1,A2=A2,m01=m01,m02=m02,text=text,scale_bar=scale_bar,show_bar_label=show_bar_label, rms1=rms1, rms2=rms2, output_image=output_image, sigma_smooth=5., SB_min1=28.)
