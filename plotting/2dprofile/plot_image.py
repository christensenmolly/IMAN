#!/usr/bin/python


# EXAMPLE: python 
#from __future__ import unicode_literals
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
from mpl_toolkits.axes_grid1 import AxesGrid
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
import re
import glob
from matplotlib.colors import LogNorm
from astropy.io import fits as pyfits
import argparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import wcs
import warnings
warnings.filterwarnings("ignore")

import numpy.linalg as la
import PIL
from PIL import Image
from astropy import wcs



matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

fsize = 25
x_off = 0.95
y_off = 0.05


 
#@xl_func("numpy_row v1, numpy_row v2: float")
def py_ang(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'    """
    cosang = np.dot(v1, v2)
    sinang = la.norm(np.cross(v1, v2))
    return degrees(np.arctan2(sinang, cosang))

def crea_scale_bar(ax,nx,ny,pix2sec,l_bar=30.,kpc_per_arc=None,angle_units='arcsec', plot_scale_size=True):
      if np.isnan(kpc_per_arc):
          kpc_per_arc = None
      
      if angle_units=='arcmin':
        # l_bar in arcses
        L_bar = l_bar * 60.
        
        xc = nx*x_off - L_bar/(2.*pix2sec)
        yc = ny*y_off

        ax.errorbar(xc, yc, xerr=L_bar/(2.*pix2sec),color='lime',capsize=3,elinewidth=2)
        if plot_scale_size==True and kpc_per_arc is None:
            ax.text(xc, yc, str(int(l_bar))+'\'', color='lime',fontsize=8, horizontalalignment='center', verticalalignment='bottom')
      else:
        # l_bar in arcses
        L_bar = l_bar
        
        xc = nx*x_off - L_bar/(2.*pix2sec)
        yc = ny*y_off
        
        ax.errorbar(xc, yc, xerr=l_bar/(2.*pix2sec),color='lime',capsize=2,elinewidth=2)
        if plot_scale_size==True and kpc_per_arc is None:
            ax.text(xc, yc, str(int(l_bar))+'\"', color='lime',fontsize=8, horizontalalignment='center', verticalalignment='bottom')
      
      if kpc_per_arc is not None:
            ax.text(xc, yc, str(round(L_bar*kpc_per_arc,1))+'kpc', color='lime',fontsize=12, horizontalalignment='center', verticalalignment='bottom')          
          
          

      #else:
      #  yc = ny*(1.-y_off)
      #  # l_bar in kpc
      #  L_bar = l_bar / kpc_per_arc        # Now in arcses
      #  xc = nx*x_off - L_bar/(2.*pix2sec)
      #  ax.errorbar(xc, yc, xerr=L_bar/(2.*pix2sec),color='lime',capsize=2,c='lime',elinewidth=2)
      #  if plot_scale_size==True:
      #    ax.text(xc, yc, str(int(l_bar))+'kpc', color='lime',fontsize=8, horizontalalignment='center', verticalalignment='bottom')
    
    
'''
def axes_direct(ax,input_image,nx,ny,pix2sec,reverse=False):
    length_of_arrow = min([nx,ny]) / 10.
    length_of_arrow_arc = length_of_arrow * pix2sec
    
    x0 = nx*(1.-x_off) + length_of_arrow
    if reverse==False:    
        y0 = ny*y_off + length_of_arrow
    else:
        y0 = ny*(1.-y_off) - length_of_arrow        
  
    referenceheader = pyfits.getheader(input_image, 0)
    try:
      w = wcs.WCS(referenceheader[0:20])
    except:
      w = wcs.WCS(referenceheader[0:19])

      
    #Repeat the same again
    pixcrd_0 = np.array([[x0, y0]], np.float_)
    world_0 = w.wcs_pix2world(pixcrd_0, 1)
    RA0 = world_0[0][0]
    DEC0 = world_0[0][1]
    
    if reverse==False:
        world_1 = np.array([[RA0, DEC0+length_of_arrow_arc/3600.]], np.float_)
    else:
        world_1 = np.array([[RA0, DEC0-length_of_arrow_arc/3600.]], np.float_)
    pixcrd_1 = w.wcs_world2pix(world_1, 1)
    xN = pixcrd_1[0][0]
    yN = pixcrd_1[0][1]   
    
    dx_N = xN - x0
    dy_N = yN - y0       


    
    world_2 = np.array([[RA0+length_of_arrow_arc/3600., DEC0]], np.float_)
    pixcrd_2 = w.wcs_world2pix(world_2, 1)
    xE = pixcrd_2[0][0]
    yE = pixcrd_2[0][1]
    dx_E = xE - x0
    dy_E = yE - y0



    ax.arrow(x0, y0, dx_N, dy_N,head_width=length_of_arrow/20., head_length=length_of_arrow/10., color='red')
    ax.text(xN-length_of_arrow/10., yN, 'N', color='red',fontsize=8, horizontalalignment='right', verticalalignment='bottom')

    dx_E = np.sign(dx_E) * sqrt(dx_N**2+dy_N**2-dy_E**2)
    ax.arrow(x0, y0, dx_E, dy_E,head_width=length_of_arrow/20., head_length=length_of_arrow/10., color='red')
    xE = x0 + dx_E
    if reverse==False:
        ax.text(xE, yE+length_of_arrow/20., 'E', color='red',fontsize=8, horizontalalignment='right', verticalalignment='bottom')
    else:
        ax.text(xE, yE-length_of_arrow/20., 'E', color='red',fontsize=8, horizontalalignment='right', verticalalignment='bottom')
'''

def find_angle(input_image):
        hdulist = pyfits.open(input_image)
        header = hdulist[0].header
        data = hdulist[0].data
        ny,nx = np.shape(data)
        
        if 'COMMENT' in header:
            del header['COMMENT']
        if 'HISTORY' in header:
            del header['HISTORY']
        if '' in header:
            del header[''] 

        w = wcs.WCS(header)

        #N
        pixcrd = np.array([[0.,0.]], np.float_)
        world = w.wcs_pix2world(pixcrd, 1)

        RA_0 = world[0][0]
        DEC_0 = world[0][1]

        pixcrd = np.array([[0.,ny]], np.float_)
        world = w.wcs_pix2world(pixcrd, 1)

        RA_1 = world[0][0]
        DEC_1 = world[0][1]
        
        vector_image = np.array([RA_1-RA_0,DEC_1-DEC_0])
        vector_wcs = np.array([0.,1.])
        angle_N = py_ang(vector_image, vector_wcs)

        #E
        pixcrd = np.array([[0.,0.]], np.float_)
        world = w.wcs_pix2world(pixcrd, 1)
        RA_0 = world[0][0]
        DEC_0 = world[0][1]

        pixcrd = np.array([[-nx,0.]], np.float_)
        world = w.wcs_pix2world(pixcrd, 1)

        RA_1 = world[0][0]
        DEC_1 = world[0][1]
        
        vector_image = np.array([RA_1-RA_0,DEC_1-DEC_0])
        vector_wcs = np.array([1.,0.])
        angle_E = py_ang(vector_image, vector_wcs)
        #print angle_E
        #print angle_N,angle_E

        if fabs(angle_N-angle_E)<180.:
          flip=False
        else:
          flip=True
        hdulist.close()
        return angle_N,flip

def axes_direct(input_fits_image, input_image, output_image):
    angle,flip = find_angle(input_fits_image)
    
    
    fig =  plt.figure(0,figsize=(5,5))
    ax = plt.axes()
    if flip==False:
      ax.arrow(0, 0, -0.5, 0, head_width=0.05, head_length=0.1, fc='lime', ec='lime',lw=3)
      if angle>90. and angle<270.:
        ax.text(-0.67, -0.05, r'$\exists$', color='lime',fontsize=50)
      else:
        ax.text(-0.67, -0.05, 'E', color='lime',fontsize=50)
      ax.arrow(0, 0, 0, 0.5, head_width=0.05, head_length=0.1, fc='lime', ec='lime',lw=3)
      ax.text(-0.04, 0.61, 'N', color='lime',fontsize=50)
      xlim(-0.7,0.05)
      ylim(-0.1,0.7)
    else:
      ax.arrow(0, 0, +0.5, 0, head_width=0.05, head_length=0.1, fc='lime', ec='lime',lw=3)
      if angle>90. and angle<270.:
        ax.text(-0.67, -0.05, r'$\exists$', color='lime',fontsize=50)
      else:
        ax.text(-0.67, -0.05, 'E', color='lime',fontsize=50)
      ax.arrow(0, 0, 0, 0.5, head_width=0.05, head_length=0.1, fc='lime', ec='lime',lw=3)
      ax.text(-0.04, 0.61, 'N', color='lime',fontsize=50)
      xlim(-0.02,0.7)
      ylim(-0.0,0.7)
      
    plt.axis('off')

    plt.savefig('EN.png', bbox_inches='tight', pad_inches=0,transparent=True)
    plt.clf()
    plt.close()
    #exit()

    img = Image.open('EN.png')
    im = img.rotate(angle, expand=2)
    im.save('EN_rot.png') 

    background = Image.open(input_image)
    width, height = background.size

    basewidth = int(ceil(min([width,height]) / 5.))
    img = Image.open('EN_rot.png')
    wpercent = (basewidth/float(img.size[0]))
    hsize = int((float(img.size[1])*float(wpercent)))
    img = img.resize((basewidth,hsize), PIL.Image.ANTIALIAS)
    img.save('EN1.png') 

    foreground = Image.open('EN1.png')
    width_foreground, height_foreground = foreground.size
    width_background, height_background = background.size

    background.paste(foreground, (0,height_background-height_foreground), foreground)
    #background.show()

    background.save(output_image)#, "PNG", dpi=(dpi,dpi)) 
    os.remove('EN1.png')
    os.remove('EN_rot.png')
    os.remove('EN.png')

def main(input_file, m0, pix2sec, mu_min, mu_max, output_file=None, text=None):
    vmin_faint = 10**(0.4*(m0 - mu_max)) * (pix2sec**2)
    vmax_faint = 10**(0.4*(m0 - mu_max+1)) * (pix2sec**2)

    vmax_bright = 10**(0.4*(m0 - mu_min)) * (pix2sec**2)
    vmin_bright = vmax_faint
    #print vmin_faint,vmax_faint,vmax_bright,vmin_bright
    
    if output_file==None:
        output_file = input_file.split('/')[-1].split('.fits')[0]+'.png'
    hdulist1 = pyfits.open(input_file)
    data = hdulist1[0].data
    ny,nx = np.shape(data)

    fig =  plt.figure(0)
    ax=plt.subplot()
    #fig.patch.set_alpha(0)

    img2 = ax.imshow(data,cmap='Blues', norm=LogNorm(vmin=vmin_faint, vmax=vmax_faint),interpolation='bicubic', origin="lower", aspect='equal',alpha=.5)
    img3 = ax.imshow(data,cmap='gray_r', norm=LogNorm(vmin=vmin_bright, vmax=vmax_bright),interpolation='bicubic', origin="lower", aspect='equal',alpha=.4) 
    if text!=None:
        plt.text(0.05, 0.93, text, fontsize=fsize, color='blue',transform=ax.transAxes)

    

    ax.set_frame_on(False)
    ax.set_xticks([]); ax.set_yticks([])
    plt.axis('off')
    crea_scale_bar(ax,nx,ny,pix2sec,l_bar=30.,kpc_per_arc=0.329)
    crea_scale_bar(ax,nx,ny,pix2sec,l_bar=30.,kpc_per_arc=None)
    #axes_direct(ax,input_file,nx,ny,pix2sec)
    plt.savefig(output_file, bbox_inches='tight', pad_inches=0)
    plt.clf()
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky background estimation")
    parser.add_argument("inputImage", help="Input fits image")
    parser.add_argument("--o", nargs='?', const=1, help="Optional: Input the name of the output file",type=str,default=None)
    parser.add_argument("--t", nargs='?', const=1, help="Optional: Input the text",type=str,default=None)
    args = parser.parse_args()
    input_image = args.inputImage
    output_image = args.o
    text = args.t
    main(input_image, m0=29.38, pix2sec=0.834, mu_min=20., mu_max=30., output_file=output_image, text=text)
