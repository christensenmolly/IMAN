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
import pyfits
import argparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy import wcs
warnings.filterwarnings("ignore")
from matplotlib.patches import Ellipse
import numpy.linalg as la
import PIL
from PIL import Image
#import Image
from astropy import wcs
import pyfits
fsize = 25
x_off = 0.95
y_off = 0.05


 
#@xl_func("numpy_row v1, numpy_row v2: float")
def py_ang(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'    """
    cosang = np.dot(v1, v2)
    sinang = la.norm(np.cross(v1, v2))
    return degrees(np.arctan2(sinang, cosang))

def roundup(x):
    return int(math.ceil(x / 30.0)) * 30

def crea_scale_bar(ax,nx,ny,pix2sec,l_bar=30.,kpc_per_arc=None,angle_units='arcsec'):
    if l_bar==None:
      l_bar = roundup(nx/1.375/10.)	# in arcsec
      if l_bar>=60.:
	l_bar = round(l_bar / 60.,1)
	angle_units = 'arcmin'
	
    if kpc_per_arc==None:
      xc = nx*x_off - l_bar/(2.*pix2sec)
      yc = ny*y_off
      if angle_units=='arcmin':
	# l_bar in arcses
	L_bar = l_bar * 60.
	ax.text(xc, yc, str(int(l_bar))+'\'', color='red',fontsize=8, horizontalalignment='center', verticalalignment='bottom',backgroundcolor='white', zorder=1)
	ax.errorbar(xc, yc, xerr=L_bar/(2.*pix2sec),color='lime',capsize=2,c='red', zorder=2)
      else:
	# l_bar in arcses
	ax.text(xc, yc, str(int(l_bar))+'\"', color='red',fontsize=8, horizontalalignment='center', verticalalignment='bottom',backgroundcolor='white', zorder=1)
	ax.errorbar(xc, yc, xerr=l_bar/(2.*pix2sec),color='lime',capsize=2,c='red', zorder=2)
    else:
      yc = ny*(1.-y_off)
      # l_bar in kpc
      L_bar = l_bar / kpc_per_arc	# Now in arcses
      xc = nx*x_off - L_bar/(2.*pix2sec)
      ax.text(xc, yc, str(int(l_bar))+'kpc', color='red',fontsize=8, horizontalalignment='center', verticalalignment='bottom',backgroundcolor='white', zorder=1)
      ax.errorbar(xc, yc, xerr=L_bar/(2.*pix2sec),color='lime',capsize=2,c='red', zorder=2)
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
      ax.arrow(0, 0, -0.5, 0, head_width=0.05, head_length=0.1, fc='r', ec='r',lw=2)
      if angle>90. and angle<270.:
	ax.text(-0.67, -0.05, r'$\exists$', color='red',fontsize=50)
      else:
	ax.text(-0.67, -0.05, 'E', color='red',fontsize=50)
      ax.arrow(0, 0, 0, 0.5, head_width=0.05, head_length=0.1, fc='r', ec='r',lw=2)
      ax.text(-0.04, 0.61, 'N', color='red',fontsize=50)
      xlim(-0.7,0.05)
      ylim(-0.1,0.7)
    else:
      ax.arrow(0, 0, +0.5, 0, head_width=0.05, head_length=0.1, fc='r', ec='r',lw=2)
      if angle>90. and angle<270.:
	ax.text(-0.67, -0.05, r'$\exists$', color='red',fontsize=50)
      else:
	ax.text(-0.67, -0.05, 'E', color='red',fontsize=50)
      ax.arrow(0, 0, 0, 0.5, head_width=0.05, head_length=0.1, fc='r', ec='r',lw=2)
      ax.text(-0.04, 0.61, 'N', color='red',fontsize=50)
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



    background.paste(foreground, (0,0), foreground)
    #background.show()
    background.save(output_image) 
    os.remove('EN1.png')
    os.remove('EN_rot.png')
    os.remove('EN.png')

def main(input_file, galaxy_image, m0, pix2sec, mu_min, mu_max, output_file=None, text=None,ellipse=None):
    vmin_faint = 10**(0.4*(m0 - mu_max)) * (pix2sec**2)
    vmax_faint = 10**(0.4*(m0 - mu_max+1)) * (pix2sec**2)

    vmax_bright = 10**(0.4*(m0 - mu_min)) * (pix2sec**2)
    vmin_bright = vmax_faint
    #print vmin_faint,vmax_faint,vmax_bright,vmin_bright
    
    if output_file==None:
        output_file = input_file.split('/')[-1].split('.fits')[0]+'.png'

    hdulist = pyfits.open(galaxy_image)
    image = hdulist[0].data

    hdulist1 = pyfits.open(input_file)
    data = hdulist1[0].data
    ny,nx = np.shape(data)

    fig =  plt.figure(0, figsize=(5, 5*ny/nx))
    ax=plt.subplot()

    data = np.fabs(data)
    image = np.fabs(image)    


    img3 = ax.imshow(image,cmap='Blues', norm=LogNorm(vmin=vmin_bright, vmax=vmax_bright),interpolation='bicubic', origin="lower", aspect='equal',alpha=.4) 

    #img2 = ax.imshow(data,cmap='Blues', norm=LogNorm(vmin=vmin_faint, vmax=vmax_faint),interpolation='bicubic', origin="lower", aspect='equal',alpha=.5)
    img3 = ax.imshow(data,cmap='gray', norm=LogNorm(vmin=vmin_bright, vmax=vmax_bright),interpolation='bicubic', origin="lower", aspect='equal',alpha=.7) 
    if text!=None:
        #plt.text(0.05, 0.93, text, fontsize=fsize, color='red',transform=ax.transAxes)
        plt.figtext(0.5, 0.85, text,color='red', horizontalalignment='center',verticalalignment='center',backgroundcolor='white')

    if ellipse!=None:
      [xc, yc, SMA, SMB, PA] = ellipse
      ellipse = Ellipse(xy=(xc, yc), width=SMA*2., height=SMB*2., angle=PA+90.,
			      edgecolor='r', fc='None', lw=2, alpha = 0.5)
      ax.add_patch(ellipse)


    ax.set_frame_on(False)
    ax.set_xticks([]); ax.set_yticks([])
    plt.axis('off')
    #crea_scale_bar(ax,nx,ny,pix2sec,l_bar=30.,kpc_per_arc=0.329)
    crea_scale_bar(ax,nx,ny,pix2sec,l_bar=None,kpc_per_arc=None)
    #axes_direct(ax,input_file,nx,ny,pix2sec)
    plt.savefig(output_file, bbox_inches='tight', pad_inches=0)
    plt.clf()
    plt.close()

'''
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
'''