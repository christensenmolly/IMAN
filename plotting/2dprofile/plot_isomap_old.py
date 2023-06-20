#!/usr/bin/python
# -*- coding:  cp1251 -*-

# Import standard modules
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
import matplotlib.ticker as ticker
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
import re
import glob
import warnings
from scipy import special
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy import optimize
import numpy.ma as ma
import scipy.ndimage as ndimage

from astropy.io import fits as pyfits
from scipy import interpolate
import argparse

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

warnings.filterwarnings("ignore")




#***************************
offset_x_factor = 0.98
offset_y_factor = 0.1
fsize = 15
width_of_scalebar = 4
#***************************

def crea_scale_bar(AX,x0,x1,y0,y1,pix2sec,orientation,scale_bar_length = 30):
  if orientation=='N':
    xc = fabs(x1)-scale_bar_length/pix2sec - (1.-offset_x_factor)*(x1-x0)
    yc = fabs(y0) + (y1-y0)* offset_y_factor
  else:
    xc = fabs(x0)+scale_bar_length/pix2sec + (1.-offset_x_factor)*(x1-x0)
    yc = fabs(y1) - (y1-y0)* offset_y_factor  

  AX.errorbar(xc, yc, xerr=scale_bar_length/pix2sec, color='black',capsize=2,c='black', elinewidth=width_of_scalebar)


def myround(x, base=5):
    return int(base * round(float(x)/base))



def main(input_image,m0,pix2sec,inner_iso,outer_iso,step=None, xc=0.,yc=0.,x0=0.,x1=0.,y0=0.,y1=0.,mask_image=None, AX=None, output_image = None, simple_map_color=None, name=None, scalebar=None, units='pix', Scale=0.1):

     
      if mask_image!=None:
        hdulist_mask = pyfits.open(mask_image)
        data_mask = hdulist_mask[0].data 
        #mask_astropy = np.zeros_like(np.array(data_mask,dtype=float),dtype=bool)
        #data = ma.array(data, mask=mask_astropy)

      if step==None:
        step = (outer_iso - inner_iso) / 25.
      else:
        step = float(step)

      # Find number of layers:
      hdu = pyfits.open(input_image)
      number_of_layers = len(hdu)    

      if number_of_layers>=5:
        model_to_fil = hdu[1].data
      else:
        model_to_fil = hdu[0].data

      nx = model_to_fil.shape[1]
      ny = model_to_fil.shape[0]

      no_fig = True
      if AX==None:
        no_fig = False
        plt.figure()
        AX = plt.subplot()


      if x0==0. and y0==0. and x1==0. and y1==0.:
	      '''
	      x0 = nx/8.
	      x1 = 7.*nx/8.
	      y0 = ny/8.
	      y1 = 7.*ny/8.
	      '''
	      x0 = 1.
	      x1 = nx
	      y0 = 1.
	      y1 = ny
	      
      if xc==0. and yc==0.:
        xc = nx/2.
        yc = ny/2.

      if number_of_layers<5:
        layers = [0]
      elif number_of_layers==5:
        layers = [0,1]#[0,4,1]
      elif number_of_layers>5:
        layers = [0,1]#[0]+range(4,number_of_layers)+[1]

      for layer in layers:
        prihdr = hdu[layer].header
        try:
            Label = prihdr['NAME_OF_LAYER']
        except:
            Label = 'galaxy'

        data = hdu[layer].data

        if layer==0 and number_of_layers>=5:
            # Fill masked areas with the model
            for k in range(ny):
                for i in range(nx):
                    if data_mask[k,i]!=0.:
                        data[k,i]=model_to_fil[k,i]

        # Smoothing the data
        data = ndimage.gaussian_filter(data, sigma=1.0, order=0)


        for k in range(ny):
            for i in range(nx):
                if data[k,i]<=0.: #### FOR NGC100: or data[k,i]<=10**(0.4*(m0-25.5)) *(pix2sec**2):
                    data[k,i] = outer_iso
                else:
                    data[k,i] = m0 - 2.5*log10(data[k,i]) + 5.*log10(pix2sec)
    

        vmin = np.min(fabs(data))
        vmax = outer_iso

        if no_fig == False:
            if simple_map_color==None:
                #im = plt.imshow(data,cmap='gray', vmin=vmin, vmax=vmax, interpolation='none', origin="lower", aspect='equal')
                z=1
            elif simple_map_color=='white':
                im = plt.imshow(data,cmap='gray', vmin=vmin/100., vmax=vmin/10.,interpolation='none', origin="lower", aspect='equal')
        else:
            im = plt.imshow(data,cmap=simple_map_color, vmin=vmin, vmax=vmax, interpolation='none', origin="lower", aspect='equal')


        levels = np.arange(inner_iso, outer_iso, step)           # levels

        if layer==0 and number_of_layers>=5:
            if layer==0:
                CS = plt.contour(data, levels,colors='black',
                    origin='lower',             # origin in lower-left corner
                    linewidths=1,               # line width
                    aspect='equal')         	# outer pixel boundaries
            else:
                CS = plt.contour(data, levels,colors='red',
                    origin='lower',             # origin in lower-left corner
                    linewidths=1,               # line width
                    aspect='equal')         	# outer pixel boundaries 
        else:
            CS = plt.contour(data, levels,
                    origin='lower',             # origin in lower-left corner
                    linewidths=1,               # line width
                    aspect='equal')         	# outer pixel boundaries
        if layer==0 and number_of_layers<5:

            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(AX)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            CB = plt.colorbar(CS, cax=cax)
            CB.set_ticks(levels[0::3])
            CB.set_ticklabels(levels[0::3])
  
            if name!=None:
                props = dict(facecolor='white', alpha=1.0)
                plt.text(0.03, 0.95, name, transform=AX.transAxes,fontsize=fsize+2, fontweight='bold', va='top', bbox=props)
            AX.axis('equal')
        else:
            if simple_map_color==None:
                #im = AX.imshow(data,cmap='gray', vmin=vmin, vmax=vmax, interpolation='none', origin="lower", aspect='equal')
                z=1
            elif simple_map_color=='white':
                im = AX.imshow(data,cmap='gray', vmin=vmin/100., vmax=vmin/10.,interpolation='none', origin="lower", aspect='equal')
            else:
                im = AX.imshow(data,cmap=simple_map_color, vmin=vmin, vmax=vmax, interpolation='none', origin="lower", aspect='equal')

        levels = np.arange(inner_iso, outer_iso, step)           # levels

        if simple_map_color=='white':
            if number_of_layers>=5:
                if layer==0:
                    CS = AX.contour(data, levels,colors='black',
					    origin='lower',             # origin in lower-left corner
					    linewidths=1,               # line width
					    aspect='equal')         	# outer pixel boundaries
                else:
                    CS = AX.contour(data, levels,colors='red',
					    origin='lower',             # origin in lower-left corner
					    linewidths=1,               # line width
					    aspect='equal')         	# outer pixel boundaries   
            else: 
                CS = AX.contour(data, levels, colors='black',
					origin='lower',             # origin in lower-left corner
					linewidths=1,               # line width
					aspect='equal')         	# outer pixel boundaries
		else:
            CS = AX.contour(data, levels,
				      origin='lower',             # origin in lower-left corner
				      linewidths=1,               # line width
				      aspect='equal')         	# outer pixel boundaries

		if layer==0:
              from mpl_toolkits.axes_grid1 import make_axes_locatable

              if simple_map_color!='white':
                divider = make_axes_locatable(AX)
                cax = divider.append_axes("right", size="5%", pad=0.05)			
                CB = plt.colorbar(CS, cax=cax)
                CB.set_label(r' $\mu$ (mag arcsec$^{-2}$) ')
                CB.set_ticks(levels)
                CB.set_ticklabels(levels)
                if units=='pix':
                    AX.set_xlabel(r'X (pix)')
                    AX.set_ylabel(r'Y (pix)')
                elif units=='arcsec':
                    AX.set_xlabel(r'X (arcsec)')
                    AX.set_ylabel(r'Y (arcsec)')
                elif units=='kpc':
                    AX.set_xlabel(r'X (kpc)')
                    AX.set_ylabel(r'Y (kpc)')

              AX.tick_params(axis='both',color='black',which='major',length=7)
              AX.tick_params(axis='both',color='black',which='minor',length=4)
              AX.axis('equal')
 


      if scalebar!=None:
        crea_scale_bar(AX,x0,x1,y0,y1,pix2sec,'N',scale_bar_length=float(scalebar))

      AX.set_xlim(x0,x1)
      AX.set_ylim(y0,y1)  

      if units=='arcsec':
        ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(myround((x-xc)*pix2sec*Scale, base=10)))
        AX.xaxis.set_major_formatter(ticks_x)

        ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(myround((x-xc)*pix2sec*Scale, base=10)))
        AX.yaxis.set_major_formatter(ticks_y)

      if units=='kpc':
        ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(myround((x-xc)*pix2sec*Scale, base=5)))
        AX.xaxis.set_major_formatter(ticks_x)

	ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(myround((x-xc)*pix2sec*Scale, base=5)))
	AX.yaxis.set_major_formatter(ticks_y)

      if no_fig == False:
            AX.get_xaxis().set_ticks([])	# To remove ticks
            AX.get_yaxis().set_ticks([])	# To remove ticks
    if output_image == None:
            output_image = 'isomap.eps'


    plt.savefig(output_image, transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.clf()
    plt.close()
    return output_image
      


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create isophote map")
    parser.add_argument("input_image", help="Input image")
    parser.add_argument("m0", help="Input Zero Point in [mag/arcsec^2]")     
    parser.add_argument("scale", help="Input scale in [arcsec/pix]")  
    parser.add_argument("inner_iso", help="Input inner isophote level in [mag/arcsec^2]")
    parser.add_argument("outer_iso", help="Input outer isophote level in [mag/arcsec^2]")
    parser.add_argument("--step", nargs='?', const=1, help="Input step in [mag/arcsec^2]",type=float,default=None)
    parser.add_argument("--xc", nargs='?', const=1, help="Optional: Input x-coordinate of the center [pixels]",type=float,default=0.)
    parser.add_argument("--yc", nargs='?', const=1, help="Optional: Input y-coordinate of the center [pixels]",type=float,default=0.)     
    parser.add_argument("--x0", nargs='?', const=1, help="Optional: Input x-coordinate of the left corner [pixels]",type=float,default=0.)
    parser.add_argument("--y0", nargs='?', const=1, help="Optional: Input y-coordinate of the left corner [pixels]",type=float,default=0.)    
    parser.add_argument("--x1", nargs='?', const=1, help="Optional: Input x-coordinate of the rigth corner [pixels]",type=float,default=0.)
    parser.add_argument("--y1", nargs='?', const=1, help="Optional: Input y-coordinate of the rigth corner [pixels]",type=float,default=0.)    
    parser.add_argument("--mask_image", nargs='?', help="Optional: Input mask image",default=None)
    parser.add_argument("--output_image", nargs='?', help="Optional: Input output image",default=None)
    parser.add_argument("--simple_map_color", nargs='?', help="Optional: Input color",default=None)
    parser.add_argument("--name", nargs='?', help="Optional: Input name of the object",default=None)
    parser.add_argument("--scalebar", nargs='?', help="Optional: Input length of the scalebar [arcsec]",default=None)
    parser.add_argument("--units", nargs='?', help="Optional: Input units [pix,arcsec,kpc]",type=str,default='pix')
    parser.add_argument("--Scale", nargs='?', help="Input Scale in [kpc/arcsec]",type=float,default=0.1)  
    
    args = parser.parse_args()
    
    input_image = args.input_image
    m0 = float(args.m0)
    scale = float(args.scale)
    inner_iso = float(args.inner_iso)
    outer_iso = float(args.outer_iso)
    step = args.step
    xc = float(args.xc)
    yc = float(args.yc)    
    x0 = float(args.x0)
    x1 = float(args.x1)
    y0 = float(args.y0)
    y1 = float(args.y1)
    mask_image = args.mask_image
    output_image = args.output_image
    simple_map_color = args.simple_map_color
    name = args.name
    scalebar = args.scalebar
    units = args.units
    Scale = float(args.Scale)
    
    main(input_image,m0,scale,inner_iso,outer_iso,step, xc=xc, yc=yc,x0=x0,x1=x1,y0=y0,y1=y1,mask_image=mask_image,
    AX=None, output_image = output_image, simple_map_color=simple_map_color, name=name, scalebar=scalebar,units=units,Scale=Scale)





