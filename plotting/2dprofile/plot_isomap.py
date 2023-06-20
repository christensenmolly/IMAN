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

from astropy.io import fits as fits
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



def main(input_image, m0, pix2sec, inner_iso, outer_iso, step=0.5, output_image = None, simple_map_color='white', no_fig=False, units='arcsec'):
      if step==None:
        step = (outer_iso - inner_iso) / 25.
      else:
        step = float(step)
      levels = np.arange(inner_iso, outer_iso, step)

      hdulist = fits.open(input_image)
      prihdr = hdulist[0].header
      data = hdulist[0].data

      ny,nx = np.shape(data)


      plt.figure()
      AX = plt.subplot()

      # Smoothing the data
      data = ndimage.gaussian_filter(data, sigma=1.0, order=0)


      for k in range(ny):
         for i in range(nx):
                if data[k,i]<=0.:
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
                im = plt.imshow(data,cmap='gray', vmin=vmin, vmax=vmax,interpolation='none', origin="lower", aspect='equal')
            else:
                im = plt.imshow(data,cmap=simple_map_color, vmin=vmin, vmax=vmax, interpolation='none', origin="lower", aspect='equal')




      if simple_map_color=='white':
                    CS = AX.contour(data, levels,
                                    origin='lower',
                                    linewidths=1,
                                    aspect='equal')

      else:
                    CS = AX.contour(data, levels,
                                    origin='lower',
                                    linewidths=1,
                                    aspect='equal')

      if True:
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
 




      #AX.set_xlim(x0,x1)
      #AX.set_ylim(y0,y1)  

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
            output_image = 'isomap.png'


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
    parser.add_argument("--step", nargs='?', const=1, help="Input step in [mag/arcsec^2]",type=float,default=0.5)
    #parser.add_argument("--xc", nargs='?', const=1, help="Optional: Input x-coordinate of the center [pixels]",type=float,default=0.)
    #parser.add_argument("--yc", nargs='?', const=1, help="Optional: Input y-coordinate of the center [pixels]",type=float,default=0.)     
    #parser.add_argument("--x0", nargs='?', const=1, help="Optional: Input x-coordinate of the left corner [pixels]",type=float,default=0.)
    #parser.add_argument("--y0", nargs='?', const=1, help="Optional: Input y-coordinate of the left corner [pixels]",type=float,default=0.)    
    #parser.add_argument("--x1", nargs='?', const=1, help="Optional: Input x-coordinate of the rigth corner [pixels]",type=float,default=0.)
    #parser.add_argument("--y1", nargs='?', const=1, help="Optional: Input y-coordinate of the rigth corner [pixels]",type=float,default=0.)    
    #parser.add_argument("--mask_image", nargs='?', help="Optional: Input mask image",default=None)
    parser.add_argument("--output_image", nargs='?', help="Optional: Input output image",default='isomap.png')
    #parser.add_argument("--simple_map_color", nargs='?', help="Optional: Input color",default=None)
    #parser.add_argument("--name", nargs='?', help="Optional: Input name of the object",default=None)
    #parser.add_argument("--scalebar", nargs='?', help="Optional: Input length of the scalebar [arcsec]",default=None)
    #parser.add_argument("--units", nargs='?', help="Optional: Input units [pix,arcsec,kpc]",type=str,default='pix')
    #parser.add_argument("--Scale", nargs='?', help="Input Scale in [kpc/arcsec]",type=float,default=0.1)  
    
    args = parser.parse_args()
    
    input_image = args.input_image
    m0 = float(args.m0)
    scale = float(args.scale)
    inner_iso = float(args.inner_iso)
    outer_iso = float(args.outer_iso)
    step = args.step
    #xc = float(args.xc)
    #yc = float(args.yc)    
    #x0 = float(args.x0)
    #x1 = float(args.x1)
    #y0 = float(args.y0)
    #y1 = float(args.y1)
    #mask_image = args.mask_image
    output_image = args.output_image
    #simple_map_color = args.simple_map_color
    #name = args.name
    #scalebar = args.scalebar
    #units = args.units
    #Scale = float(args.Scale)
    
    main(input_image, m0, scale, inner_iso, outer_iso, step, output_image=output_image)
    #, xc=xc, yc=yc,x0=x0,x1=x1,y0=y0,y1=y1,mask_image=mask_image,
    #AX=None, output_image = output_image, simple_map_color=simple_map_color, name=name, scalebar=scalebar,units=units,Scale=Scale)





