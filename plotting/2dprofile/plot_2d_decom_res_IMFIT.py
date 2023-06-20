#! /usr/bin/env python
#!/usr/bin/python
# -*- coding:  cp1251 -*-
# USAGE:
# python ~/CurrentWork/ImaPrep/IMAN/Plot/plot_2d_decom_res_IMFIT.py composed_model.fits 24.5 0.6 21.581 0,0,0,0 pdf column gray nosub
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

warnings.filterwarnings("ignore")

#***************************
fsize = 4
#***************************

'''
def crea_scale_bar(ax,nx,ny,pix2sec):
  xc = offset_x_factor*(nx*pix2sec-scale_bar_length)
  yc = offset_y_factor*ny*pix2sec
  ax.errorbar(xc, yc, xerr=scale_bar_length,color='lime',capsize=2)
'''

def discrete_cmap(N=8):
    # define individual colors as hex values
    cpool = [ '#000000', '#00EE00', '#0000EE', '#00EEEE', '#EE0000','#FFFF00', '#EE00EE', '#FFFFFF']
    cmap_i8 = mpl.colors.ListedColormap(cpool[0:N], 'i8')
    mpl.cm.register_cmap(cmap=cmap_i8)
    return cmap_i8

def define_scale_bar_length(x_extent,pix2sec):
    scale_bar = round((x_extent * pix2sec) / 12.,0)
    #scale_bar_new = int(5. * round(float(scale_bar)/5.)) # Length of the bar in arcsec --- round to the nearest 5
    scale_bar_new = int(math.ceil(float(scale_bar) / 10.0)) * 10  # Length of the bar in arcsec --- round to the nearest 10
    if scale_bar_new>=60.:
      scale_bar_new=60.

    return scale_bar_new

def crea_scale_bar(ax,x0,x1,y0,y1,pix2sec,sub_scale_bar):
  offset_x_factor = 0.98
  offset_y_factor = 0.1 
  x_extent = x1 - x0

  scale_bar_length = define_scale_bar_length(x_extent,pix2sec) / 2. #### divide by 2 !!!

  xc = fabs(x1)-scale_bar_length/pix2sec - (1.-offset_x_factor)*(x1-x0)
  yc = fabs(y0) + (y1-y0)* offset_y_factor
  ax.errorbar(xc, yc, xerr=scale_bar_length/pix2sec,color='black',capsize=1,c='black')
  
  if x_extent>100:
      fsize = 8
  
  if sub_scale_bar=='sub':
    if scale_bar_length*2.>=60.:
      ax.text(xc, yc, str(int(1))+'\'', color='black',fontsize=fsize+1, horizontalalignment='center', verticalalignment='bottom')
    else:
      ax.text(xc, yc, str(int(scale_bar_length*2.))+'\"', color='black',fontsize=fsize+1, horizontalalignment='center', verticalalignment='bottom')
    

def plot_indiv_frame(grid,composed_model_file,frame,m0,pix2sec,borders,min_level,color,scale_bar,view):
	print view
        if frame=='Data':
            layer=0
            numb_of_grid = 0
        if frame=='Model':
            layer=1
            numb_of_grid = 1
        if frame=='Residual':
            layer=3
            numb_of_grid = 2
	hdulist = pyfits.open(composed_model_file)
	data = hdulist[layer].data
	data = data #/ sum(data)

        xaxis = data.shape[1]
        yaxis = data.shape[0]
        if borders!='0,0,0,0':    
            borders = borders.split(',')
            x0 = int(borders[0])
            y0 = int(borders[1])
            x1 = int(borders[2])
            y1 = int(borders[3])
        else:       
            x0 = 0
            y0 = 0
            x1 = xaxis
            y1 = yaxis
        
        
        if numb_of_grid == 0:
            mod_data = fabs(data)
            min_data = np.min(mod_data[np.nonzero(mod_data)])
            for k in range(yaxis):
                for i in range(xaxis):
                    if data[k,i]<=0.:
                        data[k,i]=min_data

	vmax = np.max(data) #np.mean([np.max(data_ski),np.max(data_ref)])
	vmin = np.min(data) #np.mean([np.min(data_ski),np.min(data_ref)])

	
	if numb_of_grid!=2:
            min_int = 10**(0.4*(m0-min_level)) * pix2sec #m0 - 2.5*log10(min_level) + 5.*log10(pix2sec)
            if color == 'gray':
	      im = grid[numb_of_grid].imshow(data,cmap='gray_r', norm=LogNorm(vmin=min_int, vmax=vmax),interpolation='none', origin="lower", aspect='equal')
	    else:
	      im = grid[numb_of_grid].imshow(data,cmap='nipy_spectral_r', norm=LogNorm(vmin=min_int, vmax=vmax),interpolation='none', origin="lower", aspect='equal')            
        else:
	  if view=='column':
            im = grid[numb_of_grid].imshow(data,cmap=discrete_cmap(), vmin=0.001, vmax=1,interpolation='none', origin="lower", aspect='equal')
            cb = grid[numb_of_grid].cax.colorbar(im)
            #cb.set_xticklabels(labelsize=1)
            #grid[numb_of_grid].cax.toggle_label(True)
            for cax in grid.cbar_axes:
                cax.toggle_label(True)
                cax.axis[cax.orientation].set_label(' ')
                #cax.axis[cax.orientation].set_fontsize(3)
                if view=='line':
		  cax.tick_params(labelsize=4)
		else:
		  cax.tick_params(labelsize=12)
                cax.set_xlim(0,1)
                #cax.set_yticklabels([0, 0.5, 1])
	  if view=='line':                
            im = grid[numb_of_grid].imshow(data,cmap=discrete_cmap(), vmin=0.001, vmax=1,interpolation='none', origin="lower", aspect='equal')
            cb = grid[numb_of_grid].cax.colorbar(im)
            #cb.set_xticklabels(labelsize=1)
            #grid[numb_of_grid].cax.toggle_label(True)
            for cax in grid.cbar_axes:
                cax.toggle_label(True)
                cax.axis[cax.orientation].set_label(' ')
                #cax.axis[cax.orientation].set_fontsize(3)
                cax.tick_params(labelsize=4)
                cax.set_ylim(0,1)
                #cax.set_yticklabels([0, 0.5, 1])
                
            

        #grid[numb_of_grid].text(0.03, 0.95,frame, color='lime',transform=grid[numb_of_grid].transAxes,fontsize=fsize+2, fontweight='bold', va='top')
	
	if numb_of_grid==0:    
            #crea_scale_bar(grid[numb_of_grid],x0,x1,y0,y1,pix2sec,scale_bar)
            crea_scale_bar(grid[numb_of_grid],x0,x1,y0,y1,pix2sec,scale_bar)
	return x0,x1,y0,y1


def main(composed_model_file,min_level,pix2sec,m0,borders='0,0,0,0',format_out='pdf',view='line',color='color',scale_bar='sub'):
	if view=='line':
	    fig = plt.figure(0, (8, 8))
	    fig.subplots_adjust(left=0.05, right=0.95) 

	    grid = AxesGrid(fig, 111,
			    nrows_ncols=(1, 3),
			    axes_pad=0.02,
			    label_mode="L",
			    share_all=True,
			    cbar_location="right",
			    cbar_mode='single',
			    cbar_size="5%",
			    cbar_pad="1%",
			    ) # cbar_mode="single"



	    for cax in grid.cbar_axes:
		cax.toggle_label(False)

	    # Reference image:
	    x0,x1,y0,y1 = plot_indiv_frame(grid,composed_model_file,'Data',m0,pix2sec,borders,min_level,color,scale_bar,view)
	    # Model image:
	    x0,x1,y0,y1 = plot_indiv_frame(grid,composed_model_file,'Model',m0,pix2sec,borders,min_level,color,scale_bar,view)
	    # Residual image:
	    x0,x1,y0,y1 = plot_indiv_frame(grid,composed_model_file,'Residual',m0,pix2sec,borders,min_level,color,scale_bar,view)

	    # This affects all axes because we set share_all = True.
	    grid.axes_llc.set_xlim(x0,x1)
	    grid.axes_llc.set_ylim(y0,y1)
	    grid.axes_llc.set_xticklabels([])
	    grid.axes_llc.set_yticklabels([])
	    grid.axes_llc.get_xaxis().set_ticks([])	# To remove ticks
	    grid.axes_llc.get_yaxis().set_ticks([])	# To remove ticks



	    plt.draw()  
	    plt.savefig('2d_decom_res.%s' % (format_out), bbox_inches='tight', pad_inches=0.02, dpi = 300)
	    plt.clf()
	    plt.close()   
	    #plt.show()
	elif view=='column':
	    fig = plt.figure(0, (8, 8))
	    fig.subplots_adjust(left=0.05, right=0.95) 

	    grid = AxesGrid(fig, 111,
			    nrows_ncols=(3, 1),
			    axes_pad=0.02,
			    label_mode="L",
			    share_all=True,
			    cbar_location="bottom",
			    cbar_mode='single',
			    cbar_size="8%",
			    cbar_pad="1%",
			    ) # cbar_mode="single"



	    for cax in grid.cbar_axes:
		cax.toggle_label(False)

	    # Reference image:
	    x0,x1,y0,y1 = plot_indiv_frame(grid,composed_model_file,'Data',m0,pix2sec,borders,min_level,color,scale_bar,view)
	    # Model image:
	    x0,x1,y0,y1 = plot_indiv_frame(grid,composed_model_file,'Model',m0,pix2sec,borders,min_level,color,scale_bar,view)
	    # Residual image:
	    x0,x1,y0,y1 = plot_indiv_frame(grid,composed_model_file,'Residual',m0,pix2sec,borders,min_level,color,scale_bar,view)

	    # This affects all axes because we set share_all = True.
	    grid.axes_llc.set_xlim(x0,x1)
	    grid.axes_llc.set_ylim(y0,y1)
	    grid.axes_llc.set_xticklabels([])
	    grid.axes_llc.set_yticklabels([])
	    grid.axes_llc.get_xaxis().set_ticks([])	# To remove ticks
	    grid.axes_llc.get_yaxis().set_ticks([])	# To remove ticks



	    plt.draw()  
	    plt.savefig('2d_decom_res.%s' % (format_out), bbox_inches='tight', pad_inches=0.02, dpi = 300)
	    plt.clf()
	    plt.close()   
	    #plt.show()
	else:
	    zz=1
	return '2d_decom_res.%s' % (format_out)
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create 2D image using the model and the reference image")
    parser.add_argument("model_file", help="Input model image")
    parser.add_argument("min_level", help="Input minimum isophote level to highlight the structure in [mag/arcsec^2]")
    parser.add_argument("Scale", help="Input scale in [arcsec/pix]")
    parser.add_argument("m0", help="Input Zero Point in [mag/arcsec^2]")     
    parser.add_argument("borders", nargs='?', const='0,0,0,0', help="Input the borders of the frame to be plotted: x1,y1,x2,y2",type=str,default='0,0,0,0')    
    parser.add_argument("format_out", nargs='?', const=0., help="Input the format of the output file",type=str,default='eps')      
    parser.add_argument("view", nargs='?', const=0., help="Input the view of the plotted image: line or column",type=str,default='line')
    parser.add_argument("color", nargs='?', const=0., help="Input the color for the plotted image: gray or color",type=str,default='color')
    parser.add_argument("scale_bar", nargs='?', const=0., help="Subscribe the scale bar: sub or nosub",type=str,default='sub')
    
    args = parser.parse_args()
    
    composed_model_file = args.model_file
    min_level = float(args.min_level)
    pix2sec = float(args.Scale)
    m0 = float(args.m0)
    borders = str(args.borders)
    format_out = str(args.format_out)
    view = str(args.view)
    color = str(args.color)
    scale_bar = str(args.scale_bar)
    
    main(composed_model_file,min_level,pix2sec,m0,borders,format_out,view,color,scale_bar)

