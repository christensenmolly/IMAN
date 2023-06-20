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
import pyregion
from mpl_toolkits.axes_grid1 import make_axes_locatable

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
#matplotlib.use('agg',warn=False, force=True)   ####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

warnings.filterwarnings("ignore")

#***************************
fsize = 4
#***************************

'''
def discrete_cmap(N=8):
    # define individual colors as hex values
    cpool = [ '#000000', '#00EE00', '#0000EE', '#00EEEE', '#EE0000','#FFFF00', '#EE00EE', '#FFFFFF']
    cmap_i8 = mpl.colors.ListedColormap(cpool[0:N], 'i8')
    mpl.cm.register_cmap(cmap=cmap_i8)
    return cmap_i8
'''
def discrete_cmap(N=8, base_cmap='nipy_spectral'):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


def define_scale_bar_length(x_extent,pix2sec):
    scale_bar = round((x_extent * pix2sec) / 12.,0)
    #scale_bar_new = int(5. * round(float(scale_bar)/5.)) # Length of the bar in arcsec --- round to the nearest 5
    scale_bar_new = int(math.ceil(float(scale_bar) / 10.0)) * 10  # Length of the bar in arcsec --- round to the nearest 10
    if scale_bar_new>=60.:
      scale_bar_new=60.

    return scale_bar_new

def crea_scale_bar(ax,x0,x1,y0,y1,pix2sec,scale_bar,show_bar_label):
  offset_x_factor = 0.98
  offset_y_factor = 0.1 
  x_extent = x1 - x0

  #scale_bar_length = define_scale_bar_length(x_extent,pix2sec) / 2. #### divide by 2 !!!
  scale_bar_length = scale_bar / 2.

  xc = fabs(x1)-scale_bar_length/pix2sec - (1.-offset_x_factor)*(x1-x0)
  yc = fabs(y0) + (y1-y0)* offset_y_factor
  ax.errorbar(xc, yc, xerr=scale_bar_length/pix2sec,color='black',capsize=1,c='black')
  
  if x_extent>100:
      fsize = 8
  else:
      fsize = 4
  #ax.text(xc, yc, str(int(scale_bar_length*2.))+'\"', color='black',fontsize=fsize+1, horizontalalignment='center', verticalalignment='bottom')
  if show_bar_label=='yes':
    if scale_bar>=60.:
	  ax.text(xc, yc, str(int(1))+'\'', color='black',fontsize=fsize+1, horizontalalignment='center', verticalalignment='bottom')
    else:
	  ax.text(xc, yc, str(int(scale_bar_length*2.))+'\"', color='black',fontsize=fsize+1, horizontalalignment='center', verticalalignment='bottom')
    

def plot_indiv_frame(grid,k,number_of_comps,composed_model_file,mask_file,frame,m0,pix2sec,borders,min_level,color,scale_bar,view,show_bar_label,region,data_MASK=None):
        if frame=='Data':
            layer=(4+number_of_comps)*k+0
            numb_of_grid = 3*k + 0
        if frame=='Model':
            layer=(4+number_of_comps)*k+1
            numb_of_grid = 3*k + 1
        if frame=='Residual':
            layer=(4+number_of_comps)*k+3
            numb_of_grid = 3*k + 2
	
	hdulist = pyfits.open(composed_model_file)
	data = hdulist[layer].data

        xaxis = data.shape[1]
        yaxis = data.shape[0]
        
	if mask_file!=None:
	  hdulist_mask = pyfits.open(mask_file)
	  data_mask = hdulist_mask[0].data
          for k in range(yaxis):
              for i in range(xaxis):
		if data_mask[k,i]!=0.:
		  data[k,i] = 0.
	  

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

        if frame=='Data':
            data_MASK = np.ones((yaxis, xaxis))
            for k in range(yaxis):
                for i in range(xaxis):
                    if data[k,i]==0.:
                        data_MASK[k,i] = 0.
     
        
        if True:#frame=='Data': ###TODO:!!!
            mod_data = fabs(data)
            min_data = np.min(mod_data[np.nonzero(mod_data)])
            for k in range(yaxis):
                for i in range(xaxis):
                    if data[k,i]<=0.:
                        data[k,i]=min_data

	vmax_int = np.max(data)
        if min_level==None:
            vmin_int = np.min(data)
        else:
            vmin_int = 10**(0.4*(m0-min_level)) * (pix2sec**2)
	#print vmin_int,vmax_int
	#print data_MASK
	if frame!='Residual':
	    im = grid[numb_of_grid].imshow(data,cmap=color, norm=LogNorm(vmin=vmin_int, vmax=vmax_int),interpolation='none', origin="lower", aspect='equal') 
	    if True:#numb_of_grid == 0:
                im = grid[numb_of_grid].imshow(data_MASK,cmap='Reds',interpolation='none', origin="lower", aspect='equal', alpha=.1)
        else:
	  for k in range(yaxis):
             for i in range(xaxis):
                  if data_MASK[k,i]==0.:
                        data[k,i] = float('nan')
          cmap=discrete_cmap()
          cmap.set_bad((0, 0, 0, 0))
          
	  if view=='column':
            im = grid[numb_of_grid].imshow(data,cmap=cmap, vmin=0.001, vmax=1,interpolation='none', origin="lower", aspect='equal')
            cb = grid[numb_of_grid].cax.colorbar(im)
            for cax in grid.cbar_axes:
                cax.toggle_label(True)
                cax.axis[cax.orientation].set_label(' ')
                if view=='line':
		  cax.tick_params(labelsize=4)
		else:
		  cax.tick_params(labelsize=12)
                cax.set_xlim(0,1)
	  if view=='line':                
            im = grid[numb_of_grid].imshow(data,cmap=cmap, vmin=0.001, vmax=1,interpolation='none', origin="lower", aspect='equal')
            cb = grid[numb_of_grid].cax.colorbar(im)
            for cax in grid.cbar_axes:
                cax.toggle_label(True)
                cax.axis[cax.orientation].set_label(' ')
                cax.tick_params(labelsize=4)
                cax.set_ylim(0,1)

	if region!=None:            
            r2 = pyregion.open(region).as_imagecoord(hdulist[layer].header)
            patch_list, artist_list = r2.get_mpl_patches_texts()
            for p in patch_list:
                grid[numb_of_grid].add_patch(p)
            for t in artist_list:
                grid[numb_of_grid].add_artist(t)            
	
	if numb_of_grid==0 and scale_bar!=None:    
            crea_scale_bar(grid[numb_of_grid],x0,x1,y0,y1,pix2sec,float(scale_bar),show_bar_label)

        if frame=='Data':
            head = hdulist[layer].header
            grid[numb_of_grid].text(x0+(x1-x0)/10, y0+0.9*(y1-y0), head["BAND"], color='red',fontweight='bold',fontsize=fsize+3)#, backgroundcolor='white', alpha=.5)
        if frame!='Data': 
            return x0,x1,y0,y1,None,data_MASK
        else:
            return x0,x1,y0,y1,vmin_int,data_MASK

def main(composed_model_file,min_level,pix2sec,m0,mask_file=None,
	 borders='0,0,0,0',output_file='plot_2d.png',view='line',color='nipy_spectral_r',scale_bar=None,grid=None,show_bar_label='yes', region=None):
    no_fig = True
    if True:
        hdulist = pyfits.open(composed_model_file)
        header = hdulist[0].header

	if grid==None:
	  fig = plt.figure(0, (8, 8))
	  fig.subplots_adjust(left=0.05, right=0.95)
	  no_fig = False
	  
	  NUMBER_OF_BANDS = int(header['number_of_bands'])
	  NUMBER_OF_COMPS = int(header['number_of_comps'])

	  if view=='line':
	      grid = AxesGrid(fig, 111,
			      nrows_ncols=(NUMBER_OF_BANDS, 3),
			      axes_pad=0.02,
			      label_mode="L",
			      share_all=True,
			      cbar_location="right",
			      cbar_mode='single',
			      cbar_size="3%",
			      cbar_pad="1%",
			      )

        if len(min_level)==1:
            min_level = NUMBER_OF_BANDS * min_level

        if len(m0)==1:
            m0 = NUMBER_OF_BANDS * m0  

	for cax in grid.cbar_axes:
		cax.toggle_label(False)	    

        for k in range(NUMBER_OF_BANDS):
            # Reference image:
            x0,x1,y0,y1,vmin_int,data_MASK = plot_indiv_frame(grid,k,NUMBER_OF_COMPS,composed_model_file,mask_file,'Data',float(m0[k]),pix2sec,borders,float(min_level[k]),color,scale_bar,view,show_bar_label,region, data_MASK=None)

            # Model image:
            x0,x1,y0,y1,vmin_int,data_MASK = plot_indiv_frame(grid,k,NUMBER_OF_COMPS,composed_model_file,mask_file,'Model',float(m0[k]),pix2sec,borders,float(m0[k])-2.5*log10(vmin_int/(pix2sec**2)),color,scale_bar,view,show_bar_label,region,data_MASK=data_MASK)
            # Residual image:
            x0,x1,y0,y1,vmin_int,data_MASK = plot_indiv_frame(grid,k,NUMBER_OF_COMPS,composed_model_file,mask_file,'Residual',float(m0[k]),pix2sec,borders,float(min_level[k]),color,scale_bar,view,show_bar_label,region=None,data_MASK=data_MASK)

            # This affects all axes because we set share_all = True.
            grid.axes_llc.set_xlim(x0,x1)
            grid.axes_llc.set_ylim(y0,y1)
            grid.axes_llc.set_xticklabels([])
            grid.axes_llc.set_yticklabels([])
            grid.axes_llc.get_xaxis().set_ticks([])	# To remove ticks
            grid.axes_llc.get_yaxis().set_ticks([])	# To remove ticks


	if no_fig == False:
	  plt.draw()  
	  plt.savefig(output_file, bbox_inches='tight', pad_inches=0.02, dpi = 300)
	  plt.clf()
	  plt.close()
	  plt.close('all')
	  return output_file
    else:
	  plt.clf()
	  plt.close()
	  plt.close('all')        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create 2D image using the model and the reference image")
    parser.add_argument("model_file", help="Input model image")
    parser.add_argument("Scale", help="Input scale in [arcsec/pix]",type=str)
    parser.add_argument("m0", help="Input Zero Point in [mag/arcsec^2]",type=str)
    parser.add_argument("--min_level", nargs='?', const=0., help="Input minimum isophote level to highlight the structure in [mag/arcsec^2]",type=str,default=None)
    parser.add_argument("--mask_file", nargs='?', const=0., help="Name of the mask file",type=str,default=None)      
    parser.add_argument("--borders", nargs='?', const='0,0,0,0', help="Input the borders of the frame to be plotted: x1,y1,x2,y2",type=str,default='0,0,0,0')    
    parser.add_argument("--output_file", nargs='?', const=0., help="Name of the output file",type=str,default='plot_2d.png')      
    parser.add_argument("--view", nargs='?', const=0., help="Input the view of the plotted image: line or column",type=str,default='line')
    parser.add_argument("--color", nargs='?', const=0., help="Input the color for the plotted image: gray_r or nipy_spectral_r",type=str,default='nipy_spectral_r')
    parser.add_argument("--scale_bar", nargs='?', const=0., help="Create the scale bar: length [arcsec]",type=str,default=None)
    parser.add_argument("--show_bar_label", nargs='?', const=0., help="Show the scale bar label: yes",type=str,default='yes')
    parser.add_argument("--region", nargs='?', const=0., help="Region file",type=str,default=None)
    
    args = parser.parse_args()
    
    composed_model_file = args.model_file
    min_level = args.min_level
    pix2sec = float(args.Scale)
    m0 = args.m0
    mask_file = args.mask_file
    borders = str(args.borders)
    output_file = args.output_file
    view = str(args.view)
    color = str(args.color)
    scale_bar = args.scale_bar
    show_bar_label = args.show_bar_label
    region = args.region


    min_level = min_level.split(',')
    m0 = m0.split(',')   
    
    main(composed_model_file,min_level,pix2sec,m0,mask_file=mask_file,
	 borders=borders,output_file=output_file,view=view,color=color,scale_bar=scale_bar,show_bar_label=show_bar_label,region=region)
