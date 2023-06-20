#! /usr/bin/env python


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
import astropy.io.fits as pyfits
import argparse
try:
    import pyregion
except:
    z=1
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.stats import sigma_clipped_stats
import scipy.ndimage as ndimage
import warnings
from photutils import make_source_mask
from astropy.stats import sigma_clipped_stats



matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
#matplotlib.use('agg',warn=False, force=True)   ####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

LOCAL_DIR = "/plotting/2dprofile"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/rotate'))

import rotate_image

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


def discrete_cmap_neg(N=17, base_cmap='nipy_spectral'):
    from matplotlib.colors import BoundaryNorm
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:
    
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)
    '''
    # define the colormap
    cmap = plt.cm.get_cmap(base_cmap)

    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize and forcing 0 to be part of the colorbar!
    bounds = np.arange(-1,1,N)
    idx=np.searchsorted(bounds,0)
    bounds=np.insert(bounds,idx,0)
    norm = BoundaryNorm(bounds, cmap.N)
    return cmap,norm
    '''
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
    

def plot_indiv_frame(grid,composed_model_file,mask_file,frame,m0,pix2sec,borders,min_level,color,scale_bar,view,show_bar_label,region,text,data_MASK=None, show_negative=False, sigma=None, rotate_PA=0.):

        if frame=='Data':
            layer=0
            numb_of_grid = 0
        if frame=='Model':
            layer=1
            numb_of_grid = 1
        if frame=='Residual':
            layer=3
            numb_of_grid = 2
        
        if rotate_PA!=0.:
            rotate_image.main(composed_model_file, rotate_PA, xc=None, yc=None, output_image='tmp_rot.fits', hdu_inp=layer, cval=float('nan'), cropping=True, verbosity=True)
            hdulist = pyfits.open('tmp_rot.fits')
            data = hdulist[0].data
            os.remove('tmp_rot.fits')
            
        else:
            hdulist = pyfits.open(composed_model_file)
            data = hdulist[layer].data
        
        
        if sigma is not None:
            data = ndimage.gaussian_filter(data, sigma=sigma, order=0)


        xaxis = data.shape[1]
        yaxis = data.shape[0]
        
        if mask_file!=None:
          if rotate_PA==0.:
                hdulist_mask = pyfits.open(mask_file)
                data_mask = hdulist_mask[0].data
          else:
                rotate_image.main(mask_file, rotate_PA, xc=None, yc=None, output_image='tmp_mask.fits', hdu_inp=0, cval=0., cropping=True, verbosity=True)
                hdulist_mask = pyfits.open('tmp_mask.fits')
                data_mask = hdulist_mask[0].data
                os.remove('tmp_mask.fits')              
              
              
          #data_MASK = np.ones((yaxis, xaxis))
          
          
          for k in range(yaxis):
              for i in range(xaxis):
                if data_mask[k,i]!=0.:
                  data[k,i] = 0.
                  #data_MASK[k,i] = 0.
          

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
            data_MASK = np.ones((yaxis, xaxis))
            mod_data = fabs(data)
            min_data = np.min(mod_data[np.nonzero(mod_data)])
            for k in range(yaxis):
                for i in range(xaxis):
                    if data[k,i]==0.:
                        data_MASK[k,i] = 0.
                    if data[k,i]<=0.:
                        data[k,i]=min_data
        '''
        data_MASK = np.ones((yaxis, xaxis))
        for k in range(yaxis):
             for i in range(xaxis):
                  if data[k,i]==0.:
                        data_MASK[k,i] = 0.
        '''

        vmax = np.nanmax(data)
        vmin = np.nanmin(data)
           

        if numb_of_grid!=2:
            min_int = 10**(0.4*(m0-min_level)) * (pix2sec**2)
            im = grid[numb_of_grid].imshow(data,cmap=color, norm=LogNorm(vmin=min_int, vmax=vmax),interpolation='none', origin="lower", aspect='equal')
            if True:#numb_of_grid == 0:
                #if sigma is not None:
                #    data_MASK = ndimage.gaussian_filter(data_MASK, sigma=sigma, order=0)
                im = grid[numb_of_grid].imshow(data_MASK,cmap='Reds',interpolation='none', origin="lower", aspect='equal', alpha=.1)
        else:
          for k in range(yaxis):
             for i in range(xaxis):
                  if data_MASK[k,i]==0.:
                        data[k,i] = float('nan')
          if show_negative==False:
            data = np.fabs(data)
            cmap=discrete_cmap()
            vmin=0.001
            vmax=1
          else:
            cmap=discrete_cmap_neg() 
            vmin=-1
            vmax=1
          cmap.set_bad((0, 0, 0, 0))
          
          if view=='column':
            im = grid[numb_of_grid].imshow(data,cmap=cmap, vmin=vmin, vmax=vmax,interpolation='nearest', origin="lower", aspect='equal')
            cb = grid[numb_of_grid].cax.colorbar(im)
            for cax in grid.cbar_axes:
                cax.toggle_label(True)
                cax.axis[cax.orientation].set_label(' ')
                if view=='line':
                  cax.tick_params(labelsize=12)
                else:
                  cax.tick_params(labelsize=12)
                cax.set_xlim(vmin,vmax)
          if view=='line':                
            im = grid[numb_of_grid].imshow(data,cmap=cmap, vmin=vmin, vmax=vmax,interpolation='nearest', origin="lower", aspect='equal')
            cb = grid[numb_of_grid].cax.colorbar(im)
            for cax in grid.cbar_axes:
                cax.toggle_label(True)
                cax.axis[cax.orientation].set_label(' ')
                cax.tick_params(labelsize=12)
                cax.set_ylim(vmin,vmax)
          '''
          data_MASK = np.ones((yaxis, xaxis))
          for k in range(yaxis):
             for i in range(xaxis):
                  if data_MASK[k,i]==0.:
                        data_MASK[k,i] = float('nan')
          '''
          #cmap = plt.cm.Reds
          #cmap.set_bad((1, 0, 0, 1))
          #im = grid[numb_of_grid].imshow(data_MASK,cmap=cmap,interpolation='none', origin="lower", aspect='equal', alpha=.1)

        
        if region!=None:            
            r2 = pyregion.open(region).as_imagecoord(hdulist[layer].header)
            patch_list, artist_list = r2.get_mpl_patches_texts()
            for p in patch_list:
                grid[numb_of_grid].add_patch(p)
            for t in artist_list:
                grid[numb_of_grid].add_artist(t)            
        
        if numb_of_grid==0 and scale_bar!=None:    
            crea_scale_bar(grid[numb_of_grid],x0,x1,y0,y1,pix2sec,float(scale_bar),show_bar_label)

        if text!=None:
                  [GRID_NUMB,letter, hor_pos, vert_pos, text_size] = text
                  if GRID_NUMB==numb_of_grid:
                    grid[numb_of_grid].text(hor_pos, vert_pos, letter, fontsize=text_size, color='red',transform=grid[numb_of_grid].transAxes, horizontalalignment='left',verticalalignment='top',multialignment='left')#,multialignment.horizontalalignment='top')

        return x0,x1,y0,y1,data_MASK


def main(composed_model_file,min_level,pix2sec,m0,mask_file=None,
         borders='0,0,0,0',output_file='plot_2d.png',view='line',color='nipy_spectral_r',scale_bar=None,grid=None,show_bar_label='yes', region=None, text=None,show_negative=True,sigma=5.,rotate_PA=0., verbosity=True):
    no_fig = True
    if True:
        
        if grid==None:
          fig = plt.figure(0, (8, 8))
          fig.subplots_adjust(left=0.05, right=0.95)
          no_fig = False

          if view=='line':
              grid = AxesGrid(fig, 111,
                              nrows_ncols=(1, 3),
                              axes_pad=0.02,
                              label_mode="L",
                              share_all=True,
                              cbar_location="right",
                              cbar_mode='single',
                              cbar_size="5%",
                              cbar_pad="1%",
                              )

          elif view=='column':
              grid = AxesGrid(fig, 111,
                              nrows_ncols=(3, 1),
                              axes_pad=0.02,
                              label_mode="L",
                              share_all=True,
                              cbar_location="bottom",
                              cbar_mode='single',
                              cbar_size="8%",
                              cbar_pad="1%",
                              )

        for cax in grid.cbar_axes:
                cax.toggle_label(False)            


        if min_level==None:
            hdulist = pyfits.open(composed_model_file)
            data = hdulist[0].data
            mask_I = make_source_mask(data, nsigma=2, npixels=5, dilate_size=11)
            mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask_I)
            
            #mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5)
            #print std
            min_level = m0 - 2.5*log10(3.*std/(pix2sec**2))
            hdulist.close()
            if verbosity: print('We estimate the outer isophote at the 3-sigma level: %.2f' % (min_level))            
            
            
        #print min_level
        #exit()
        # Reference image:
        x0,x1,y0,y1,data_MASK = plot_indiv_frame(grid,composed_model_file,mask_file,'Data',m0,pix2sec,borders,min_level,color,scale_bar,view,show_bar_label,region,text,sigma=sigma,rotate_PA=rotate_PA)
        # Model image:
        x0,x1,y0,y1,data_MASK = plot_indiv_frame(grid,composed_model_file,mask_file,'Model',m0,pix2sec,borders,min_level,color,scale_bar,view,show_bar_label,region,text,data_MASK=data_MASK,sigma=sigma,rotate_PA=rotate_PA)
        # Residual image:
        x0,x1,y0,y1,data_MASK = plot_indiv_frame(grid,composed_model_file,mask_file,'Residual',m0,pix2sec,borders,min_level,color,scale_bar,view,show_bar_label,region=None,text=text,data_MASK=data_MASK,show_negative=show_negative,sigma=sigma,rotate_PA=rotate_PA)

        # This affects all axes because we set share_all = True.
        grid.axes_llc.set_xlim(x0,x1)
        grid.axes_llc.set_ylim(y0,y1)
        grid.axes_llc.set_xticklabels([])
        grid.axes_llc.set_yticklabels([])
        grid.axes_llc.get_xaxis().set_ticks([])        # To remove ticks
        grid.axes_llc.get_yaxis().set_ticks([])        # To remove ticks
        #print output_file
        #exit()
        #plt.show()
        #exit()

        
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
    parser.add_argument("min_level", help="Input minimum isophote level to highlight the structure in [mag/arcsec^2]")
    parser.add_argument("Scale", help="Input scale in [arcsec/pix]")
    parser.add_argument("m0", help="Input Zero Point in [mag/arcsec^2]")
    parser.add_argument("--mask_file", nargs='?', const=0., help="Name of the mask file",type=str,default=None)      
    parser.add_argument("--borders", nargs='?', const='0,0,0,0', help="Input the borders of the frame to be plotted: x1,y1,x2,y2",type=str,default='0,0,0,0')    
    parser.add_argument("--output_file", nargs='?', const=0., help="Name of the output file",type=str,default='plot_2d.png')      
    parser.add_argument("--view", nargs='?', const=0., help="Input the view of the plotted image: line or column",type=str,default='column')
    parser.add_argument("--color", nargs='?', const=0., help="Input the color for the plotted image: gray_r or nipy_spectral_r",type=str,default='nipy_spectral_r')
    parser.add_argument("--scale_bar", nargs='?', const=0., help="Create the scale bar: length [arcsec]",type=str,default=None)
    parser.add_argument("--show_bar_label", nargs='?', const=0., help="Show the scale bar label: yes",type=str,default='yes')
    parser.add_argument("--region", nargs='?', const=0., help="Region file",type=str,default=None)
    parser.add_argument("--PA", nargs='?', const=0., help="Totate image by PA (deg)",type=float,default=0.)
    
    args = parser.parse_args()
    
    composed_model_file = args.model_file
    min_level = float(args.min_level)
    pix2sec = float(args.Scale)
    m0 = float(args.m0)
    mask_file = args.mask_file
    borders = str(args.borders)
    output_file = args.output_file
    view = str(args.view)
    color = str(args.color)
    scale_bar = args.scale_bar
    show_bar_label = args.show_bar_label
    region = args.region
    rotate_PA = args.PA
    
    main(composed_model_file,min_level,pix2sec,m0,mask_file=mask_file,
         borders=borders,output_file=output_file,view=view,color=color,scale_bar=scale_bar,show_bar_label=show_bar_label,region=region,rotate_PA=rotate_PA)
