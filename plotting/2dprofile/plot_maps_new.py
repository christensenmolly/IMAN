#!/usr/bin/python
# -*- coding:  cp1251 -*-
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

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
warnings.filterwarnings("ignore")

#***************************
scale_bar_length = 30
offset_x_factor = 0.98
offset_y_factor = 0.1
format_of_output_file = 'png'
#***************************


def plot_profiles(AX,data_ref,data_ski,pix2sec,m0,profile,Label,c,MARKER,xc,Rmin,Rmax):
      ySize, xSize = data_ref.shape  
      #print Label,np.sum(data_ref),np.sum(data_ski)
      if profile!='vert':
		I_ref = []; I_ski = []
		r = []
		'''
		for x in range(0,xSize,1):
		  I_ref.append(sum(data_ref[:,x]))
		  I_ski.append(sum(data_ski[:,x]))
		  r.append(x-xSize/2.)
		'''
		I_ref = data_ref[int(ySize/2.),:]
		I_ski = data_ski[int(ySize/2.),:]
		for x in range(0,len(I_ref),1):
		  r.append(x-xSize/2.)
		
		
		r = np.array(r)
		I_ref = np.array(I_ref)
		I_ski = np.array(I_ski)
      else:
              if xc==0. and Rmin==0. and Rmax==0.:
		I_ref = []; I_ski = []
		r = []
		for y in range(0,ySize,1):
		  II_ref = []; II_ski = []
		  for x in range(0,xSize,1):
		    if data_ref[y,x]>0.:
		      II_ref.append(data_ref[y,x])
		      II_ski.append(data_ski[y,x])
		  I_ref.append(sum(II_ref))
		  I_ski.append(sum(II_ski))
		  r.append(y-ySize/2.)
		r = np.array(r)
		I_ref = np.array(I_ref)
		I_ski = np.array(I_ski)   
              else:
		I_ref = []; I_ski = []
		r = []
		for y in range(0,ySize,1):
		  II_ref = []; I_ski = []
		  ranges = range(int(xc-Rmax),int(xc-Rmin),1)+range(int(xc+Rmin),int(xc+Rmax),1)
		  for x in ranges:
		    if data_ref[y,x]>0.:
		      II_ref.append(data_ref[y,x])
		      II_ski.append(data_ski[y,x])
		  I_ref.append(sum(II_ref))
		  I_ski.append(sum(II_ski))
		  r.append(y-ySize/2.)
		r = np.array(r)
		I_ref = np.array(I_ref)
		I_ski = np.array(I_ski)
	    
  
      m0 = 2.5*log10(np.max(I_ski)) - 5.*log10(pix2sec) - m0
      if profile=='summed' or profile=='vert':
		AX.plot(r*pix2sec, m0 - 2.5*log10(I_ref)+ 5.*log10(pix2sec),MARKER,color=c,markeredgecolor=c,markersize=5,label = Label)
		min_I = m0 - 2.5*log10(I_ref)+ 5.*log10(pix2sec)
		min_mag = np.min(min_I[~np.isnan(min_I)])
		max_mag = np.max(min_I[~np.isnan(min_I) & ~np.isinf(min_I)])
		
		AX.plot(r*pix2sec, m0 - 2.5*log10(I_ski)+ 5.*log10(pix2sec),'-',color='black',markeredgecolor='black',lw=2)
     
      #AX.legend(loc=2, borderaxespad=0.,fontsize=10,numpoints=1)
      
      min_r = -max(fabs(r*pix2sec))
      max_r = max(fabs(r*pix2sec))
      return min_r,max_r,min_mag,max_mag


















'''
def demo_grid_with_each_cbar(fig):
    """
    A grid of 2x2 images. Each image has its own colorbar.
    """

    grid = AxesGrid(fig, 143,  # similar to subplot(143)
                    nrows_ncols=(2, 2),
                    axes_pad=0.1,
                    label_mode="1",
                    share_all=True,
                    cbar_location="top",
                    cbar_mode="each",
                    cbar_size="7%",
                    cbar_pad="2%",
                    )
    Z, extent = get_demo_image()
    for i in range(4):
        im = grid[i].imshow(Z, extent=extent, interpolation="nearest")
        grid.cbar_axes[i].colorbar(im)

    for cax in grid.cbar_axes:
        cax.toggle_label(False)

    # This affects all axes because we set share_all = True.
    grid.axes_llc.set_xticks([-2, 0, 2])
    grid.axes_llc.set_yticks([-2, 0, 2])
'''

fsize = 2
def find_ski_image(name):
  return glob.glob('./plot_simulations/reduced_ski_images_masked/'+name+'_*.fits')[0]
  
def find_ref_image(name,band):
  return glob.glob('./plot_simulations/reduced_images_masked/'+name+'_'+band+'_*.fits')[0]

def crea_scale_bar(ax,nx,ny,pix2sec):
  xc = offset_x_factor*(nx*pix2sec-scale_bar_length)
  yc = offset_y_factor*ny*pix2sec
  ax.errorbar(xc, yc, xerr=scale_bar_length,color='lime',capsize=2)
  #ax.annotate('30\"', (xc*0.98,yc*1.08), color='lime',fontsize=5)

def split_seq(seq, p):
    newseq = []
    n = len(seq) / p    # min items per subsequence
    r = len(seq) % p    # remaindered items
    b,e = 0, n + min(1, r)  # first split
    for i in range(p):
        newseq.append(seq[b:e])
        r = max(0, r-1)  # use up remainders
        b,e = e, e + n + min(1, r)  # min(1,r) is always 0 or 1

    return newseq  
'''
def find_min_and_max(ref_image,ski_image):
    ySize, xSize = ref_image.shape
    for k in range(ySize):
      for i in range(xSize):
	if ref_image[k,i]>0.0000001:
	  su = max_v + 
'''    

def main(NUMBERS,NAMES,PIX2SEC,fluxes_file,sed_file,split_plots=1):
    Numbers = split_seq(NUMBERS, split_plots)
    Names = split_seq(NAMES, split_plots)
    Pix2sec = split_seq(PIX2SEC, split_plots)

    # STEP Additional calibration:
    SURVEY,BAND,WAVELENGTH,FLUX,FLUX_ERR = loadtxt(fluxes_file, usecols=[0,1,2,3,4],dtype=str, unpack=True,skiprows=1,delimiter='\t')
    WAVELENGTH = np.array(WAVELENGTH,float)
    FLUX = np.array(FLUX,float)
    FLUX_ERR = np.array(FLUX_ERR,float)
    
    WAVELENGTH_SKIRT,FLUX_SKIRT = loadtxt(sed_file, usecols=[0,1],dtype=float, unpack=True,skiprows=7,delimiter=' ')


    #idx = find_nearest(array,value)
    #print NUMBERS,NAMES,PIX2SEC
    #exit()
    for k in range(len(PIX2SEC)):
	    if PIX2SEC[k]==1.:
		  #print NAMES[k]
		  ski_image = './plot_simulations/reduced_ski_images_masked/'+str(NUMBERS[k])+'_model_conv_rebin_red.fits'
		  hdulist_ski = pyfits.open(ski_image)
		  data_ski = hdulist_ski[0].data 
		  xaxis = data_ski.shape[1]
		  yaxis = data_ski.shape[0]
		  XAXIS = xaxis
		  YAXIS = yaxis
		  break
    #exit()
    for k in range(len(Numbers)):
	numbers = Numbers[k]
	names = Names[k]
	pix2sec = Pix2sec[k]
      
	fig = plt.figure(k, (4, 20))
	fig.subplots_adjust(left=0.05, right=0.95) 

	grid = AxesGrid(fig, 111,
			nrows_ncols=(len(numbers), 2),
			axes_pad=0.02,
			label_mode="L",
			share_all=True,
			cbar_location="bottom",
			cbar_mode=None,
			cbar_size="7%",
			cbar_pad="2%",
			) # cbar_mode="single"
	count = 0
	for i in range(len(numbers)):
	    #ski_image = './reduced_ski_images/' + str(numbers[i]) + '_model_conv_rebin_red.fits'
	    ski_image = find_ski_image(str(numbers[i]))

	    #ref_image = './reduced_images/' + names[i].split('_')[0]+'_'+names[i].split('_')[1]+'_a_r_ext_wcs_clean_rebin_red.fits'
	    ref_image = find_ref_image(names[i].split('_')[0],names[i].split('_')[1])

	    hdulist_ski = pyfits.open(ski_image)
	    data_ski = hdulist_ski[0].data 
	    xaxis = data_ski.shape[1]
	    yaxis = data_ski.shape[0]

	    '''
	    if pix2sec[i]==1.:
	      XAXIS = xaxis
	      YAXIS = yaxis
	    '''

	    hdulist_ref = pyfits.open(ref_image)
	    data_ref = hdulist_ref[0].data

	    xaxis = data_ref.shape[1]
	    yaxis = data_ref.shape[0]


	    # Right normalization
	    #data_ref = data_ref * FLUX[int(numbers[i])-1]
	    data_ref = data_ref / np.sum(data_ref)
	    #data_ski = data_ski * FLUX_SKIRT[int(numbers[i])-1]
	    data_ski = data_ski / np.sum(data_ski)

	    #print FLUX[int(numbers[i])-1],FLUX_SKIRT[int(numbers[i])-1]
	    #data_ref = data_ref * FLUX[int(numbers[i])-1] / np.sum(data_ref)
	    #data_ski = data_ski * FLUX_SKIRT[int(numbers[i])-1] / np.sum(data_ski)
	    

	    for l in range(yaxis):
		for m in range(xaxis):
		    if data_ref[l,m]<=0.000000001:
			data_ref[l,m]=0.
			data_ski[l,m]=0.



	    #mod_data = fabs(data_ref)
            #min_data = np.min(mod_data[np.nonzero(fabs(mod_data))])
            #for l in range(yaxis):
            #    for m in range(xaxis):
            #        if data_ref[l,m]<=0.0000001:
	    #		data_ref[l,m]=min_data
	    #		data_ski[l,m]=min_data

			
	    #vmax = np.max(data_ref)
	    #vmin = np.min(data_ski)
	    vmax = np.max(data_ref)
	    #vmin = np.min(data_ski)
	    vmin = np.max([np.max(data_ski[:,0]),np.max(data_ski[:,xaxis-1])])
	    
	    
	    
	    #print i,vmin,vmax
	    #print names[i],vmin,vmax
	    im = grid[count].imshow(log10(data_ref),cmap='nipy_spectral_r',vmin=-7, vmax=log10(vmax), extent=(0,xaxis*pix2sec[i],0,yaxis*pix2sec[i]),interpolation='none', origin="lower")	#hot or gray_r
	    #im = grid[count].imshow(data_ref,cmap='jet',norm=LogNorm(vmin=vmin, vmax=vmax), extent=(0,xaxis*pix2sec[i],0,yaxis*pix2sec[i]),interpolation='none')
	    #grid.cbar_axes[count].colorbar(im)
	    if names[i].split('_')[1]=='w1' and SURVEY[int(numbers[i])-1]=='Spitzer':
	      text_label = names[i].split('_')[0] + '[3.6]'
	    else:
	      text_label = names[i].split('_')[0] + '['+names[i].split('_')[1]+'] '
	    grid[count].text(0.03, 0.88,text_label, color='black',transform=grid[count].transAxes,fontsize=fsize+2, fontweight='bold', va='top', bbox=dict(facecolor='white', edgecolor='black',linewidth=0.5))
	    
	    #if i==0:
	    #  crea_scale_bar(grid[count],xaxis,yaxis,pix2sec[i])

	    #vmax = np.max(data_ski)
	    #vmin = np.min(data_ski)
	    #print vmin,vmax, '\n'
	    #print i, np.sum(data_ref),np.sum(data_ski)
	    im = grid[count+1].imshow(log10(data_ski),cmap='nipy_spectral_r',vmin=-7, vmax=log10(vmax), extent=(0,xaxis*pix2sec[i],0,yaxis*pix2sec[i]),interpolation='none', origin="lower")	#hot or gray_r
	    #im = grid[count+1].imshow(data_ski,cmap='jet',norm=LogNorm(vmin=vmin, vmax=vmax), extent=(0,xaxis*pix2sec[i],0,yaxis*pix2sec[i]),interpolation='none')
	    #grid.cbar_axes[count+1].colorbar(im)
	    count = count + 2
	    #im.set_xticklabels([])
	    #im.set_yticklabels([])


	for cax in grid.cbar_axes:
	    cax.toggle_label(False)

	# This affects all axes because we set share_all = True.
	grid.axes_llc.set_xlim(0,XAXIS)
	grid.axes_llc.set_ylim(0,YAXIS)
	grid.axes_llc.set_xticklabels([])
	grid.axes_llc.set_yticklabels([])
	grid.axes_llc.get_xaxis().set_ticks([])	# To remove ticks
	grid.axes_llc.get_yaxis().set_ticks([])	# To remove ticks

	plt.draw()  
	plt.savefig('map_%i_of_%i.eps' % (k+1,len(Numbers)), bbox_inches='tight', pad_inches=0.01, dpi = 300)
	plt.clf()
	plt.close()   


    
    
    
    
    
    
    
    
    
    
    
    
    color=iter(cm.rainbow(np.linspace(0,1,len(WAVELENGTH))))
    
    Labels = list(SURVEY)
    filled_markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd','o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd']
    Markers = filled_markers[:len(Labels)]
    

    used_labels = []    
    
    
    

    #### 1D profiles
    count = 0
    #print Labels
    #exit()
    for k in range(len(Numbers)):
	numbers = Numbers[k]
	names = Names[k]
	pix2sec = Pix2sec[k]
      
	
	
	fig1 = plt.figure(k+1,figsize=(5, 7))
	AX1 = fig1.add_subplot(111)

	fig2 = plt.figure(k+2,figsize=(5, 7))
	AX2 = fig2.add_subplot(111)


	AX1.set_xlabel(r' z (arcsec) ', fontsize=15)
	AX2.set_xlabel(r' r (arcsec) ', fontsize=15)
	AX1.set_ylabel(r' $\mu$ (mag arcsec$^{-2}$) ', fontsize=15)
	AX2.set_ylabel(r' $\mu$ (mag arcsec$^{-2}$) ', fontsize=15)
	  
	  
	
	#AX2.gca().invert_yaxis()

	output_file1 = 'vertical_profiles_%i_of_%i.pdf' % (k+1,len(Numbers))
	output_file2 = 'horizontal_profiles_%i_of_%i.pdf' % (k+1,len(Numbers))

	MIN_z=[]; MAX_z=[]; MIN_r=[]; MAX_r=[]; MIN_mag_z=[]; MAX_mag_z=[]; MIN_mag_r=[];MAX_mag_r=[]  

	m0 = 0.

	for i in range(len(numbers)):
	    c=next(color)
	    MARKER = Markers[count]
	    
	    
	    ski_image = find_ski_image(str(numbers[i]))
	    ref_image = find_ref_image(names[i].split('_')[0],names[i].split('_')[1])
	    hdulist_ski = pyfits.open(ski_image)
	    data_ski = hdulist_ski[0].data 
	    #data_ski = data_ski / np.sum(data_ski)

	    xaxis = data_ski.shape[1]
	    yaxis = data_ski.shape[0]

	    
	    if pix2sec[i]==1.:
	      XAXIS = xaxis
	      YAXIS = yaxis


	    hdulist_ref = pyfits.open(ref_image)
	    data_ref = hdulist_ref[0].data
	    #data_ref = data_ref / np.sum(data_ref)
	    xaxis = data_ref.shape[1]
	    yaxis = data_ref.shape[0]

	    data_ref = data_ref * FLUX[int(numbers[i])-1]
	    data_ref = data_ref / np.sum(data_ref)
	    data_ski = data_ski * FLUX_SKIRT[int(numbers[i])-1]
	    data_ski = data_ski / np.sum(data_ski)
	    #print WAVELENGTH[int(numbers[i])-1],np.sum(data_ref),np.sum(data_ski)
	    

	    for l in range(yaxis):
		for m in range(xaxis):
		    if data_ref[l,m]<=0.000000001:
			data_ref[l,m]=0.
			data_ski[l,m]=0.



		    

	    if names[i].split('_')[1]=='w1':
	      text_label = names[i].split('_')[0] + '[3.6]'
	    else:
	      text_label = names[i].split('_')[0] + '['+names[i].split('_')[1]+']'

	    min_z,max_z,min_mag_z,max_mag_z = plot_profiles(AX1,data_ref,data_ski,pix2sec[i],m0,'vert',text_label,c,MARKER,0.,0.,0.)#Rmin=int(50./pix2sec[i]),Rmax=int(150./pix2sec[i]))
	    min_r,max_r,min_mag_r,max_mag_r = plot_profiles(AX2,data_ref,data_ski,pix2sec[i],m0,'summed',text_label,c,MARKER,0.,0.,0.)#Rmin=int(50./pix2sec[i]),Rmax=int(150./pix2sec[i]))
	    MIN_z.append(min_z)
	    MAX_z.append(max_z)
	    MIN_r.append(min_r)
	    MAX_r.append(max_r)
	    
	    MIN_mag_z.append(min_mag_z)
	    MAX_mag_z.append(max_mag_z)
	    MIN_mag_r.append(min_mag_r)
	    MAX_mag_r.append(max_mag_r)	    
  
	    count = count + 1
	    m0 = m0 -2.



	AX1.set_xlim([min(MIN_z),max(MAX_z)])
	AX2.set_xlim([min(MIN_r),max(MAX_r)])

	AX1.set_ylim([min(MIN_mag_z),max(MAX_mag_z)])
	AX2.set_ylim([min(MIN_mag_r),max(MAX_mag_r)])

	AX1.invert_yaxis()
	AX2.invert_yaxis()
	fig1.savefig(output_file1, bbox_inches='tight', pad_inches=0.01, dpi = 300)
	fig2.savefig(output_file2, bbox_inches='tight', pad_inches=0.01, dpi = 300)
	plt.clf()
	plt.close()   

    
    
  
  
'''  
file_with_ref_data = 'info.dat'  
  
numbers,wavelengths,names,kernels,FWHM,pix2sec,input_units = loadtxt(file_with_ref_data, usecols=[0,1,2,3,4,5,6],dtype=str, unpack=True,skiprows=1,delimiter='\t')
main(numbers,names,np.array(pix2sec,float))
'''