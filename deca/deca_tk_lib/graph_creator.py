#! /usr/bin/env python

from math import radians, sin

import Tkinter as Tk
import tkFileDialog
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from scipy.optimize import fmin
import numpy
from numpy import arange, linspace, abs, array, zeros_like, concatenate, log
from numpy import sum as npsum
import weakref

import time

import Image
from PIL import ImageTk

from copy import copy

import random as random_number
import sys
import math
import numpy as np
from scipy import stats
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.gridspec as gridspec
from matplotlib.ticker import NullFormatter
from mpl_toolkits.axes_grid1 import AxesGrid
from numpy import *
from pylab import *
import os
import shutil
import subprocess
from os.path import exists
import fileinput
import pyfits
import re
from astropy.stats import sigma_clipped_stats

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
sys.path.append(PATH_TO_PACKAGE.split('/DECA')[0]+'/Plot')
import plot_isomap
import plot_profile
import plot_2d_profile

fsize=10
xsize=8
ysize=8
  
def fig2data(fig):
    """
    @brief Convert a Matplotlib figure to a 4D numpy array with RGBA channels and return it
    @param fig a matplotlib figure
    @return a numpy 3D array of RGBA values
    """
    # draw the renderer
    fig.canvas.draw()
 
    # Get the RGBA buffer from the figure
    w,h = fig.canvas.get_width_height()
    buf = numpy.fromstring(fig.canvas.tostring_argb(), dtype=numpy.uint8)
    buf.shape = (w, h, 4)
 
    # canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
    buf = numpy.roll(buf, 3, axis=2)
    return buf
		  
def fig2img(fig):
    """
    @brief Convert a Matplotlib figure to a PIL Image in RGBA format and return it
    @param fig a matplotlib figure
    @return a Python Imaging Library ( PIL ) image
    """
    # put the figure pixmap into a numpy array
    buf = fig2data(fig)
    w, h, d = buf.shape
    return Image.fromstring("RGBA", (w, h), buf.tostring())

    
    
    
def Plot_isomap(galaxy_image,mask_image, mainGraph, canvas, m0, scale,xc=0.,yc=0., min_mag=20., max_mag=25.5, step=0.1, color=None, units=1, Scale=0.1): 
	if units==1:
	  units='pix'
	elif units==2:
	  units='arcsec'
	elif units==3:
	  units='kpc'
	  
	f = figure(0, figsize=(14, 10))
	mainGraph.subplots_adjust(hspace=0.25, wspace=0.22, left=0.06,right=0.99,top=0.95)
	ax1 = mainGraph.add_subplot(111)
	ax1.clear(); ax1.cla(); plt.clf();f.clf(); plt.close(f)
	#plot_isomap.main(galaxy_image,m0,scale,min_mag,max_mag,step,simple_map_color=color,AX=ax1)
	plot_isomap.main(galaxy_image,m0,scale,min_mag,max_mag,step,
	 xc=xc,yc=yc,x0=0.,x1=0.,y0=0.,y1=0.,mask_image=mask_image,
	 AX=ax1, output_image = None, simple_map_color=color, name=None, scalebar=None, units=units, Scale=Scale)

        imageWithGraph = fig2img(mainGraph)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)

	return ImageTk.PhotoImage(imageWithGraph)
        

def Plot_Prof1D(galaxy_image, mask_image, mainGraph, canvas, m0, scale, type_prof='rad', xc=0., yc=0., PA=0., Rmin=0., Rmax=0., step=1., zmin=0., zmax=0.,geom_units='arcsec',SB_units='mag/arcsec2',Scale=0.1):
	f = figure(1, figsize=(14, 10))
	mainGraph.subplots_adjust(hspace=0.25, wspace=0.22, left=0.06,right=0.99,top=0.95)
	ax1 = mainGraph.add_subplot(111)
	ax1.clear(); ax1.cla(); plt.clf();f.clf(); plt.close(f)
	plot_profile.main(galaxy_image,m0,scale,mask_image=mask_image,profile = type_prof,xc=xc,yc=yc,PA=PA,Rmin=Rmin,Rmax=Rmax,
		   step=step,zmin=zmin,zmax=zmax,AX=ax1,geom_units=geom_units,SB_units=SB_units,Scale=Scale)  

        imageWithGraph = fig2img(mainGraph)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)

	return ImageTk.PhotoImage(imageWithGraph)     

'''
def Plot_Prof2D(galaxy_image, mainGraph, canvas, m0, scale, type_prof='rad', xc=0., yc=0., PA=0., Rmin=0., Rmax=0., step=1.):
	f = figure(0, figsize=(14, 10))
	mainGraph.subplots_adjust(hspace=0.25, wspace=0.22, left=0.06,right=0.99,top=0.95)
	ax1 = mainGraph.add_subplot(111)
	ax1.clear(); ax1.cla(); plt.clf();f.clf(); plt.close(f)
	plot_profile.main(galaxy_image,m0,scale,profile = type_prof,xc=xc,yc=yc,PA=PA,Rmin=Rmin,Rmax=Rmax,outp_format='eps',AX=ax1)


        imageWithGraph = fig2img(mainGraph)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)

	return ImageTk.PhotoImage(imageWithGraph)     

'''      
def Plot_Prof2D(composed_model_file, mainGraph2D, canvas2D, min_level,pix2sec,m0,mask_file,view,color,scale_bar):
	ff = figure(33, figsize=(8, 8))
	mainGraph2D.subplots_adjust(hspace=0.25, wspace=0.22, left=0.06,right=0.99,top=0.95)
	#ax1 = mainGraph.add_subplot(111, aspect='equal')
	plt.clf();ff.clf()

	if view==1:
	      View='line'
	      grid = AxesGrid(ff, 111,
			      nrows_ncols=(1, 3),
			      axes_pad=0.02,
			      label_mode="L",
			      share_all=True,
			      cbar_location="right",
			      cbar_mode='single',
			      cbar_size="5%",
			      cbar_pad="1%",
			      )
	       	      
	elif view==2:
	       View='column'
	       grid = AxesGrid(ff, 111,
			      nrows_ncols=(3, 1),
			      axes_pad=0.02,
			      label_mode="L",
			      share_all=True,
			      cbar_location="bottom",
			      cbar_mode='single',
			      cbar_size="8%",
			      cbar_pad="1%",
			      )



	plot_2d_profile.main(composed_model_file,min_level,pix2sec,
			      m0,mask_file=mask_file,view=View,color=color,scale_bar=scale_bar,grid=grid)
        
        imageWithGraph2D = fig2img(mainGraph2D)
        canvas2D.show()
        canvas2D.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)

	return ImageTk.PhotoImage(imageWithGraph2D)   



def Plot_Alignment(galaxy_image, mainGraph, canvas, profile, xc, yc, Rmax, PA=None, Rmin=None, zmin=None, zmax=None): 
	ff = figure(3, figsize=(14, 10))
	mainGraph.subplots_adjust(hspace=0.25, wspace=0.22, left=0.06,right=0.99,top=0.95)
	ax1 = mainGraph.add_subplot(111, aspect='equal')
	ax1.clear(); ax1.cla(); plt.clf();ff.clf(); plt.close(ff)
	
	hdulist = pyfits.open(galaxy_image)
	data = hdulist[0].data
	

	vmin = np.min(fabs(data))
	vmax = data[int(yc),int(xc)]
	im = ax1.imshow(data,cmap='nipy_spectral_r', vmin=vmin, vmax=vmax, interpolation='none', origin="lower", aspect='equal')
	
	# Show center
	ax1.plot([xc],[yc],'x',markersize=10,color='black')
	
	if profile=='cut':
	  # Draw line
	  endy1 = yc + Rmax * math.sin(math.radians(PA))
	  endx1 = xc + Rmax * math.cos(math.radians(PA))
	  endy2 = yc + Rmax * math.sin(math.radians(PA+180.))
	  endx2 = xc + Rmax * math.cos(math.radians(PA+180.))
	  ax1.plot([endx1,endx2], [endy1,endy2], color='black')
	
	if profile=='azim' or profile=='rad':
	  circ = patches.Circle((xc,yc), radius=Rmax, color='black', hatch='/', fill=False)
	  ax1.add_patch(circ)

	if profile=='summed':
	  if zmin==0.:
	    rect = patches.Rectangle((xc-Rmax,yc-zmax), 2.*Rmax, 2.*zmax, hatch='/', fill=False)
	    ax1.add_patch(rect)
	  else:
	    rect1 = patches.Rectangle((xc-Rmax,yc+zmin), 2.*Rmax,zmax-zmin, hatch='/', fill=False)
	    ax1.add_patch(rect1)    
	    rect2 = patches.Rectangle((xc-Rmax,yc-zmax), 2.*Rmax,zmax-zmin, hatch='/', fill=False)
	    ax1.add_patch(rect2)
	    
	if profile=='vert':
	  if Rmin==0.:
	    rect = patches.Rectangle((xc-Rmax,yc-zmax), 2.*Rmax, 2.*zmax, hatch='/', fill=False)
	    ax1.add_patch(rect)
	  else:
	    rect1 = patches.Rectangle((xc-Rmax,yc-zmax), Rmax-Rmin, 2.*zmax, hatch='/', fill=False)
	    ax1.add_patch(rect1)    
	    rect2 = patches.Rectangle((xc+Rmin,yc-zmax), Rmax-Rmin, 2.*zmax, hatch='/', fill=False)
	    ax1.add_patch(rect2)    


        imageWithGraph = fig2img(mainGraph)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)

	return ImageTk.PhotoImage(imageWithGraph) 

def Plot_IRAF(input_files,labels, mainGraph, canvas,m0,scale,Scale=0.1,geom_units='arcsec',SB_units='mag/arcsec2'):
	f = figure(44, figsize=(14, 10))
	mainGraph.subplots_adjust(hspace=0.25, wspace=0.22, left=0.06,right=0.99,top=0.95)

	ax1 = mainGraph.add_subplot(221)
	ax2 = mainGraph.add_subplot(222)
	ax3 = mainGraph.add_subplot(223)
	ax4 = mainGraph.add_subplot(224)


	ax1.clear(); ax2.clear(); ax3.clear(); ax4.clear()
	ax1.cla(); ax2.cla(); ax3.cla(); ax4.cla(); f.clf(); plt.close(f)

	if geom_units=='arcsec':
	  ax1.set_xlabel(r'r (arcsec)', fontsize=fsize)
	  ax2.set_xlabel(r'r (arcsec)', fontsize=fsize)
	  ax3.set_xlabel(r'r (arcsec)', fontsize=fsize)
	  ax4.set_xlabel(r'r (arcsec)', fontsize=fsize)
	elif geom_units=='pix':
	  ax1.set_xlabel(r'r (pix)', fontsize=fsize)
	  ax2.set_xlabel(r'r (pix)', fontsize=fsize)
	  ax3.set_xlabel(r'r (pix)', fontsize=fsize)
	  ax4.set_xlabel(r'r (pix)', fontsize=fsize)
	elif geom_units=='kpc':
	  ax1.set_xlabel(r'r (kpc)', fontsize=fsize)
	  ax2.set_xlabel(r'r (kpc)', fontsize=fsize)
	  ax3.set_xlabel(r'r (kpc)', fontsize=fsize)
	  ax4.set_xlabel(r'r (kpc)', fontsize=fsize)

	ax1.set_ylabel(r'PA (deg)', fontsize=fsize)
	ax2.set_ylabel(r'$\epsilon$ ', fontsize=fsize+3)
	ax3.set_ylabel(r'B$_4$ ', fontsize=fsize)
	if SB_units=='mag/arcsec2':
	  ax4.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
	elif SB_units=='ADU/pix2':
	  ax4.set_ylabel(r' Intensity (ADU) ', fontsize=fsize)

	ax1.yaxis.set_tick_params(labelsize=ysize)
	ax1.xaxis.set_tick_params(labelsize=xsize)
	ax2.yaxis.set_tick_params(labelsize=ysize)
	ax2.xaxis.set_tick_params(labelsize=xsize)
	ax3.yaxis.set_tick_params(labelsize=ysize)
	ax3.xaxis.set_tick_params(labelsize=xsize)
	ax4.yaxis.set_tick_params(labelsize=ysize)
	ax4.xaxis.set_tick_params(labelsize=xsize)


  
	for kkk in range(len(input_files)):
	      file = input_files[kkk]
	      label = labels[kkk]

	      sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = loadtxt(file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 6, dtype='str')

	      for k in range(len(sma)):
			if sma[k]=='INDEF': sma[k]=0
			if inten[k]=='INDEF': inten[k]=0
			if inten_err[k]=='INDEF': inten_err[k]=0
			if ell[k]=='INDEF': ell[k]=0
			if errell[k]=='INDEF': errell[k]=0
			if PA[k]=='INDEF': PA[k]=0
			if errPA[k]=='INDEF': errPA[k]=0
			if x0[k]=='INDEF': x0[k]=0
			if y0[k]=='INDEF': y0[k]=0
			if B4[k]=='INDEF': B4[k]=0
			if errB4[k]=='INDEF': errB4[k]=0
	      sma = np.array(sma,dtype='float')
	      inten = np.array(inten,dtype='float')
	      inten_err = np.array(inten_err,dtype='float')
	      ell = np.array(ell,dtype='float')
	      errell = np.array(errell,dtype='float')
	      PA = np.array(PA,dtype='float')
	      errPA = np.array(errPA,dtype='float')
	      x0 = np.array(x0,dtype='float')
	      y0 = np.array(y0,dtype='float')
	      B4 = np.array(B4,dtype='float')
	      errB4 = np.array(errB4,dtype='float')
	      
	      if SB_units=='mag/arcsec2':
		mag = m0-2.5*log10(inten)+5.*log10(scale)
		mag_err = fabs((2.5/log(10.0)) * inten_err/inten)
	      elif SB_units=='ADU/pix2':
		mag = inten
		mag_err = inten_err
		
	      for kk in range(len(PA)):
		if PA[kk]<0.:
		  PA[kk] = PA[kk] + 180.

	      if geom_units=='arcsec':
		Radius = sma*scale
	      elif geom_units=='pix':
		Radius = sma
	      if geom_units=='kpc':
		Radius = sma*scale*Scale



      
	      ax1.errorbar(Radius,  PA, yerr=errPA, fmt='o',color='black', markersize=3)
	      ax2.errorbar(Radius,  ell, yerr=errell, fmt='o',color='black', markersize=3)
	      ax3.errorbar(Radius,  B4, yerr=errB4, fmt='o',color='black', markersize=3)
	      ax4.errorbar(Radius, mag, mag_err, color="white", marker='o',markersize=6, label=label)
	      #print kk
	      if kkk==0:
		  mean_mag, median_mag, std_mag = sigma_clipped_stats(mag, sigma=3.0, iters=5)
		  mean_B4, median_B4, std_B4 = sigma_clipped_stats(B4, sigma=3.0, iters=5)

		  ax1.set_ylim(np.min(fabs(PA)-errPA)*0.99,np.max(fabs(PA)+errPA)*1.01)
		  ax2.set_ylim(np.min(fabs(ell)-errell)*0.99,np.max(fabs(ell)+errell)*1.01)
		  ax3.set_ylim(median_B4-5*std_B4,median_B4+5*std_B4)
		  if SB_units=='mag/arcsec2':
		    ax4.set_ylim(max(mag)+std_mag,min(mag)-std_mag)
		  elif SB_units=='ADU/pix2':
		    ax4.set_ylim(min(mag)-std_mag,max(mag)+std_mag)
	      
	ax4.legend(loc=0, numpoints=1, prop={'size':fsize},markerscale=1)

	    
        imageWithGraph = fig2img(mainGraph)
        canvas.show()
        canvas.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)
	return ImageTk.PhotoImage(imageWithGraph)
          
