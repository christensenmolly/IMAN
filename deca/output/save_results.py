#!/usr/bin/env python
# Module to save output files and create plots using the results of decomposition
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **        © Astronomical Observatory, Ghent University         **
# *****************************************************************

# Import standard modules
import pylab
import sys
import os
import shutil
import math
import numpy as np
import scipy as sp
from numpy import *
from pylab import *
import subprocess
import argparse
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import glob

# Import astro packages
import pyfits
from photutils import datasets
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import EllipticalAperture
from astropy.table import hstack


tmp_out = sys.stdout

# Import DECA modules
import deca_setup

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0].split('DECA')[0]

sys.path.append(PATH_TO_PACKAGE+'Decomposition')
sys.path.append(PATH_TO_PACKAGE+'Plot')
import make_model_ima_galfit
import make_model_ima_imfit
import plot_profile
import plot_2d_profile


# -----------------------------------------------------------------
# Function to copy the images from the certain directory to the new one
def copy_images(IMAGES,directory):
    new_images = []
    for k in range(len(IMAGES)):
      #print IMAGES
      new_images.append(IMAGES[k].split('/')[-1])
      if IMAGES[k]!='none' and os.path.exists(IMAGES[k]):
	shutil.copy(IMAGES[k],directory+IMAGES[k].split('/')[-1])
    return new_images
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO CREATE TEX FILE 
def create_tex(self):
      [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,Sky_level,SkySubtr,sampling,convolution_box,FILTER] = observation_info
      [xc,yc,name] = object_info
      [number,fitting_proc,mode,model,region_to_fit] = add_info
      [input_image,sigma_image,psf_image,mask_image] = new_images
      [control,delay] = keys
    
      #input_files,observation_info,object_info,settings = self.inp_f()
      #[input_image,sigma_image,psf_image,mask_image,sampling] = input_files
      #[nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,Filter] = observation_info
      #[RA,DEC,galaxy_name,Dist,Ext,Kcorr] = object_info
      
      [ChooseCode,ChooseGeom,ChooseSB,ChooseLum] = settings
      Scale = distance * 1000. / 206265.
      
      # Define results file and composed image
      composed_image = self.SaveModel.get()+'/'+'composed_model.fits'
      if ChooseCode=='IMFIT':
	results_file = self.SaveModel.get()+'/'+'imfit.01'
      elif ChooseCode=='GALFIT':
	results_file = self.SaveModel.get()+'/'+'galfit.01'
      
      hdu = pyfits.open(composed_image)
      number_of_layers = len(hdu)
      fluxes = []
      for k in range(number_of_layers):
	  if k==1:
	    Flux = np.sum(hdu[k].data)
	  if k>=4:
	    fluxes.append(np.sum(hdu[k].data))
      luminosities = []
      for k in range(len(fluxes)):
	luminosities.append(fluxes[k]/Flux)
      luminosities.append(Flux)
      
      # Create pictures
      pictures = []
      pictures.append(plot_2d_profile.main(composed_image,25.5,scale,m0,mask_file=mask_image))
      print pictures
      tex_creator.main(galaxy_name,results_file,ChooseCode,luminosities,pictures,m0,scale,Scale,Filter,Ext,Kcorr,self.ChooseGeomOut.get(),self.ChooseLumOut.get(),self.ChooseSBOut.get())
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# MAIN FUNCTION
def main(images,pix2sec,m0,new_dir=None):
 #input_image,mask_image,sigma_image,psf_image = new_images
 [input_image,sigma_image,psf_image,mask_image] = images

 try:
  if deca_setup.code=='GALFIT':
	#1. Create new folder for this model
	if new_dir==None:
            arr = glob.glob('galfit_*')
            dir_numbers = []
            for ar in arr:
                dir_numbers.append(int(ar.split('_')[-1]))
            if dir_numbers==[]:
                new_dir = './galfit_1'
            else:
                new_dir = './galfit_' + str(max(dir_numbers)+1)           
	os.makedirs(new_dir)

	#2. Create log file if it hasn't been created
	if not os.path.exists("galfit.log"):
	  shutil.move("fit.log","galfit.log")
	else:
	  with open("fit.log") as f:
	    with open("galfit.log", "a") as f1:
		for line in f:
			f1.write(line) 

	#3. Create composed model
	make_model_ima_galfit.main(input_image, 'galfit.01', subtract_sky=True, galfitPath=deca_setup.galfitPath)
        #print 'HEREE'

  elif deca_setup.code=='IMFIT':
	#1. Create new folder for this model
	if new_dir==None:
            arr = glob.glob('imfit_*')
            dir_numbers = []
            for ar in arr:
                dir_numbers.append(int(ar.split('_')[-1]))
            if dir_numbers==[]:
                new_dir = './imfit_1'
            else:
                new_dir = './imfit_' + str(max(dir_numbers)+1)
	os.makedirs(new_dir)

	#2. Create log file if it hasn't been created
	if not os.path.exists("imfit.log"):
	  shutil.copy("imfit.01","imfit.log")
	else:
	  with open("imfit.01") as f:
	    with open("imfit.log", "a") as f1:
		for line in f:
			f1.write(line) 
	
	#3. Create composed model
	make_model_ima_imfit.main(input_image, 'imfit.01', psf_image, imfitPath=deca_setup.imfitPath)
  elif deca_setup.code=='GALFITM':
	#1. Create new folder for this model
	if new_dir==None:
            arr = glob.glob('galfitm_*')
            dir_numbers = []
            for ar in arr:
                dir_numbers.append(int(ar.split('_')[-1]))
            if dir_numbers==[]:
                new_dir = './galfitm_1'
            else:
                new_dir = './galfitm_' + str(max(dir_numbers)+1)           
	os.makedirs(new_dir)

	#2. Create log file if it hasn't been created
	if not os.path.exists("galfitm.log"):
	  shutil.move("fit.log","galfitm.log")
	else:
	  with open("fit.log") as f:
	    with open("galfitm.log", "a") as f1:
		for line in f:
			f1.write(line) 

	#3. Create composed model
	make_model_ima_galfitm.main('model.galfit.01.band', subtract_sky=True, galfitPath=deca_setup.galfitmPath)
        #print 'HEREE'






  # Check the model:
  edge_disk = False
  hdulist = pyfits.open('composed_model.fits')
  for k in range(len(hdulist)):
    header = hdulist[k].header
    if header['NAME_OF_LAYER'] == 'edgedisk':
      edge_disk==True
  

  if deca_setup.code=='IMFIT':
    res_outp_file = 'imfit.01'
    inp_file = 'imfit.inp'
  elif deca_setup.code=='GALFIT':
    res_outp_file = 'galfit.01'  
    inp_file = 'galfit.inp'
  elif deca_setup.code=='GALFITM':
    res_outp_file = 'model.galfit.01.band'  
    inp_file = 'galfitm.inp'
  
  
  
  
  if deca_setup.code!='GALFITM':
        hdulist = pyfits.open('composed_model.fits')
        data = hdulist[1].data
        min_level = m0-2.5*log10(np.max(data)) + 5.
        
        #4. Create 1D and 2D profiles
        if deca_setup.plot_2D==True:
            #image_2d = plot_2d_decom_res_IMFIT.main('composed_model.fits',25.5,pix2sec,m0,borders='0,0,0,0',format_out='eps',view='column',color='color',scale_bar='sub')
            image_2d = plot_2d_profile.main('composed_model.fits',min_level,pix2sec,m0,mask_file=None,
                borders='0,0,0,0',output_file='plot_2d.png',view=deca_setup.plot_2d_allignment,color='gray_r',scale_bar=None,grid=None)
        else:
            image_2d = 'none'
        
        if deca_setup.plot_hor==True:
            zz = plot_profile.main('composed_model.fits',m0,pix2sec,profile = 'cut',xc=0.,Rmin=0.,Rmax=0.)
            image_hor = 'composed_model_prof_cut.eps'
        else:
            image_hor = 'none'
            
        if edge_disk==True:
            zz = plot_profile.main('composed_model.fits',m0,pix2sec,profile = 'vert',xc=0.,yc=0.,PA=0.,Rmin=0.,Rmax=0.)
            image_vert = 'composed_model_prof_ver.eps'
        else:
            image_vert = 'none'
        
        if deca_setup.plot_azim_aver == True:
            zz = plot_profile.main('composed_model.fits',m0,pix2sec,profile = 'azim',xc=0.,yc=0.,PA=0.,Rmin=0.,Rmax=0.)
            image_azim = 'composed_model_prof_azim.eps'
        else:
            image_azim = 'none'

        if deca_setup.plot_pix == True:
            zz = plot_profile.main('composed_model.fits',m0,pix2sec,profile = 'rad',xc=0.,yc=0.,PA=0.,Rmin=0.,Rmax=0.)
            image_pix = 'composed_model_prof_rad.eps'
        else:
            image_pix = 'none'
  else:
        if deca_setup.plot_2D==True:
            image_2d = plot_2d_profile_gm.main('composed_model.fits',None,1.,20.,mask_file=None,
                borders='0,0,0,0',output_file='plot_2d.png',view='line',color='gray_r',scale_bar=None,grid=None)      
        image_pix = 'none'
        image_azim = 'none'
        image_vert = 'none'
        image_hor = 'none'
        

  FILES = ['composed_model.fits',image_2d ,image_hor,image_vert,image_azim,image_pix,res_outp_file,inp_file,'constraint.txt']
  
  #5. Place files in the folder
  copy_images(FILES, new_dir + '/')
  
  for file in [image_2d ,image_hor,image_vert,image_azim,image_pix]:
    if os.path.exists(file): 
	os.remove(file)
  
  return new_dir,0
 except:
  return new_dir,9
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO REMOVE FILES
def remove_files_func(images,del_all=True):
  #input_image,mask_image,sigma_image,psf_image = new_images
  [input_image,sigma_image,psf_image,mask_image] = images
  

  
  if deca_setup.code=='IMFIT':
    res_outp_file = 'imfit.01'
  elif deca_setup.code=='GALFIT':
    res_outp_file = 'galfit.01'
  elif deca_setup.code=='GALFIT':
    res_outp_file = 'model.galfit.01'
    
  #7. Final remove
  if del_all==True:
    for file in ['composed_model.fits','2d_decom_res.eps',res_outp_file,'fit.log',\
      'constraint.txt',input_image,mask_image,sigma_image,psf_image,'model.galfit.01.band']:
      if isinstance(file, np.ndarray)==False:
        if os.path.exists(file): 
            os.remove(file)
      else:
          for fil in file:
            if os.path.exists(fil): 
                os.remove(fil)              
  else:
    for file in ['composed_model.fits','2d_decom_res.eps',res_outp_file,'fit.log','model.galfit.01.band']:
      if os.path.exists(file): 
	os.remove(file)      
# -----------------------------------------------------------------
