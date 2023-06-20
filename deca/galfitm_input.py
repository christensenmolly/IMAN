#!/usr/bin/python
# Module to work with GalfitM
# *****************************************************************
# **       DECA -- DEComposition Analysis of galaxy images       **
# **          Astronomical Observatory, Ghent University         **
# *****************************************************************

# Import standard packages
import random as random_number
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
from numpy import *
from pylab import *
import os
import shutil
import signal

# Import astro packages
import pyfits

# Import DECA modules
import deca_setup

# -----------------------------------------------------------------
# FUNCTION TO GET WAVELENGTH BY A BAND NAME
def band_to_wavelength(band):
    if band=='J':
      wavelength = 1.24
    if band=='H':
      wavelength = 1.66
    if band=='Ks':
      wavelength = 2.16   
    if band=='FUV':
      wavelength = 0.1528
    if band=='NUV':
      wavelength = 0.2271
    if band=='w1':
      wavelength = 3.35
    if band=='w2':
      wavelength = 4.60
    if band=='w3':
      wavelength = 11.56
    if band=='w4':
      wavelength = 22.09
    if band=='Z':
      wavelength = 0.882  
    if band=='Y':
      wavelength = 1.031  
      
    if band=='K':
      wavelength = 2.201      
    if band=='u':
	wavelength = 0.354
    if band=='g':
	wavelength = 0.475
    if band=='r':
	wavelength = 0.622
    if band=='i':
	wavelength = 0.763
    if band=='z':
	wavelength = 0.905
	
    if band=='B':
	wavelength = 0.445
    if band=='R':
	wavelength = 0.658
    if band=='V':
	wavelength = 0.551

    if band=='I1':
       wavelength = 3.550
    if band=='I2':
       wavelength = 4.493       
    if band=='I3':
       wavelength = 5.731  
    if band=='I4':
       wavelength = 7.872
     
    if band=='M1':
       wavelength = 23.68
    if band=='M2':
       wavelength = 71.42
    if band=='M3':
       wavelength = 155.9

    if band == 'P100':
      wavelength = 100.
    if band == 'P160':
      wavelength = 160.      
    if band == 'S250':
      wavelength = 250.
    if band == 'S350':
      wavelength = 350.
    if band == 'S500':
      wavelength = 500.
    return str(int(wavelength*10000.))	# in A
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION GET CHEBYSHEV'S ORDER
def get_cheb_order(comp_name, par_numb, number_of_bands):
    cheb_order = 1
    
    if comp_name=='sersic':
        if par_numb=='3)':
            cheb_order = number_of_bands
        if par_numb=='4)':
            cheb_order = 3#max([2,number_of_bands - 1])
        if par_numb=='5)':
            cheb_order = 3#max([2,number_of_bands - 1])           
        if par_numb=='9)':
            cheb_order = deca_setup.chb_bulge_q
        if par_numb=='10)':
            cheb_order = deca_setup.chb_bulge_PA      
            
    if comp_name=='expdisk':
        if par_numb=='3)':
            cheb_order = number_of_bands
        if par_numb=='4)':
            cheb_order = deca_setup.chb_disc_rs         
        if par_numb=='9)':
            cheb_order = deca_setup.chb_disc_q
        if par_numb=='10)':
            cheb_order = deca_setup.chb_disc_PA     
            
    if comp_name=='edgedisk':
        if par_numb=='3)':
            cheb_order = number_of_bands
        if par_numb=='4)':
            cheb_order = deca_setup.chb_disc_z0
        if par_numb=='5)':
            cheb_order = deca_setup.chb_disc_rs              
        if par_numb=='10)':
            cheb_order = deca_setup.chb_disc_PA   
    
    if comp_name=='sky':
	 cheb_order = 0
    return cheb_order
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO CONVERT BAND TO WAVELENGTH
def bands_to_wavelengths(bands):
    wavelengths = []
    for band in bands:
      wavelengths.append(band_to_wavelength(band))
    return wavelengths  
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO CREATE HEADER FOR GALFITM INPUT FILE
def header(file_Images,bands,file_constraint,
	   file_galfit_outimage,xmin,xmax,ymin,ymax,magZPs,generalScaleValue,sampling,convolution_box):
	
	GALAXY_IMAGES = []; SIGMA_IMAGES = []; PSF_IMAGES = []; MASK_IMAGES = []
	for k in range(len(file_Images)):
	  GALAXY_IMAGE,SIGMA_IMAGE,PSF_IMAGE,MASK_IMAGE = file_Images[k]
	  GALAXY_IMAGES.append(GALAXY_IMAGE)
	  SIGMA_IMAGES.append(SIGMA_IMAGE)
	  PSF_IMAGES.append(PSF_IMAGE)
	  MASK_IMAGES.append(MASK_IMAGE)
	  
	  
	file_Images =  ','.join(GALAXY_IMAGES)
	wavelengths = bands_to_wavelengths(bands)
	wavelengths =  ','.join(wavelengths)
	bands =  ','.join(bands)
	file_Sigmas =  ','.join(SIGMA_IMAGES)
	file_Masks =  ','.join(MASK_IMAGES)
	magZPs =  ','.join(magZPs)

	#How to choose convolution box:
	psf_box_x = []
	psf_box_y = []
	for k in range(len(PSF_IMAGES)):
	  hdulist = pyfits.open(PSF_IMAGES[k])
	  image = hdulist[0].data    
	  ySize_psf, xSize_psf = image.shape
	  psf_box_x.append(xSize_psf)
	  psf_box_y.append(ySize_psf)
	psf_box = min([min(psf_box_x),min(psf_box_y)])
	file_Psfs =  ','.join(PSF_IMAGES)
	
	if file_constraint==None:
	  file_constraint = 'none'
	
	#if convolution_box!=None:
        #    psf_box = convolution_box
	#psf_box = min([min([ySize_psf, xSize_psf]),ySize])
	#psf_box = min([ySize_psf, xSize_psf,xmax-xmin+5])
	print "\n==============================================================================="
	print "# IMAGE and GALFIT CONTROL PARAMETERS"
	print "A) %s                  # Input data images (FITS file)" % (file_Images)
	print "A1) %s                  # Band labels " % (bands)
	print "A2) %s                  # Band wavelengths " % (wavelengths)
	print "B) %s                  # Output data image block" % (file_galfit_outimage)
	print "C) %s                # Sigma image names (made from data if blank or none)" % (file_Sigmas)
	print "D) %s                # Input PSF images and (optional) diffusion kernel" % (file_Psfs)
	print "E) %s                   # PSF fine sampling factor relative to data" % (sampling[0])
	print "F) %s                # Bad pixel masks (FITS images)" % (file_Masks)
	print "G) %s                # File with parameter constraints (ASCII file)" % (file_constraint)
	print "H) %i    %i   %i    %i   # Image region to fit (xmin xmax ymin ymax)" % (xmin, xmax, ymin, ymax)
	print "I) %i    %i              # Size of the convolution box (x y)" % (psf_box, psf_box)
	print "J) %s              # Magnitude photometric zeropoints" % (magZPs)
	print "K) %.3f  %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (generalScaleValue[0],generalScaleValue[0])
	print "O) regular             # Display type (regular, curses, both)"
	print "P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps"
	print "U) %s                   # Choose: 0=standard parametric" % (deca_setup.param_method)
	print "V) %s                   # Use standard Levenburg-Marquardt algorithm" % (deca_setup.fit_method)
	print "W) default             # == blank,input,model,residual and assumed if omitted\n"
	
	print "# INITIAL FITTING PARAMETERS"
	print "#"
	print "#   For object type, the allowed functions are:" 
	print "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat," 
	print "#       ferrer, powsersic, sky, and isophote." 
	print "#"  
	print "#   Hidden parameters will only appear when they're specified:"
	print "#       C0 (diskyness/boxyness)," 
	print "#       Fn (n=integer, Azimuthal Fourier Modes),"
	print "#       R0-R10 (PA rotation, for creating spiral structures)."
	print "#" 
	print "# -----------------------------------------------------------------------------"
	print "#   par)    par value(s)    fit toggle(s)    # parameter description" 
	print "# -----------------------------------------------------------------------------\n"
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO MERGE COMPONENTS IN DIFFERENT BANDS
def merge_components_data(galfit_input_files, number_of_bands):
  # Check if all files have the same set of components
  Components = []
  for k in range(len(galfit_input_files)):
    galfit_file = galfit_input_files[k]
    f = open(galfit_file,'r')
    lines = f.readlines()
    components = []
    first_lines = []
    last_lines = []
    for i in range(len(lines)):
      if "0)" in lines[i] and "1)" in lines[i+1]:
	first_lines.append(i)
      if "Z)" in lines[i]:
	last_lines.append(i)
    for i in range(len(first_lines)):
      components.append(lines[first_lines[i]:last_lines[i]])
    Components.append(components)

  # Names of components:
  names_of_comps = []
  for i in range(len(Components[0])):
      line = Components[0][i]
      for ii in range(len(line)-1):
	if '0)' in line[ii] and '1)' in line[ii+1]:
	  names_of_comps.append(line[ii].split()[1])
  
  #print names_of_comps,len(Components[0])
  #exit()
  # Coordinates:
  x_coord = []
  y_coord = []
  for k in range(len(Components)):
    for i in range(len(Components[k])):
      line = Components[k][i]
      for ii in range(len(line)-1):
	if "1)" in line[ii]:
	 if "sky" not in line[ii-1]:
	  x_coord.append(float(line[ii].split()[1]))
	  if line[ii].split()[3] == '0' or line[ii].split()[3] == '1':
	    y_coord.append(float(line[ii].split()[2]))
	  if "2)" in line[ii+1]:
	    y_coord.append(float(line[ii+1].split()[1]))  
  #print y_coord
  #exit()
  if len(x_coord)>=3:
    x_med = np.median(x_coord)
    y_med = np.median(y_coord)
  elif len(x_coord)==2:
    x_med = np.mean(x_coord)
    y_med = np.mean(y_coord)
  else:
    x_med = x_coord[0]
    y_med = y_coord[0]
  
  #print x_coord
  #print y_coord
  #exit()

  # Merge content of each component's description
  FINAL_COMPONENTS = []
  FINAL_COMMENTS = []
  FINAL_NAMES = []
  for k in range(len(Components)):
    COMPONENTS = []
    NAMES = []
    COMMENTS = []
    for i in range(len(Components[k])):
      line = Components[k][i]
      for ii in range(len(line)):
	if '3)' in line[ii]:
	  if 'sky' not in line[ii-3]:
	    first_line = ii
	  else:
	    first_line = ii - 2
      COmponents = []
      NAmes = []
      COmments = []
      for ii in range(first_line,len(line),1):
	COmponents.append(float(line[ii].split()[1]))
	NAmes.append(line[ii].split()[0])
	COmments.append(line[ii].split('#')[-1].split('\n')[0])
      COMPONENTS.append(COmponents)
      NAMES.append(NAmes)
      COMMENTS.append(COmments)
    FINAL_COMPONENTS.append(COMPONENTS)
    FINAL_COMMENTS.append(COMMENTS)
    FINAL_NAMES.append(NAMES)
  

  # Final preparation of the component's part
  Number_of_files = len(FINAL_COMPONENTS)
  Number_of_components = len(FINAL_COMPONENTS[0])
  Number_of_pars_in_each_comp = []
  for k in range(len(FINAL_COMPONENTS[0])):
      Number_of_pars_in_each_comp.append(len(FINAL_COMPONENTS[0][k])) 
  #print Number_of_pars_in_each_comp

  NEW_COMPONENTS = []
  NEW_COMMENTS = []
  NEW_NAMES = []
  for k in range(Number_of_files):
    NEW_COMPONENTS.append(list(flatten(FINAL_COMPONENTS[k])))  # -> convert matrix to one line
    NEW_COMMENTS.append(list(flatten(FINAL_COMMENTS[k])))
    NEW_NAMES.append(list(flatten(FINAL_NAMES[k])))
    
  
  F_COMPONENTS = []


  for i in range(len(NEW_COMPONENTS[0])):
    F = ''
    for k in range(len(NEW_COMPONENTS)):
      if k==len(NEW_COMPONENTS)-1:
	F = F + str(NEW_COMPONENTS[k][i])
      else:
	F = F + str(NEW_COMPONENTS[k][i])+','
    F_COMPONENTS.append(F)
  #print F_COMPONENTS
  
  TWO_ARRAYS_COMP = []
  TWO_ARRAYS_COMMENTS = []
  TWO_ARRAYS_NAMES = []
  start = 0
  for k in range(len(Number_of_pars_in_each_comp)):
    end = start + Number_of_pars_in_each_comp[k]
    TWO_ARRAYS_COMP.append(F_COMPONENTS[start:end])
    TWO_ARRAYS_COMMENTS.append(NEW_COMMENTS[0][start:end])
    TWO_ARRAYS_NAMES.append(NEW_NAMES[0][start:end])
    start = start + Number_of_pars_in_each_comp[k]
  
  
  # Printing
  for k in range(len(TWO_ARRAYS_COMP)):
        #if names_of_comps[k]!='sky':
	print "# Component number: %i" % (k+1)
        print "0) %s               #  Component type" % (names_of_comps[k])
        if names_of_comps[k]!='sky':
	  print "1) %.3f  %i    # position x [pixel]" % (x_med,deca_setup.chb_coord_x)
	  print "2) %.3f  %i    # position y [pixel]" % (y_med,deca_setup.chb_coord_y)
        for i in range(len(TWO_ARRAYS_COMP[k])):
            cheb_order = get_cheb_order(names_of_comps[k], TWO_ARRAYS_NAMES[k][i], number_of_bands)             
            if TWO_ARRAYS_NAMES[k][i]=='9)' or TWO_ARRAYS_NAMES[k][i]=='10)':
		TWO_ARRAYS_COMP[k][i] = str(np.mean(np.array(TWO_ARRAYS_COMP[k][i].split(','),dtype=float)))
            print "%s %s %i # %s" % (TWO_ARRAYS_NAMES[k][i],TWO_ARRAYS_COMP[k][i],cheb_order,TWO_ARRAYS_COMMENTS[k][i])
        print "Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"
        print "\n\n"
	#else:
	#  print TWO_ARRAYS_COMP[k]
	#  print TWO_ARRAYS_NAMES
	
    
  
  return Components
# -----------------------------------------------------------------
