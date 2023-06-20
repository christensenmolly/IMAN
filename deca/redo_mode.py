#!/usr/bin/env python
# -*- coding: utf8 -*-
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

import pyfits
from photutils import datasets
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import EllipticalAperture
from astropy.table import hstack
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

tmp_out = sys.stdout

import deca_setup
import models
import galfit_input
import initial_guess
import save_results
FNULL = open(os.devnull, 'w')

galfitPath = deca_setup.galfitPath

def run_galfit(file):
  if show_galf_run==True:
    subprocess.call(galfitPath+'galfit ' + file, shell=True)
  else:
    subprocess.call(galfitPath+'galfit ' + file, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

# -----------------------------------------------------------------
# FUNCTION TO FIND FLUXES WITHIN AN ELLIPSE TO ESTIMATE TOTAL LUMINOSITY AND EFFECTIVE RADIUS
def find_eff_radius(data,mask_data,xc,yc,sma,smb,theta):
  Sma = [max([sma/30.,2.]),sma/20.,sma/10.,sma/2.,sma]
  Sma = np.array(Sma)
  Smb = Sma * smb/sma
  fluxes = []

  for k in range(len(Sma)):
    apertures = EllipticalAperture([(xc,yc)], Sma[k],Smb[k],radians(theta))
    fluxes.append(aperture_photometry(data, apertures, mask=mask_data))
  phot_table = hstack(fluxes)

  sum1 = float(phot_table['aperture_sum_1'])
  sum2 = float(phot_table['aperture_sum_2'])
  sum3 = float(phot_table['aperture_sum_3'])
  sum4 = float(phot_table['aperture_sum_4'])
  sum5 = float(phot_table['aperture_sum_5'])
  Sum = np.array([sum1,sum2,sum3,sum4,sum5])

  f2 = interp1d(Sma, Sum)

  #plt.plot(Sma,Sum, 'o')
  x = np.linspace(max([sma/30.,2.]), sma, num=1000)
  
  re = x[find_nearest(f2(x),sum5/2.)]
  
  #plt.plot(x, f2(x))
  #plt.plot([re],[f2(re)],'o')
  #plt.show()
  return re,sum5,data[int(yc),int(xc)]
  
  #return 
# -----------------------------------------------------------------  
  


# -----------------------------------------------------------------
# MAIN FUNCTION
def main(file_with_galaxies,new_images,observation_info,object_info,add_info):
  #input_image,mask_image,sigma_image,psf_image = new_images
  #[scale,m0,gain,read_out_noise,ncombine,exptime,fwhm,SkySubtr] = ADD_INFO
  #[number,FILTER,name,sampling,distance,fitting_proc] = INFO

  [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,Sky_level,SkySubtr,sampling,convolution_box,FILTER] = observation_info
  [xc,yc,name,distance,Aext,Kext] = object_info
  [number,fitting_proc,mode,model,region_to_fit] = add_info
  [input_image,sigma_image,psf_image,mask_image] = new_images
  
  if np.isnan(exptime)==True:
    exptime = 1.

  if not os.path.exists("galfit.inp"):
        print 'File galfit.inp does not exist'
        return 1
  else:
        os.chmod(r"galfit.inp",0777)

        repeat = True
        while repeat==True:
	    if os.path.exists('model.fits'):
		os.remove('model.fits')
	    run_galfit('galfit.inp')

	    while not os.path.exists('model.fits'):
		if mode=='control':
		      print 'Galfit crashed. Please change the input file!'
		      subprocess.call(deca_setup.text_editor + ' galfit.inp' , shell=True)
		      raw_input("Press Enter to continue...")
		      run_galfit('galfit.inp')
		else:
		      return 1

            # Read results:
            hdulist_model = pyfits.open('model.fits')
            model_header = hdulist_model[2].header
            
            #PARS = models.read_results_First_Ser(1,PARS,model_header)
            #PARS = models.read_results_First_Disc(2,PARS,model_header)

            new_dir = save_results.main(new_images,scale,m0)              
            
            if mode=='control':
                            if os.path.exists(new_dir + '/' + '2d_decom_res.eps'):
                                subprocess.call(deca_setup.image_viewer + ' ' + new_dir + '/' + '2d_decom_res.eps' , shell=True) 
                            if os.path.exists(new_dir + '/' + 'composed_model_prof_hor.eps'):
                                subprocess.call(deca_setup.image_viewer + ' ' + new_dir + '/' + 'composed_model_prof_hor.eps' , shell=True)                 
                            if os.path.exists(new_dir + '/' + 'composed_model_prof_vert.eps'):
                                subprocess.call(deca_setup.image_viewer + ' ' + new_dir + '/' + 'composed_model_prof_vert.eps' , shell=True)
                            Repeat = 'yes'
                            Repeat = str(raw_input('Are you happy with the results?  (%s): ' % (str(Repeat)))) or str(Repeat)
			    if Repeat == 'no':
				  if os.path.exists('galfit.inp'):
				      subprocess.call(deca_setup.text_editor + ' galfit.inp' , shell=True)
				      raw_input("Press Enter to continue...")
				  repeat = True
				  save_results.remove_files_func(new_images, del_all=False)
			    else:
				  repeat = False
				  save_results.remove_files_func(new_images, del_all=True)
            else:
                            repeat = False
                            save_results.remove_files_func(new_images, del_all=True)
  return 0
