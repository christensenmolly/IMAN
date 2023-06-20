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

from astropy import units as u
from astropy.coordinates import UnitSphericalRepresentation
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord

tmp_out = sys.stdout

import deca_setup
import models
import galfit_input
import imfit_input
import initial_guess
import save_results
import threading
FNULL = open(os.devnull, 'w')
import signal
from astropy import wcs


galfitPath = deca_setup.galfitPath
imfitPath = deca_setup.imfitPath
show_code_run = deca_setup.show_code_run
timeout = deca_setup.timeout

class RunCmd(threading.Thread):
    #http://stackoverflow.com/questions/4789837/how-to-terminate-a-python-subprocess-launched-with-shell-true
    def __init__(self, cmd, timeout):
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.timeout = timeout

    def run(self):
	if show_code_run==True:
	  self.p = subprocess.Popen(self.cmd, shell=True,preexec_fn=os.setsid)#, stdout=FNULL, stderr=subprocess.STDOUT)
	else:
	  self.p = subprocess.Popen(self.cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT,preexec_fn=os.setsid)
        self.p.wait()

    def Run(self):
        self.start()
        self.join(self.timeout)

        if self.is_alive():
            #self.p.kill()#p.terminate()      #use self.p.kill() if process needs a kill -9
            os.killpg(os.getpgid(self.p.pid), signal.SIGTERM)
            self.join()
'''
def RunCmd(cmd):
  from __future__ import print_function
  from threading import Timer
  #from subprocess import Popen, PIPE, STDOUT
  #from subprocess import call

  def terminate(process):
      if process.poll() is None:
	  subprocess.call('taskkill /F /T /PID ' + str(process.pid))
  # start process
  process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
		  bufsize=1, universal_newlines=True)

  # terminate process in 15 seconds
  timer = Timer(timeout, terminate, args=[process])
  timer.start()

  # print output
  for line in iter(process.stdout.readline, ''):
      print(line, end='')
  process.stdout.close()
  process.wait() # wait for the child process to finish
  timer.cancel()
'''

def ellipse_borders(xy,width,height,angle):
  width = width 
  height = height 
  angle = np.radians(angle)
  
  X = math.sqrt( (width*math.cos(angle))**2 + (height*math.sin(angle))**2 )
  Y = math.sqrt( (width*math.sin(angle))**2 + (height*math.cos(angle))**2 )
  x_max = int(math.ceil(xy[0] + X))
  x_min = int(math.ceil(xy[0] - X))
  y_max = int(math.ceil(xy[1] + Y))
  y_min = int(math.ceil(xy[1] - Y))
  return [[x_min,y_min],[x_max,y_max]]

def check_belong(xl,yl,xr,yr,nx,ny):
  if xr<0. or yr<0. or xl>nx-1 or yl>ny-1:
    return False

  XA1 = xl; YA1 = yl; XA2 = xr; YA2 = yr
  XB1 = 0.; YB1 = 0.; XB2 = nx-1.; YB2 = ny-1.
  if XA1<0.:
    XA1 = XA1 - XA1
    XA2 = XA2 - XA1
    XB1 = XB1 - XA1
    XB2 = XB2 - XA1    

  if YA1<0.:
    YA1 = YA1 - YA1
    YA2 = YA2 - YA1
    YB1 = YB1 - YA1
    YB2 = YB2 - YA1 
  #print XA1,YA1,XA2,YA2
  #print XB1,YB1,XB2,YB2
  #SI= max([0., min([XA2, XB2]) - max([XA1, XB1]) * max([0., min([YA2, YB2]) - max([YA1, YB1])])])	# Area of intersection
  x_intersection = min([XA2, XB2]) - max([XA1, XB1]) + 1.
  y_intersection = min([YA2, YB2]) - max([YA1, YB1]) + 1.
 
  if x_intersection <= 0. or y_intersection <= 0.:
    return False
  else:
    SI = x_intersection * y_intersection
  #print SI  
  if floor(SI)>100.:	# If the intersection area is larger than 100 pixels
    return True
  else:
    return False
  
def run_galfit(file):
  if show_code_run==True:
    #subprocess.call(galfitPath+'galfit ' + file, shell=True)
    RunCmd([galfitPath+'galfit '+file], timeout).Run()
  else:
    #subprocess.call(galfitPath+'galfit ' + file, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    RunCmd([galfitPath+'galfit '+file], timeout).Run()

def run_imfit(string):
  if show_code_run==True:
    #subprocess.call(string, shell=True)
    RunCmd([string], timeout).Run()
  else:
    #subprocess.call(string, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    RunCmd([string], timeout).Run()

# -----------------------------------------------------------------
# FUNCTION TO FIND THE NEAREST ELEMENT IN THE ARRAY
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
# -----------------------------------------------------------------



# -----------------------------------------------------------------
# Function to find WCS rotation angle at the position of a SkyCoord coordinate.
def angle_at_skycoord(skycoord, wcs, offset=1. * u.arcsec):
    """
    Calculate the pixel scale and WCS rotation angle at the position of
    a SkyCoord coordinate.
    Parameters
    ----------
    skycoord : `~astropy.coordinates.SkyCoord`
        The SkyCoord coordinate.
    wcs : `~astropy.wcs.WCS`
        The world coordinate system (WCS) transformation to use.
    offset : `~astropy.units.Quantity`
        A small angular offset to use to compute the pixel scale and
        position angle.
    Returns
    -------
    scale : `~astropy.units.Quantity`
        The pixel scale in arcsec/pixel.
    angle : `~astropy.units.Quantity`
        The angle (in degrees) measured counterclockwise from the
        positive x axis to the "North" axis of the celestial coordinate
        system.
    Notes
    -----
    If distortions are present in the image, the x and y pixel scales
    likely differ.  This function computes a single pixel scale along
    the North/South axis.
    """

    # We take a point directly "above" (in latitude) the input position
    # and convert it to pixel coordinates, then we use the pixel deltas
    # between the input and offset point to calculate the pixel scale and
    # angle.

    # Find the coordinates as a representation object
    coord = skycoord.represent_as('unitspherical')

    # Add a a small perturbation in the latitude direction (since longitude
    # is more difficult because it is not directly an angle)
    coord_new = UnitSphericalRepresentation(coord.lon, coord.lat + offset)
    coord_offset = skycoord.realize_frame(coord_new)

    # Find pixel coordinates of offset coordinates and pixel deltas
    x_offset, y_offset = skycoord_to_pixel(coord_offset, wcs, mode='all')
    x, y = skycoord_to_pixel(skycoord, wcs, mode='all')
    dx = x_offset - x
    dy = y_offset - y

    scale = offset.to(u.arcsec) / (np.hypot(dx, dy) * u.pixel)
    angle = (np.arctan2(dy, dx) * u.radian).to(u.deg)

    return angle

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
  x = np.linspace(max([sma/30.,2.]), sma, num=1000)
  re = x[find_nearest(f2(x),sum5/2.)]
  return re,sum5,data[int(yc),int(xc)]
# -----------------------------------------------------------------  
  
  
  
# -----------------------------------------------------------------  
# Function to retrieve results from the imfit output file
def read_galfit(galf_output_file):
	f = open(galf_output_file,'r')
	lines = f.readlines()
	
	XC = []; YC = []; Mag = []; Re = []; N = []; Q = []; PAA = []
	for kk in range(len(lines)):
	  line = lines[kk]
	  if '#  Position x, y' in line:
	   if ' 0) sersic                 #  Component type' in lines[kk-1]:
	    xcc = float(line.split()[1])
	    ycc = float(line.split()[2])
	    XC.append(xcc)
	    YC.append(ycc)
	    Mag.append(float(lines[kk+1].split()[1]))
	    Re.append(float(lines[kk+2].split()[1]))
	    N.append(float(lines[kk+3].split()[1]))
	    Q.append(float(lines[kk+7].split()[1]))
	    PAA.append(float(lines[kk+8].split()[1]))
	f.close()
	return XC,YC,Mag,Re,N,Q,PAA
# -----------------------------------------------------------------  



# -----------------------------------------------------------------
# MAIN FUNCTION
def main(file_with_galaxies,new_images,observation_info,object_info,add_info,log_text,keys,del_files=True): 
  #*******************************************************************          
  [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,Sky_level,SkySubtr,sampling,convolution_box,FILTER] = observation_info
  [xc,yc,name,distance,Aext,Kext] = object_info
  [number,fitting_proc,mode,model,region_to_fit] = add_info
  [input_image,sigma_image,psf_image,mask_image] = new_images
  [control,delay] = keys  
  #*******************************************************************
  
  # Read in region file with selected objects
  [xc,yc],[xc_star,yc_star] = misc_functions.READ_REGION_OBJECTS(file_with_galaxies, input_image)
  

  ## If the exposure time is not specified, it will be set to 1.
  if np.isnan(exptime)==True:
    exptime = 1.

  ## Read the input galaxy image
  hdulist = pyfits.open(input_image)
  data = hdulist[0].data

  if mask_image!='none' and os.path.exists(mask_image):
      ## Read the input mask and convert it to astropy mask
      hdulist_mask = pyfits.open(mask_image)
      data_mask = hdulist_mask[0].data

      mask_astropy = np.zeros_like(np.array(data,dtype=float),dtype=bool)

      for k in range(ny):
	for i in range(nx):
	  if data_mask[k,i]!=0.:
	    mask_astropy[k,i] = True
  else:
      mask_astropy = np.zeros_like(np.array(data,dtype=float),dtype=bool) 


  try:      
        # -----------------------------------------------------------------
        #  ***************FIT WITH A SINGLE SERSIC FUNCTION****************
	if deca_setup.code=='GALFIT':
	    ## Create input galfit file
	    f = open('galfit.inp', "w") 
	    sys.stdout = f  
	    galfit_input.header(input_image,sigma_image,psf_image,mask_image,'none','model.fits',1,nx,1,ny,m0-2.5*log10(exptime),scale,int(sampling))
	    
	    for k in range(len(xc)):
		## Coarsely estimate the single sersic parameters
		aperture = EllipticalAperture([(xc[k],yc[k])], ellA[k],ellB[k],radians(ellPA[k]))
		phot_table = aperture_photometry(data, aperture, mask=mask_astropy)
		if SkySubtr==1:
		  lum = phot_table['aperture_sum'][0]
		else:
		  bkg_sum = Sky_level * aperture.area()
		  lum = phot_table['aperture_sum'][0] - bkg_sum

		## Define the single Sersic model:
		galfit_input.Sersic(component=k+1,xc=xc[k],yc=yc[k],meb=float('nan'),
		      Mb = m0 - 2.5*log10(lum),reb=ellA[k]/5.,
		      n=2.,q=ellB[k]/ellA[k],PA=ellPA[k])

	    for k in range(len(xc_stars)):    
		## Coarsely estimate the single sersic parameters
		aperture = EllipticalAperture([(xc_stars[k],xc_stars[k])], 2.*fwhm,2.*fwhm,0.)
		phot_table = aperture_photometry(data, aperture, mask=mask_astropy)
		if SkySubtr==1:
		  lum = phot_table['aperture_sum'][0]
		else:
		  bkg_sum = Sky_level * aperture.area()
		  lum = phot_table['aperture_sum'][0] - bkg_sum

		## Define the single Sersic model:
		galfit_input.Agn(component=len(xc)+k+1,xc=xc_stars[k],yc=yc_stars[k],mag=m0 - 2.5*log10(lum))
		
	    galfit_input.Sky(component=len(xc)+len(xc_stars)+1,sky_level=Sky_level)
	    sys.stdout = tmp_out
	    f.close()

	    os.chmod(r"galfit.inp",0777)
	    
	    ## Fit data with a single Sersic model
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
		  return 1,log_text
	    
	    #  ***************FIT A MORE COMPLEX MODEL TO THE DATA****************
	    # TODO: Write the same for IMFIT!
	    if os.path.exists('galfit.01') and model!="sersic":
                XC,YC,Mag,Re,N,Q,PAA = read_galfit('galfit.01')

                # Create new input galfit file with a more complex model
                f = open('galfit.inp', "w") 
                sys.stdout = f  
                galfit_input.header(input_image,sigma_image,psf_image,mask_image,'none','model.fits',1,nx,1,ny,m0-2.5*log10(exptime),scale,int(sampling))

                
                for kk in range(len(XC)):
                    #Galaxies:                       
                    if model=='exp_disc':  
                        galfit_input.ExpDisc(component=kk+1,xc=XC[kk],yc=YC[kk],m0d=float('nan'),Md=Mag[kk],h=Re[kk] / 1.68,q=Q[kk],PA=PAA[kk])

                    if model=='eon_disc':
                        galfit_input.EdgeDisc(component=kk+1,xc=XC[kk],yc=YC[kk],m0d=initial_guess.m0_disc_edge_on_f(Mag[kk],Re[kk] / 1.68,Re[kk] / 1.68 / 4.),z0=Re[kk] / 1.68 / 4.,h=Re[kk] / 1.68,PA=PAA[kk])
                        
                    if model=='sersic+exp_disc':
                        galfit_input.Sersic(component=kk+1,xc=XC[kk],yc=YC[kk],meb=float('nan'),Mb = Mag[kk]+0.75,reb=0.1*Re[kk],n=2.,q=0.8,PA=PAA[kk])
                        galfit_input.ExpDisc(component=kk+2,xc=XC[kk],yc=YC[kk],m0d=float('nan'),Md=Mag[kk]+0.65,h=Re[kk] / 1.68,q=Q[kk],PA=PAA[kk])
                        
                    if model=='sersic+eon_disc':
                        galfit_input.Sersic(component=kk+1,xc=XC[kk],yc=YC[kk],meb=float('nan'),Mb = Mag[kk]+0.75,reb=0.1*Re[kk],n=2.,q=0.8,PA=PAA[kk])
                        galfit_input.EdgeDisc(component=kk+2,xc=XC[kk],yc=YC[kk],m0d=Mag[kk] + 0.65 + 2.5*log10(2.*math.pi*(Re[kk] / 1.68)*scale*(Re[kk] / 1.68 / 4.)*scale) + 5.*log10(scale),z0=Re[kk] / 1.68 / 4.,h=Re[kk] / 1.68,PA=PAA[kk]) 
                        
                    if model=='agn+exp_disc':
                        galfit_input.Agn(component=kk+1,xc=XC[kk],yc=YC[kk],mag=100. - 2.5*log10( (10**(0.4*(100.-Mag[kk]))/100.) ))
                        galfit_input.ExpDisc(component=kk+2,xc=XC[kk],yc=YC[kk],m0d=float('nan'),Md=Mag[kk],h=Re[kk] / 1.68,q=Q[kk],PA=PAA[kk])   
                        
                    if model=='agn+eon_disc':
                        galfit_input.Agn(component=kk+1,xc=XC[kk],yc=YC[kk],mag=100. - 2.5*log10( (10**(0.4*(100.-Mag[kk]))/100.) ))
                        galfit_input.EdgeDisc(component=kk+2,xc=XC[kk],yc=YC[kk],m0d=initial_guess.m0_disc_edge_on_f(Mag[kk],Re[kk] / 1.68,Re[kk] / 1.68 / 4.),z0=Re[kk] / 1.68 / 4.,h=Re[kk] / 1.68,PA=PAA[kk])         
                        
                    if model=='sersic+sersic':
                        galfit_input.Sersic(component=kk+1,xc=XC[kk],yc=YC[kk],meb=float('nan'),Mb = Mag[kk]+0.75,reb=0.1*Re[kk],n=2.,q=0.8,PA=PAA[kk])
                        galfit_input.Sersic(component=kk+2,xc=XC[kk],yc=YC[kk],meb=float('nan'),Mb = Mag[kk]+0.65,reb=Re[kk],n=1.,q=Q[kk],PA=PAA[kk])                                               

                for kk in range(len(XC_stars)):    
                    #Stars:
                    galfit_input.Agn(component=len(XC)+kk+1,xc=XC_stars[kk],yc=YC_stars[kk],mag=Mag_stars[kk])
                    
                galfit_input.Sky(component=len(XC)+len(XC_stars)+1,sky_level=Sky_level)
                sys.stdout = tmp_out
                f.close()

                os.chmod(r"galfit.inp",0777)
                
                ## Fit data with the single Sersic model
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
                        return 1,log_text
                    
	    
	    # Read results:
	    #hdulist_model = pyfits.open('model.fits')
	    #model_header = hdulist_model[2].header
	    
	    #PARS = models.read_results_First_Ser(1,PARS,model_header)
	
	elif deca_setup.code=='IMFIT':
	    ## Create input imfit file
	    run_line = imfit_input.imfit_exec_line(new_images,observation_info,object_info,'imfit.inp')

	    f = open('imfit.inp', "w") 
	    sys.stdout = f  
	    imfit_input.header(run_line,0.,0.)

	    for k in range(len(xc)):
		## Coarsely estimate the single sersic parameters
		aperture = EllipticalAperture([(xc[k],yc[k])], ellA[k],ellB[k],radians(ellPA[k]))
		phot_table = aperture_photometry(data, aperture, mask=mask_astropy)
		if SkySubtr==1:
		  lum = phot_table['aperture_sum'][0]
		else:
		  bkg_sum = Sky_level * aperture.area()
		  lum = phot_table['aperture_sum'][0] - bkg_sum

		## Define the single Sersic model:
		imfit_input.Sersic(component=k+1,xc=xc[k],yc=yc[k],meb=float('nan'),
		      Mb = m0 - 2.5*log10(lum),reb=ellA[k]/5.,
		      n=2.,q=ellB[k]/ellA[k],PA=ellPA[k],m0=m0,pix2sec=scale,write_coords=False)

	    
	    imfit_input.Sky(component=len(xc)+1,sky_level=Sky_level)
	    sys.stdout = tmp_out
	    f.close()

	    os.chmod(r"imfit.inp",0777)
	    
	    ## Fit data with the single Sersic model
	    if os.path.exists('model.fits'):
	      os.remove('model.fits')
	    run_imfit(run_line)

	    while not os.path.exists('model.fits'):
	      if mode=='control':
		    print 'IMFIT crashed. Please change the input file!'
		    subprocess.call(deca_setup.text_editor + ' imfit.inp' , shell=True)
		    raw_input("Press Enter to continue...")
		    run_imfit(run_line)
	      else:
		  return 1,log_text

	    # Read results:    
	    #PARS = models.read_results_First_Ser_IMFIT(1,PARS,'imfit.01',m0,scale)
	
	save_results.main(new_images,scale,m0)
	if del_files==True:
	  save_results.remove_files_func(new_images, del_all=False)
	for file in ['fit.log','galfit.01','imfit.01']:
	  if os.path.exists(file): 
	    os.remove(file)
	if model=='sersic':
	  return 0,log_text
        # -----------------------------------------------------------------


  except:
    return 1,log_text


