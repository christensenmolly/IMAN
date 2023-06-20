#!/usr/bin/python
# -*- coding:  cp1251 -*-

#import sys
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
import shutil
import subprocess
import random
from numpy import fft
import pyfits
import re
import os

from photutils import datasets
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import EllipticalAperture
from photutils import CircularAperture

import imp_setup
import imp_center
import rebin_image
import crop_psf_star
import imp_psf

FNULL = open(os.devnull, 'w')
box_psf = imp_setup.box_psf
beta = imp_setup.beta
window = imp_setup.window


tmp_out = sys.stdout

def backgr_around_star(scidata,xc,yc,Rin,Rout=None):
  nx, ny =scidata.shape[1], scidata.shape[0]
  if Rout==None:
    Rin = int(ceil(1.1*Rin))
    Rout = Rin + max([imp_setup.Npix_around_star,int(ceil(1.3*Rin))])
  else:
    Rin = int(ceil(1.1*Rin))
    Rout = int(ceil(Rout))    
  backgr = []
  for y in range(int(yc)-Rout,int(yc)+Rout,1):
    for x in range(int(xc)-Rout,int(xc)+Rout,1):
      if (x-int(xc))**2+(y-int(yc))**2>Rin**2 and (x-int(xc))**2+(y-int(yc))**2<=Rout**2 and x>=0 and y>=0 and x<nx and y<ny:
        backgr.append(scidata[y,x])
  mean, median, std = sigma_clipped_stats(backgr, sigma=3.0, iters=5)
  return median,std


def bad_pix_mask(image_file,file_segm,file_badpix,xc,yc,radius,out):
            print('Masking contaminants...')
            # FUNCTION TO CREATE MASK FILE FOR GALFIT
            hdulist1 = pyfits.open(image_file)
            data = hdulist1[0].data

            shutil.copy(file_segm,file_badpix) 
            hdulist3 = pyfits.open(file_badpix, do_not_scale_image_data=True, mode='update')
            img3 = hdulist3[0].data
            for l in range(len(xc)):
                Radius = int(radius[l]*1.2)
                for k in range(-Radius,Radius):
                        for i in range(-Radius,Radius):
                            try:
                                if img3[yc[l]+k,xc[l]+i]==img3[yc[l],xc[l]]:
                                    img3[yc[l]+k,xc[l]+i] = data[yc[l]+k,xc[l]+i]
                            except:
                                zz = 1

            hdulist3.flush()

def sky_subtration(input_image, output_image, sky_level):
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    header = hdulist[0].header
    outData = data - sky_level
    outHDU = pyfits.PrimaryHDU(outData, header=header)
    outHDU.writeto(output_image,clobber=True)



# ********************************* Simulation with GALFIT *********************************
def mod_psf(window,fwhm_psf,M_psf,ell,PA,m0,beta,pix2secx,pix2secy):
	nx_psf = box_psf
	ny_psf = box_psf
	if nx_psf%2==0:
		xc_psf = int(nx_psf/2. + 1)
	else:
		xc_psf = int(nx_psf/2. + 0.5)
	if ny_psf%2==0:
		yc_psf = int(ny_psf/2. + 1)
	else:
		yc_psf = int(ny_psf/2. + 0.5)

	#========================================
	#pix2secx = pix2sec
	#pix2secy = pix2sec
	#========================================


	f = open(r"modelPSF.txt", "w") 
	sys.stdout = f
	print "\n==============================================================================="
	print "# IMAGE and GALFIT CONTROL PARAMETERS"
	print "A) none                # Input data image (FITS file)"
	print "B) psf.fits         # Output data image block"
	print "C) none                # Sigma image name (made from data if blank or none)" 
	print "D) none                # Input PSF image and (optional) diffusion kernel"
	print "E) 1                   # PSF fine sampling factor relative to data" 
	print "F) none                # Bad pixel mask (FITS image or ASCII coord list)"
	print "G) none                # File with parameter constraints (ASCII file)" 
	print "H) 1    %i   1    %i   # Image region to fit (xmin xmax ymin ymax)" % (nx_psf, ny_psf)
	print "I) %.3f    %.3f        # Size of the convolution box (x y)" % (0, 0)
	print "J) %.3f              # Magnitude photometric zeropoint" % (m0)
	print "K) %.3f    %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (pix2secx,pix2secy)
	print "O) regular             # Display type (regular, curses, both)"
	print "P) 1                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n"

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

	if window=='gauss':
		print "# Gaussian function\n"
		print "0) gaussian           # object type"
		print "1) %.3f  %.3f  0 0   # position x, y        [pixel]" % (xc_psf,yc_psf)
		print "3) %.3f       0        # total magnitude" % (M_psf)     
		print "4) %.3f       0        #   FWHM               [pixels]" % (fwhm_psf)
		print "9) %.3f        0       # axis ratio (b/a)" % (1.-ell)  
		print "10) %.3f         0       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA)
		print "Z) 0                  # leave in [1] or subtract [0] this comp from data"
		print "\n================================================================================"

	if window=='moffat':
		print "# Moffat function\n"
		print "0) moffat           # object type"
		print "1) %.3f  %.3f  0 0   # position x, y        [pixel]" % (xc_psf,yc_psf)
		print "3) %.3f       0        # total magnitude" % (M_psf)     
		print "4) %.3f       0        #   FWHM               [pixels]" % (fwhm_psf)
		print "5) %.3f        0       # powerlaw" % (beta)
		print "9) %.3f        0       # axis ratio (b/a)" % (1.-ell)  
		print "10) %.3f         0       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA)
		print "Z) 0                  # leave in [1] or subtract [0] this comp from data"
		print "\n================================================================================"          

	sys.stdout = tmp_out
	f.close()
	os.chmod(r"modelPSF.txt",0777)
	subprocess.call("galfit modelPSF.txt", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)







def Header(file_image,file_out,bad_pix_mask,xc,yc,m0,pix2secx,pix2secy,RR):
	print "\n==============================================================================="
	print "# IMAGE and GALFIT CONTROL PARAMETERS"
	print "A) %s                # Input data image (FITS file)" % (file_image)
	print "B) %s  	              # Output data image block" % (file_out)
	print "C) none                # Sigma image name (made from data if blank or none)" 
	print "D) none                # Input PSF image and (optional) diffusion kernel"
	print "E) 1                   # PSF fine sampling factor relative to data" 
	print "F) %s                # Bad pixel mask (FITS image or ASCII coord list)" % (bad_pix_mask)
	print "G) none                # File with parameter constraints (ASCII file)" 
	print "H) %i    %i   %i    %i   # Image region to fit (xmin xmax ymin ymax)" % (xc-int(floor(RR)),xc+int(floor(RR)),yc-int(floor(RR)),yc+int(floor(RR)))
	print "I) 0    0        # Size of the convolution box (x y)"
	print "J) %.3f              # Magnitude photometric zeropoint" % (m0)
	print "K) %.3f    %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (pix2secx,pix2secy)
	print "O) regular             # Display type (regular, curses, both)"
	print "P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n"

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

def moffat (xc,yc,mtot,fwhm,ell,PA):
	 print "Moffat function"
 	 print "0) moffat           # object type"
	 print "1) %.3f  %.3f  1 1  # position x, y        [pixel]" % (xc,yc)
	 print "3) %.3f       1       # total magnitude" % (mtot)     
	 print "4) %.3f       1       #   FWHM               [pixels]" % (fwhm)
	 print "5) %.3f        1       # powerlaw" % (beta) 
	 print "9) %.3f      1       # axis ratio (b/a)" % (1.-ell)
	 print "10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA)
	 print "Z) 0                  # leave in [1] or subtract [0] this comp from data?"

def gauss (xc,yc,mtot,fwhm,ell,PA):
	 print "Gaussian function"
 	 print "0) gaussian           # object type"
	 print "1) %.3f  %.3f  1 1  # position x, y        [pixel]" % (xc,yc)
	 print "3) %.3f       1       # total magnitude" % (mtot)     
	 print "4) %.3f       1       #   FWHM               [pixels]" % (fwhm)
	 print "9) %.3f      1       # axis ratio (b/a)" % (1.-ell)
	 print "10) %.3f         1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA)
	 print "Z) 0                  # leave in [1] or subtract [0] this comp from data?"


def PSF_model_sph(image_file,reg_file,pix2sec,factor,mode):
	## Open region file with stars:
	f = open(reg_file, "r") 
	lines = f.readlines()
	
	## Read in the input image
	hdulist = pyfits.open(image_file)
	header1 = hdulist[0].header
	scidata = hdulist[0].data
	
	## Re-calculate m0 (for GALFIT) if needed
	if 'EXPTIME' in header1:
	  EXPTIME = float(header1['EXPTIME'])
	  m0 = 22.5-2.5*log10(EXPTIME)
	else:
	  m0 = 22.5
	#m0 = m0 - 5.*log10(factor)	#### ???
	  
	mtot_star=[]
	fwhm_star=[]
	beta_star=[]
	theta_star=[]
	Ell_star=[]
	chi2_star=[]
	fwhm_gaus=[]
	xc_star1 = []
	yc_star1 = []

	xcstars = []
	ycstars = []
	a_stars = []
	Backgr = []
	## Read in the PSF stars from the region file
	for k in range(len(lines)):
	    if 'circle(' in lines[k] or 'ellipse(' in lines[k]:
	      star = lines[k].split(',')
	      xcstars.append(float(star[0].split('(')[1]))
	      ycstars.append(float(star[1]))
	      try:
		a_stars.append(float(star[2]))
	      except:
		a_stars.append(float(star[2].split(')')[0]))

	## Masking contaminants except for the PSF stars
	bad_pix_mask(image_file,'segm.fits',"bad_pix.fits",xcstars,ycstars,a_stars,'fits')	####
	bad_pix_mask_fits = "bad_pix.fits"	####

	## Fitting each star to the gaussian or moffat function
	for k in range(len(lines)):

	  if 'circle(' in lines[k] or 'ellipse(' in lines[k] or 'annulus(' in lines[k]:
	    star = lines[k].split(',')
	    xc_star = float(star[0].split('(')[1])
	    yc_star = float(star[1])
	    if 'circle(' in lines[k] or 'ellipse(' in lines[k]:
	      try:
		if '{|' in lines[k]:
		  a_star = float(star[2])
		  b_star = float(star[3])
		  PA_star = float(star[4].split(')')[0])
		  m_star = float(star[5])
		  fwhm_st = float(star[6])
		  backgr_level = float(star[7].split('}')[0])
		else:
		  a_star = float(star[2])
		  b_star = float(star[3])
		  PA_star = float(star[4].split(')')[0])
		  m_star = float(star[4].split('{')[1])
		  fwhm_st = float(star[5])
		  backgr_level = float(star[6].split('}')[0])
	      except:
		  a_star = float(star[2].split(')')[0])
		  b_star = a_star
		  PA_star = 0.
		  apertures = CircularAperture([(xc_star,yc_star)], r = a_star)
		  phot_table = aperture_photometry(scidata, apertures)
		  Flux_circ = float(phot_table['aperture_sum'])
		  m_star = 22.5 - 2.5*log10(Flux_circ)
		  fwhm_st = a_star/2.
		  backgr_level,backgr_std = backgr_around_star(scidata,xc_star,yc_star,a_star)
	    if 'annulus(' in lines[k]:
		star = lines[k].split(',')
		xc_star = float(star[0].split('(')[1])
		yc_star = float(star[1])
		r_in = float(star[2])
		r_out = float(star[3].split(')')[0])
		a_star = r_in
		b_star = a_star
		PA_star = 0.
		apertures = CircularAperture([(xc_star,yc_star)], r = a_star)
		phot_table = aperture_photometry(scidata, apertures)
		Flux_circ = float(phot_table['aperture_sum'])
		m_star = 22.5 - 2.5*log10(Flux_circ)
		fwhm_st = a_star/2.
		backgr_level,backgr_std = backgr_around_star(scidata,xc_star,yc_star,a_star,r_out)
	  
	    fwhm_gaus.append(fwhm_st)

	    ## Sky subtraction:
	    sky_subtration(image_file,'sky_subtr.fits',backgr_level)

	      
	    ## Creating GALFIT input file:
	    try:
		f = open(r"modelIN.txt", "w") 
		sys.stdout = f
		Header('sky_subtr.fits','psf.fits',bad_pix_mask_fits,xc_star,yc_star,m0,pix2sec,pix2sec,max([1.2*a_star,1.2*b_star]))
		if window=='moffat':	moffat (xc_star,yc_star,m_star,fwhm_st,1.0-b_star/a_star,PA_star)
		if window=='gauss':	gauss (xc_star,yc_star,m_star,fwhm_st,1.0-b_star/a_star,PA_star)
		sys.stdout = tmp_out
		f.close()
		os.chmod(r"modelIN.txt",0777)

		## GALFIT RUNING
		subprocess.call("galfit modelIN.txt", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
		
		## Reading the results
		hdulist = pyfits.open('psf.fits')
		prihdr = hdulist[2].header
		xc_star1.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_XC']))[0])
		yc_star1.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_YC']))[0])
		mtot_star.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_MAG']))[0])
		fwhm_star.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_FWHM']))[0])
		if window=='moffat':
		  beta_star.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_C']))[0])
		else:
		  beta_star.append(float('nan'))
		theta_star.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_PA']))[0])
		Ell_star.append(1.-map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_AR']))[0])
		chi2_star.append(prihdr['CHI2NU'])
		Backgr.append(backgr_level)
		os.remove(r"modelIN.txt") 

		os.remove(r"fit.log") 
		os.remove(r"galfit.01")
		#try:
		#  os.remove(r"bad_pix.fits")
		#except:
		#    zz=1
		os.remove(r"psf.fits")
	    except:
		xc_star1.append(float('nan'))
		yc_star1.append(float('nan'))
		mtot_star.append(float('nan'))
		fwhm_star.append(float('nan'))
		beta_star.append(float('nan'))
		theta_star.append(float('nan'))
		Ell_star.append(float('nan'))
		chi2_star.append(float('nan'))
		Backgr.append(float('nan'))

	mtot_star = np.array(mtot_star)
	theta_star = np.array(theta_star)
	Ell_star = np.array(Ell_star)
	fwhm_star = np.array(fwhm_star)
	xc_star1 = np.array(xc_star1)
	yc_star1 = np.array(yc_star1)
	chi2_star = np.array(chi2_star)
	beta_star = np.array(beta_star)
	Backgr = np.array(Backgr)
	
	good_fits = np.where(np.isnan(mtot_star)==False)[0]
	

	
	## PRINT RESULTS
	#mtot_star = mtot_star[~np.isnan(mtot_star)]
	#theta_star = theta_star[~np.isnan(theta_star)]
	#Ell_star = Ell_star[~np.isnan(Ell_star)]
	#fwhm_star = fwhm_star[~np.isnan(fwhm_star)]
	#chi2_star = chi2_star[~np.isnan(chi2_star)]
	
	mtot_star = mtot_star[good_fits]
	theta_star = theta_star[good_fits]
	Ell_star = Ell_star[good_fits]
	fwhm_star = fwhm_star[good_fits]
	chi2_star = chi2_star[good_fits]
	xc_star1 = xc_star1[good_fits]
	yc_star1 = yc_star1[good_fits]
	beta_star = beta_star[good_fits]
	Backgr = Backgr[good_fits]
	
	try:
	  theta = np.mean(theta_star)
	  ell = np.mean(Ell_star)
	  fwhm = np.mean(fwhm_star)
	  mtot = np.mean(mtot_star)
	  CHI2 = np.mean(chi2_star)
	  if window=='moffat':
	    beta = np.median(beta_star)
	  else:
	    beta = float('inf')
	except:  
	  theta = np.mean(theta_star)
	  ell = np.mean(Ell_star)
	  fwhm = np.mean(fwhm_star)
	  mtot = np.mean(mtot_star)
	  CHI2 = np.mean(chi2_star)
	  if window=='moffat':
	    beta = np.mean(beta_star)
	  else:
	    beta = float('inf')

	chi2_star = list(chi2_star)
	num_min_chi2 = int(chi2_star.index(min(chi2_star)))    
	## Create image of the best model PSF
	if np.std(theta_star)>5.:
	    mod_psf(window,fwhm_star[num_min_chi2],mtot_star[num_min_chi2],0.,0.,22.5,beta_star[num_min_chi2],1.,1.)
	else:
	    mod_psf(window,fwhm_star[num_min_chi2],mtot_star[num_min_chi2],Ell_star[num_min_chi2],theta_star[num_min_chi2],22.5,beta_star[num_min_chi2],1.,1.)
	
	out_file = image_file.split('/')[-1].split('.fits')[0] + '.psf'
	fff = open(out_file,'w')
	print >> fff, '%s: Number of stars and resolution (arcsec) are: %i %.4f' % (image_file,len(Ell_star),pix2sec)
	print >> fff, 'SExtractor Gauss FWHM(pix):\t%.3f\t%.3f' % (mean(fwhm_gaus),std(fwhm_gaus))
	if window == 'moffat':
	  print >> fff, 'Moffat FWHM(pix),mag_total,ell,PA,beta:\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' % (fwhm,mtot,ell,theta,beta)
	  print >> fff, 'Moffat errors:\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f' % (std(fwhm_star),std(mtot_star),std(Ell_star),std(theta_star),std(beta_star))
	  print >> fff, 'CHI2:\t%.3f' % (CHI2)
	else:
	  print >> fff, 'Gauss FWHM(pix),mag_total,ell,PA:\t%.3f\t%.3f\t%.3f\t%.3f' % (fwhm,mtot,ell,theta)
	  print >> fff, 'Gauss errors:\t%.3f\t%.3f\t%.3f\t%.3f' % (std(fwhm_star),std(mtot_star),std(Ell_star),std(theta_star))	  
	  print >> fff, 'CHI2:\t%.3f' % (CHI2)
	  
	print >> fff, '\nNo\tMag\tFWHM\tEll\tTheta\tBeta\tChi2\tx0\ty0\tBackgr'
	No = 0
	for k in range(len(fwhm_star)):
	  No = No + 1
	  if window == 'gauss':
	    beta = float('inf')
	  else:
	    beta = beta_star[k]
	  if k==num_min_chi2:
	    print >> fff, '%i BEST\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.5f' % (No,mtot_star[k],fwhm_star[k],Ell_star[k],theta_star[k],beta,chi2_star[k],xc_star1[k],yc_star1[k],Backgr[k])
	  else:
	    print >> fff, '%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.5f' % (No,mtot_star[k],fwhm_star[k],Ell_star[k],theta_star[k],beta,chi2_star[k],xc_star1[k],yc_star1[k],Backgr[k])
	fff.close()

	## Rebin the model image according to the factor:
	if factor!=1 and 'rot' not in image_file:
	  rebin_image.downsample('psf.fits', factor,'psf.fits',set_wcs=False)
	else:
	  print 'There is no need to rebin psf image!'
	
	if 'rot' not in image_file:
	  ## Rotate the PSF image to match rotated frames (if needed)
	  if os.path.exists('./reg/rot.reg'):
	    crop_psf_star.rotate_frame('psf.fits',set_wcs=False)
	    shutil.move('psf.fits',image_file.split('/')[-1].split('.fits')[0] +'_'+window+'.fits')
	    shutil.move(image_file.split('/')[-1].split('.fits')[0] +'_'+window+'.fits','./psf/'+image_file.split('/')[-1].split('.fits')[0] +'_'+window+'.fits')
	  else:
	    shutil.move('psf.fits','./psf/'+image_file.split('/')[-1].split('.fits')[0] +'_'+window+'.fits')
	else:
	  print 'There is no need to rotate psf image!'
	  
	  
	## Cropping the best PSF star:
	if mode=='a':
	    best_star_reg = image_file.split('/')[-1].split('.fits')[0] + "_best_star.reg"
	elif mode=='sa':
	    best_star_reg = image_file.split('/')[-1].split('.fits')[0] + "_best_star.reg"
	    xc_best = xc_star1[num_min_chi2]
	    yc_best = yc_star1[num_min_chi2]
	    R_best = int(6.*fwhm_star[num_min_chi2])
    
	    ## Create tmp region file:
	    fout = open(best_star_reg, "w")
	    fout.truncate(0)
	    fout.write("# Region file format: DS9 version 4.1\n")
	    fout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 ')
	    fout.write('normal roman" select=1 highlite=1 dash=0 fixed=0 ')
	    fout.write('edit=1 move=1 delete=1 include=1 source=1\n')
	    fout.write('image\n')
	    fout.write("circle(%1.1f,%1.1f,%1.1f)\n" % (xc_best, yc_best, R_best))
	    fout.close()
	    imp_psf.psf_manual(image_file, best_star_reg)
	else:
	    ## Create tmp region file:
	    best_star_reg = image_file.split('/')[-1].split('.fits')[0] + "_best_star.reg"
	    imp_psf.psf_manual(image_file, best_star_reg)

	output_file = crop_psf_star.crop_star(image_file,best_star_reg,star_numb=1,factor=factor)
	shutil.move(output_file,'./psf/'+output_file)
	shutil.move(out_file,'./psf/'+out_file)
	shutil.move(reg_file,'./psf/'+reg_file)
	shutil.move(best_star_reg,'./psf/'+best_star_reg)
	


def main(image_input,reg_input,pix2sec,factor,mode):
  PSF_model_sph(image_input,reg_input,pix2sec,factor,mode)
  
'''
import imp_rebin
pix2sec,note = imp_rebin.resolution('K_final.fits')
print pix2sec
PSF_model_sph('K_final.fits','K_stars.reg',float(pix2sec),1.,'m')
'''
