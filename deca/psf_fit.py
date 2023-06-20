#!/usr/bin/python
# -*- coding:  cp1251 -*-

import sys
import math
from pylab import *
import subprocess
import pyfits
import re
import os

import deca_setup

tmp_out = sys.stdout
FNULL = open(os.devnull, 'w')

def Header(file_image,file_out,bad_pix_mask,xc,yc,m0,pix2secx,pix2secy,nx,ny):
	print "\n==============================================================================="
	print "# IMAGE and GALFIT CONTROL PARAMETERS"
	print "A) %s                # Input data image (FITS file)" % (file_image)
	print "B) %s  	              # Output data image block" % (file_out)
	print "C) none                # Sigma image name (made from data if blank or none)" 
	print "D) none                # Input PSF image and (optional) diffusion kernel"
	print "E) 1                   # PSF fine sampling factor relative to data" 
	print "F) %s                # Bad pixel mask (FITS image or ASCII coord list)" % (bad_pix_mask)
	print "G) none                # File with parameter constraints (ASCII file)" 
	print "H) %i    %i   %i    %i   # Image region to fit (xmin xmax ymin ymax)" % (1,nx,1,ny)
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
	 print "5) %.3f        1       # powerlaw" % (3.) 
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





def main(psf_image, xc_star, yc_star, nx, ny, window):
		print 'Estimating FWHM...'
		f = open(r"modelIN.txt", "w") 
		sys.stdout = f
		Header(psf_image,'psf_tmp.fits','none',xc_star,yc_star,20.,1.,1.,nx,ny)
		if window=='moffat':	moffat (xc_star,yc_star,10.,5.,0.,0.)
		if window=='gauss':	gauss (xc_star,yc_star,10.,5.,0.,0.)
		sys.stdout = tmp_out
		f.close()
		os.chmod(r"modelIN.txt",0777)

		## GALFIT RUNING
		subprocess.call("%sgalfit modelIN.txt" % (deca_setup.galfitPath), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
		## Reading the results
		hdulist = pyfits.open('psf_tmp.fits')
		prihdr = hdulist[2].header
		#xc_star1.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_XC']))[0])
		#yc_star1.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_YC']))[0])
		#mtot_star.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_MAG']))[0])
		fwhm_star = map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_FWHM']))[0]
		'''
		if window=='moffat':
		  beta_star.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_C']))[0])
		else:
		  beta_star.append(float('nan'))
		'''
		#theta_star.append(map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_PA']))[0])
		#Ell_star.append(1.-map(float, re.findall(r"[-+]?\d*\.\d+|\d+", prihdr['1_AR']))[0])
		#chi2_star.append(prihdr['CHI2NU'])
		#Backgr.append(backgr_level)
		os.remove(r"modelIN.txt") 

		os.remove(r"fit.log") 
		os.remove(r"galfit.01")
		os.remove(r"psf_tmp.fits")
		print 'FWHM is %.2f pix' % (fwhm_star)
		return fwhm_star
