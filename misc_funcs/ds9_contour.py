#! /usr/bin/env python
import pylab
import random
from random import gauss
import sys
import os
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
tmp_out = sys.stdout
from scipy import stats
import re
from scipy.optimize import fsolve
import pyfits
import subprocess
from scipy.odr.odrpack import *
import argparse
import heapq
from scipy.interpolate import splprep, splev
from scipy.interpolate import interp1d

from photutils import datasets
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import EllipticalAperture
import shutil
import glob
from scipy import interpolate

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
sys.path.append(PATH_TO_PACKAGE+'/Sextractor')


import polygon_reg
import superell_fit
import sextr
import read_sextr_res

warnings.filterwarnings("ignore")

# OPTIONS:
do_interpolation = 'yes'	#### Choose if you want to do intrpolation of the image, when the masked regions will be filled with some values 
repeat_interp = 'yes'	#### In case if the interpolated image has been created and you do not want to repeat interpolation again
simple_interp = 'no'	#### If 'yes' then convolution will be done to fill the masked pixels



# Colors to highlight the output text
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

def concatenate_broken_regions(x_1,y_1,x_2,y_2):
    # Find out which of the polygons is left and right or top and bottom:
    xc_1 = mean(x_1)
    xc_2 = mean(x_2)
    
    yc_1 = mean(y_1)
    yc_2 = mean(y_2)    
    
    if fabs(xc_1-xc_2)>fabs(yc_1-yc_2):
	# => There are left and right broken polygons:
	x_r = []; y_r = []; x_l = []; y_l = []
	
	if xc_1>xc_2:
	  x_r = x_1
	  y_r = y_1
	  
	  x_l = x_2
	  y_l = y_2
	else:
	  x_r = x_2
	  y_r = y_2
	  
	  x_l = x_1
	  y_l = y_1
	  
	# Concatanate
	x = np.concatenate((x_l, x_r))
	y = np.concatenate((y_l, y_r))

    else:
	# => There are top and bottom broken polygons:
	x_t = []; y_t = []; x_b = []; y_b = []
	
	if yc_1>yc_2:
	  x_t = x_1
	  y_t = y_1
	  
	  x_b = x_2
	  y_b = y_2
	else:
	  x_t = x_2
	  y_t = y_2
	  
	  x_b = x_1
	  y_b = y_1
	
	# Top:
	x1 = x_t[0:np.argmin(x_t)]
	y1 = y_t[0:np.argmin(x_t)]
	
	# Bottom:
	x2 = x_b[np.argmin(x_b):-1]
	y2 = y_b[np.argmin(x_b):-1]
	
	x3 = x_b[0:np.argmax(x_b)]
	y3 = y_b[0:np.argmax(x_b)]
	
	# Top:
	x4 = x_t[np.argmax(x_t):-1]
	y4 = y_t[np.argmax(x_t):-1]
	
	# Concatanate
	x = np.concatenate((x1, x2, x3, x4))
	y = np.concatenate((y1, y2, y3, y4))      

    return x,y
  
def replace_nans(array, max_iter, tol, kernel_size=1, method='localmean'):
    # Initialize arrays
    filled = np.empty( [array.shape[0], array.shape[1]], dtype=np.float64)
    kernel = np.empty( (2*kernel_size+1, 2*kernel_size+1), dtype=np.float64)

    # Indices where array is NaN
    inans, jnans = np.nonzero( np.isnan(array) )
    
    # Number of NaN elements
    n_nans = len(inans)
    
    # Arrays which contain replaced values to check for convergence
    replaced_new = np.zeros( n_nans, dtype=np.float64)
    replaced_old = np.zeros( n_nans, dtype=np.float64)

    # Depending on kernel type, fill kernel array
    if method == 'localmean':
        for i in range(2*kernel_size+1):
            for j in range(2*kernel_size+1):
                kernel[i,j] = 1.

    elif method == 'idw':
        kernel = np.array([[0, 0.5, 0.5, 0.5,0],
                  [0.5,0.75,0.75,0.75,0.5],
                  [0.5,0.75,1,0.75,0.5],
                  [0.5,0.75,0.75,0.5,1],
                  [0, 0.5, 0.5 ,0.5 ,0]])
        #print kernel, 'kernel'

    else:
        raise ValueError( 'method not valid. Should be one of `localmean`.')
    
    # Fill new array with input elements
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            filled[i,j] = array[i,j]

    # Make several passes
    # until we reach convergence
    for it in range(max_iter):
        
        # for each NaN element
        for k in range(n_nans):
            i = inans[k]
            j = jnans[k]
            
            # Initialize to zero
            filled[i,j] = 0.0
            n = 0
            
            # Loop over the kernel
            for I in range(2*kernel_size+1):
                for J in range(2*kernel_size+1):
                   
                    # if we are not out of the boundaries
                    if i+I-kernel_size < array.shape[0] and i+I-kernel_size >= 0:
                        if j+J-kernel_size < array.shape[1] and j+J-kernel_size >= 0:
                                                
                            # if the neighbour element is not NaN itself.
                            if filled[i+I-kernel_size, j+J-kernel_size] == filled[i+I-kernel_size, j+J-kernel_size] :
                                
                                # do not sum itself
                                if I-kernel_size != 0 and J-kernel_size != 0:
                                    
                                    # convolve kernel with original array
                                    filled[i,j] = filled[i,j] + filled[i+I-kernel_size, j+J-kernel_size]*kernel[I, J]
                                    n = n + 1

            # Divide value by effective number of added elements
            if n != 0:
                filled[i,j] = filled[i,j] / n
                replaced_new[k] = filled[i,j]
            else:
                filled[i,j] = 0.0
                
        # Check if mean square difference between values of replaced
        # elements is below a certain tolerance
        if np.mean( (replaced_new-replaced_old)**2 ) < tol:
            break
        else:
            for l in range(n_nans):
                replaced_old[l] = replaced_new[l]
    
    return filled

def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])

def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])

def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))



def main(input_image,mask_image,mag_level,m0,pix2sec,FWHM,isophote_file=None,what_to_do='all'):
	'''
	INPUT:
	input_image - image with the centered galaxy
	mask_image - image with masked objects(Intensity>0)
	mag_level - level of the outer isophote in [mag/arcsec**2]
	pix2sec - pixel scale in [arcsec/pixel]
	m0 - magnitude zero point in [mag/arcsec**2]
	FWHM - PSF FWHM in [pixels]
	isophote_file - region file with the outer isophote. In this case mag_level is not important, and you can give it any value.
	
	OUTPUT:
	For the ellipse describing the isophote: 
	  PARS_OUT['Ellipse_flux'] = [Flux with the masked objects [DN],Flux with the interpolated masked objects [DN]]
	  PARS_OUT['Isophote'] = [x-ccordinate of the center [pix],z-coordinate of the center [pix],semimajor axis [pix],semiminor axis [pix],angle of the major axis [deg]]
	For the isophote:
	  PARS_OUT['Iso_masked'] = [Total flux inside the isophote with masked objects [DN],mean SB [DN],median SB [DN],std [DN]]
	  PARS_OUT['Iso_interp'] = [Total flux with interpolated masked objects,mean SB [DN],median SB [DN],std [DN]]
	  PARS_OUT['Total_image_flux'] = [Total flux integrated over the whole image except for the masked objects [DN]]
	region file - region file with the polygon describing outer isophote
	error - perhaps, the outer isophote cannot cover the galaxy, i.e. the built polygon is wrong
	''' 
	  
	print '\t\t\t',bcolors.OKBLUE+ "*****************************************" + bcolors.ENDC
	print '\t\t\t',bcolors.OKBLUE+ "         APERTURE PHOTOMETRY 2015" + bcolors.ENDC
	print '\t\t\t',bcolors.OKBLUE+ "*****************************************" + bcolors.ENDC
	
	# Output parameters:
	PARS_OUT = {}
	PARS_OUT['Ellipse_flux'] = [float(nan),float(nan)]
	PARS_OUT['Iso_masked'] = [float(nan),float(nan),float(nan),float(nan)]
	PARS_OUT['Iso_interp'] = [float(nan),float(nan),float(nan),float(nan)]
	PARS_OUT['Isophote'] = [float(nan),float(nan),float(nan),float(nan),float(nan)]
	PARS_OUT['Total_image_flux'] = [float(nan)]
	
	error = 0	# Output value to tell that the program crashed somewhere, or the results are unrelaible. 0 is OK.

	hdulist = pyfits.open(input_image)
	inframe = hdulist[0].data
	prihdr = hdulist[0].header
	nx, ny =inframe.shape[1], inframe.shape[0]

	nan_frame = np.copy(inframe)
	zero_frame = np.copy(inframe)
	
	
	if mask_image==None:
	  mask_infr = np.zeros_like(inframe)
	  hdu = pyfits.PrimaryHDU(mask_infr)
	  hdu.writeto('zero_mask.fits',clobber=True)
	  mask_image = 'zero_mask.fits'
	else:
	  hdulist_infr = pyfits.open(mask_image)
	  mask_infr = hdulist_infr[0].data

	mask_astropy = np.zeros_like(inframe, dtype=bool)

	for k in range(ny):
	  for i in range(nx):
	    if mask_infr[k,i]>0:
	      nan_frame[k,i] = float('nan')
	      zero_frame[k,i] = 0.
	      mask_astropy[k,i] = True

	# INTERPOLATION
	#if not os.path.exists('interp.fits'):
	if do_interpolation=='yes':
	  if os.path.exists('interp.fits') and repeat_interp=='no':
	    hdulist_i = pyfits.open('interp.fits')
	    interpolated = hdulist_i[0].data  
	    print "WARNING: Using existing interp.fits file!"    
    
	  else:  
	    if simple_interp=='yes':
	      interpolated = replace_nans(nan_frame, 5, 0.5, int(ceil(FWHM)), "localmean")
	      hdu = pyfits.PrimaryHDU(interpolated,prihdr)
	      hdu.writeto('interp.fits',clobber=True)
	    else:
	      hdu = pyfits.PrimaryHDU(zero_frame,prihdr)
	      hdu.writeto('zero.fits',clobber=True)
	      sextr.run_sextr('zero.fits',m0,4.0,pix2sec,FWHM,'mask_contam.sex')
	      os.rename('models.fits','model_galaxy.fits')
	      hdulist_obj = pyfits.open('model_galaxy.fits')
	      mask_obj = hdulist_obj[0].data
	      hdulist_obj.close
	      
	      
	      sextr.run_sextr(input_image,m0,4.0,pix2sec,FWHM,'mask_contam.sex')
	      '''
	      try:
		os.remove("cont.fits")
	      except:
		zz = 1
	      try:
		os.remove("back.fits")
	      except:
		zz = 1      
	      subprocess.call("imarith models.fits model_galaxy.fits sub cont.fits", shell=True)
	      subprocess.call("imarith %s cont.fits sub back.fits" % (input_image), shell=True)
	      hdulist_back = pyfits.open('back.fits')
	      mask_back = hdulist_back[0].data    
	      
	      hdulist_interp = pyfits.open('interp.fits')
	      mask_interp = hdulist_interp[0].data   
	      
	      inframe1 = np.copy(inframe)

	      hdulist_mod = pyfits.open('models.fits')
	      models_inf = hdulist_mod[0].data
	      '''
	      hdulist_rms = pyfits.open('rms.fits')
	      rms = hdulist_rms[0].data
	      
	      
	      int_level = 10**(0.4*(m0-(mag_level-2.)+5.*log10(pix2sec)))

	      for k in range(ny):
		for i in range(nx):
		  if mask_infr[k,i]>0:
		    if mask_obj[k,i]>1.*rms[0,0]:        #### NOTE: Change rms level (by default, 1.)              
			inframe[k,i] = mask_obj[k,i] #(mask_back[k,i]+mask_interp[k,i]+mask_obj[k,i])/3.
		    else:
			inframe[k,i] = gauss(0.0, rms[0,0])#(mask_back[k,i]+mask_interp[k,i]+mask_obj[k,i])/3. #gauss(0.0, rms[0,0]) #(mask_back[k,i]+mask_interp[k,i]+mask_obj[k,i])/3.  


	      hdu = pyfits.PrimaryHDU(inframe,prihdr)
	      hdu.writeto('interp.fits',clobber=True)
	else:
          shutil.copy(input_image,'interp.fits')
	  hdulist_i = pyfits.open('interp.fits')
	  interpolated = hdulist_i[0].data  
	  print "WARNING: Interpolation was not done. Just coppied the input file to interp.fits!"

        #subprocess.call('find . -type f -not -name \'interp.fits\' -not -name \'%s\' | xargs rm' % (input_image), shell=True)
	method_fit = 'ds9'
	if method_fit=='ds9':
		  
		  # ISOPHOTE PLOTING
		  if isophote_file==None:
		    fix_iso = 0
		    int_level = 10**(0.4*(m0-mag_level+5.*log10(pix2sec)))
		    print bcolors.OKGREEN+"Intensity level corresponding to the given magnitude of %.1f mag arcsec^-2" % (mag_level) + bcolors.ENDC+": %.3f DN" % (int_level)
		    #exit()
		    subprocess.call("ds9 %s -scale log -cmap b -contour yes -contour limits %.1f %.1f -contour smooth 10 -contour nlevels 1 -contour convert -regions save con.reg -exit" % ('interp.fits',int_level,int_level), shell=True)
		    reg = open('con.reg','r')

		    lines = reg.readlines()
		    line = heapq.nlargest(2, list(lines), key=len)[0]

		    isophote_file = 'isophote_' + str(mag_level) + '.reg'
		    iso = open(isophote_file,'w')
		    print >>iso, line

		    os.remove('con.reg')
		    iso.close()
		  else:
		    print "WARNING: Using given region file with %s isophote!" % (str(mag_level))
		    fix_iso = 1


		  # READING THE POLYGON OF THE ISOPHOTE FROM THE FILE
		  reg = open(isophote_file,'r')
		  x = []
		  y = []
		  for Line in reg:
		    if 'polygon' in Line:
		      coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]
		      for k in range(0,len(coords)-1,2):
			x.append(float(coords[k]))
			y.append(float(coords[k+1]))
			
		  x = np.array(x)
		  y = np.array(y)
		  

		  # CHECKING IF THE ISOPHOTE IS CROPPED
		  xc_mean = mean(x)
		  yc_mean = mean(y)
		  BREAK=1
		  for kkk in range(len(x)-1):
		      if fabs(x[kkk] - x[kkk+1])>xc_mean:
			  BREAK=0

		  if fix_iso==0:            
		    if fabs(y[0]-y[-1])>yc_mean or BREAK==0:
			  print 'Perhaps, isophote is toren open!'
			  # Trying to concatenate two biggest polygons:
			  x1 = x 
			  y1 = y
					  
			  line_2 = heapq.nlargest(2, list(lines), key=len)[1]

			  x2 = []; y2 = []
			  coords = line_2.replace('(',',').replace(')',',').split(',')[1:-1]
			  for k in range(0,len(coords)-1,2):
			      x2.append(float(coords[k]))
			      y2.append(float(coords[k+1]))
			      
			  x2 = np.array(x2)
			  y2 = np.array(y2)
			  try:
			      x,y = concatenate_broken_regions(x1,y1,x2,y2)
			      isophote_file = 'isophote_' + str(mag_level) + '.reg'
			      iso = open(isophote_file,'w')
			      iso.write('polygon(')
			      for k in range(len(x)):
				  if k!=len(x)-1:
				      iso.write('%.3f,%.3f,' % (x[k],y[k]))
				  else:
				      iso.write('%.3f,%.3f' % (x[k],y[k]))
			      iso.write(')')
			      iso.close()
			      error = 2
			  except:
			      error = 3
			      
			      
		  # FITTING THE ISOPHOTE WITH ELLIPSE
		  try:
		    a = superell_fit.fitEllipse(x,y)
		    center = ellipse_center(a)
		    axes = ellipse_axis_length(a)
		    phi = degrees(ellipse_angle_of_rotation(a))
		    xc = center[0]
		    yc = center[1]
		    sma = axes[0]
		    smb = axes[1]
		  except:
		    sma = float(nan)
		    smb = float(nan)
		    xc = float(nan)
		    yc = float(nan)

			      
		  '''  
		  Additional interpolation - slow!
		  if fix_iso==0:
		    #x = np.round(x)
		    #y = np.round(y)
		    I = np.arange(len(x))

		    # 5x the original number of points
		    interp_i = np.linspace(0, I.max(), 5 * I.max())

		    X = interp1d(I, x, kind='cubic')(interp_i)
		    Y = interp1d(I, y, kind='cubic')(interp_i)


		    x = X
		    y = Y
		    
		    isophote_file = 'isophote.reg'
		    iso = open(isophote_file,'w')
		    iso.write('polygon(')
		    for k in range(len(x)):
		      if k!=len(x)-1:
			iso.write('%.3f,%.3f,' % (x[k],y[k]))
		      else:
			iso.write('%.3f,%.3f' % (x[k],y[k]))
		    iso.write(')')
		    iso.close()  	
		  '''
		  xc_iso = xc
		  yc_iso = yc
		  sma_iso = sma
		  smb_iso = smb
		  
		  if xc==float(nan) or yc==float(nan) or sma==float(nan) or smb==float(nan):
		    print "WARNING: Ellipse fitting failed! Using rough estimation of isophote shape."
		    xc_iso = (max(x) - min(x)) / 2. + min(x)
		    yc_iso = (max(y) - min(y)) / 2. + min(y)
		    sma_iso = max(fabs(x-xc_iso))
		    smb_iso = max(fabs(y-yc_iso))

	elif method_fit=='iraf':
	      import iraf_fit_ellipse
	      iraf_fit_ellipse.main('interp.fits',nx/2.,ny/2.,m0,pix2sec,minsma=1.,step=1.)
	      arr = glob.glob('*_ell.txt')
	      sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = loadtxt(arr[0], usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 6, dtype='str')
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
	      
	      tck = interpolate.splrep(sma*pix2sec,inten,s=0)
	      r_interp = np.arange(min(sma*pix2sec),max(sma*pix2sec),0.5)

	      I_interp = interpolate.splev(r_interp,tck,der=0)
	      mag_interp = m0 - 2.5*log10(I_interp) + 5.*log10(pix2sec)
	      
	      def find_nearest(array,value):
		  idx = (np.abs(array-value)).argmin()
		  return idx
	      
	      idx = find_nearest(mag_interp,mag_level)
	      r_mag_level = r_interp[idx] / pix2sec	# in pix
	      
	      idx = find_nearest(sma, r_mag_level)
	      
	      sma_iso = sma[idx] * ell[idx] 
	      smb_iso = sma_iso*(1.-ell[idx])
	      phi = PA[idx]-90.
	      xc_iso = x0[idx]
	      yc_iso = y0[idx]
	elif method_fit=='sextractor':
	      sextr.run_sextr('interp.fits',m0,4.0,pix2sec,FWHM,)
	      xc,yc,RA_c,DEC_c,PA,a,b,flu_rad,ma_auto,ell,kron_r,C31,rmax = read_sextr_res.read_sextr_find(nx/2.,ny/2.,'field.cat',coords=1)
	      #print xc,yc,RA_c,DEC_c,PA,a,b,flu_rad,ma_auto,ell,kron_r,C31,rmax
	      sma_iso = a*kron_r
	      smb_iso = b*kron_r
	      phi = PA
	      xc_iso = xc
	      yc_iso = yc

	      
	try:
	  print bcolors.OKGREEN+"Estimations of semimajor- and semiminor-axes (pix and arcsec): "+ bcolors.ENDC+ "%.2f, %.2f\t%.2f, %.2f" % (sma_iso,smb_iso,sma_iso*pix2sec,smb_iso*pix2sec)

	  apertures = EllipticalAperture([(xc_iso,yc_iso)], sma_iso,smb_iso,radians(phi))
	  phot_table = aperture_photometry(inframe, apertures,mask=mask_astropy)

	  Flux_ell_masked = float(phot_table['aperture_sum'])

	  phot_table = aperture_photometry(inframe, apertures)
	  Flux_ell_interp = float(phot_table['aperture_sum'])

	  print bcolors.OKGREEN+'Total flux within the ellipse (with applied masked): '+ bcolors.ENDC+ ' %.3f DN,\t %.3f mag' % (Flux_ell_masked,-2.5*log10(Flux_ell_masked)+m0)
	  print bcolors.OKGREEN+'Total flux within the ellipse (no mask, interpolated): '+ bcolors.ENDC+ ' %.3f DN,\t %.3f mag' % (Flux_ell_interp,-2.5*log10(Flux_ell_interp)+m0)
	except:
	  zz = 1
	  
	try:
	  Flux_iso_masked,mean_iso_masked,median_iso_masked,std_iso_masked = polygon_reg.stats(input_image,isophote_file,mask_image=mask_image)
	  print bcolors.OKGREEN+'Flux, mean, median and std.dev inside the given isophote (with mask): ' + bcolors.ENDC+ '%.3f DN\t %.3f mag\t %.3f DN\t%.3f DN\t%.3f DN' % (Flux_iso_masked,-2.5*log10(Flux_iso_masked)+m0,mean_iso_masked,median_iso_masked,std_iso_masked)
	except:
	  zz = 1
	  
	try:
	  Flux_iso_interp,mean_iso_interp,median_iso_interp,std_iso_interp = polygon_reg.stats('interp.fits',isophote_file)
	  print bcolors.OKGREEN+'Flux, mean, median and std.dev inside the given isophote (interpolated): '+ bcolors.ENDC+  '%.3f DN\t %.3f mag\t %.3f DN\t%.3f DN\t%.3f DN' % (Flux_iso_interp,-2.5*log10(Flux_iso_interp)+m0,mean_iso_interp,median_iso_interp,std_iso_interp)
	except:
	  zz = 1
	 
	try:  
	  PARS_OUT['Ellipse_flux'] = [Flux_ell_masked,Flux_ell_interp]
	except:
	  zz = 1

	try:  
	  PARS_OUT['Isophote'] = [xc_iso,yc_iso,sma_iso,smb_iso,phi]
	except:
	  zz = 1
  
	try:  
	  PARS_OUT['Iso_masked'] = [Flux_iso_masked,mean_iso_masked,median_iso_masked,std_iso_masked]
	except:
	  zz = 1
	  
	try:  
	  PARS_OUT['Iso_interp'] = [Flux_iso_interp,mean_iso_interp,median_iso_interp,std_iso_interp]
	except:
	  zz = 1
	
	TOTAL = 0.
	for k in range(ny):
	  for i in range(nx):
	    if mask_infr[k,i]==0.:
	      TOTAL = TOTAL + inframe[k,i]

	PARS_OUT['Total_image_flux'] = TOTAL
	print bcolors.OKGREEN+'Total flux of the whole image with the mask applied: '+ bcolors.ENDC+  '%.3f DN\t %.3f mag' % (TOTAL,-2.5*log10(TOTAL)+m0)
	
	try:
	  # PLOTTING THE IMAGE WITH THE ISOPHOTE
	  isophote_image = 'ds9_contour_' + str(mag_level) +'.png'
	  zoom = 1.
	  input_image = 'interp.fits'
	  max_i = inframe[int(yc_iso),int(xc_iso)]
	  command1 = "ds9 -geometry %ix%i -zoom %.2f " % (1.01*nx*zoom,260+ny*zoom,zoom)
	  command1 = command1 + " %s -scale log -scale limits 0 %.1f -invert " % (input_image,max_i)
	  ell = "-regions command \"ellipse %f %f %f %f %f # color=blue width=2\" " % (xc_iso,yc_iso,sma_iso,smb_iso,phi)
	  point = "-regions command \"point %f %f # point=x color=blue width=1\" " % (xc_iso,yc_iso)
	  
	  if method_fit=='ds9':
	    command1 = command1 + point + ell + "-regions load %s -colorbar no -saveimage png %s -exit" % (isophote_file,isophote_image)
	  else:
	    command1 = command1 + point + ell + "-colorbar no -saveimage png %s -exit" % (isophote_image)
	  subprocess.call(command1, shell=True)
	except:
	  isophote_file = None
	  isophote_image = None
	  error = 4
	try:
	  if Flux_ell_masked<2.*TOTAL/3.:
	    print bcolors.FAIL+ "Perhaps, your image is cropped at the input isophote! " + bcolors.ENDC
	    error = 1
	except:
	  zz=1
	try:
	  os.remove('zero_mask.fits')
	except:
	  zz = 1
	
	if do_interpolation!='yes':
	  os.remove('interp.fits')
	print 'Done!'
	return PARS_OUT,isophote_file,isophote_image,error


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Isophote analysis")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("maskImage", nargs='?', const=1, help="Input mask image (with masked objects which have intensities >0",type=str,default=None)
    parser.add_argument("magLevel", help="Input the outer isophote level in [mag/arcsec^2]")
    parser.add_argument("ZeroPoint", help="Input Zero Point in [mag/arcsec^2]") 
    parser.add_argument("Scale", help="Input scale in [arcsec/pix]")
    parser.add_argument("FWHM", help="Input PSF FWHM in [pixels]") 
    parser.add_argument("regFile", nargs='?', const=1, help="Optional: Input region file with the outer isophote. In this case mag_level is not important, and you can give it any value.",type=str,default=None) 
    parser.add_argument("what_to_do", nargs='?', const=1, help="Optional: Choose how much information you want to receive: all,.",type=str,default='all') 
    args = parser.parse_args()

    input_image = args.inputImage
    mask_image = args.maskImage
    mag_level = float(args.magLevel)
    m0 = float(args.ZeroPoint)
    pix2sec = float(args.Scale)
    FWHM = float(args.FWHM)
    regFile = args.regFile
    what_to_do = args.what_to_do

    main(input_image,mask_image,mag_level,m0,pix2sec,FWHM,isophote_file=regFile,what_to_do=what_to_do)      
    #main('cropped_r.fits','mask.fits',26.3,28.240,1.,3.25,0.7)
    
    