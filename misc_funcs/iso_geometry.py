#! /usr/bin/env python
import pylab
import random as random_number
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
from scipy import stats
import re
from scipy.optimize import fsolve
import pyfits
import subprocess
from scipy.odr.odrpack import *
import argparse
from photutils import aperture_photometry
from photutils import EllipticalAperture
from photutils import datasets
from astropy.stats import sigma_clipped_stats

import polygon_reg
import superell_fit

tmp_out = sys.stdout

sys.path.append('/home/amosenko/CurrentWork/ImaPrep')
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

def crea_reg(f,xc,zc,sma,smb,teta):
    print >>f, '# Region file format: DS9 version 4.1'
    print >>f, '# Filename: cropped_r.fits'
    print >>f, 'global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
    print >>f, 'image'
    print >>f, 'ellipse(%.1f,%.1f,%.1f,%.1f,%.1f) # color=black width=2  ' % (xc,zc,sma,smb,teta)


def bps_bulge(B, x):
    return B[1]* (1.-fabs(x/B[0])**(B[2]+2))**(1./(B[2]+2.))

def odr_func(x,y,a0,b0,c0):
  	dataToFit = RealData(x,y)
	mo = Model(bps_bulge)
	fitting = ODR(dataToFit, mo, [a0,b0,c0])
	fitting.set_job()
	fit_result = fitting.run()
	a1 = fit_result.beta[0]
	b1 = fit_result.beta[1]
	c1 = fit_result.beta[2]	
	return a1,b1,c1


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
                filled[i,j] = np.nan
                
        # Check if mean square difference between values of replaced
        # elements is below a certain tolerance
        if np.mean( (replaced_new-replaced_old)**2 ) < tol:
            break
        else:
            for l in range(n_nans):
                replaced_old[l] = replaced_new[l]
    
    return filled




def main(input_image,mag_level,R_disc_dom,R_bulge_dom,m0,pix2sec,xc=0.,zc=0.):
	'''
	Script to find disc and bulge isophote parameters by the outer isophote.
	INPUT:
	input_image - image with the centered galaxy.
	mag_level - level of the outer isophote in [mag/arcsec**2].
	R_disc_dom - disc dominated radius [pix].
	R_bulge_dom - bulge dominated radius [pix].
	pix2sec - pixel scale in [arcsec/pixel]
	m0 - magnitude zero point in [mag/arcsec**2]
	Optional:
	xc,zc - coordinates of the center.
	
	Output:
	PARS_OUT['Bulge_iso'] = [semimajor axis of the bulge [pix],semiminor axis of the bulge [pix],ellipse index]
	PARS_OUT['Disc_iso'] = [semimajor axis of the disc [pix],semiminor axis of the disc [pix],ellipse index]
	PARS_OUT['Galaxy_iso'] = [center x-coordinate [pix],center z-coordinate [pix],semimajor axis of the galaxy [pix],semiminor axis of the galaxy [pix],ellipse index]
	'''
 
	print '\t\t\t',bcolors.OKBLUE+ "*****************************************" + bcolors.ENDC
	print '\t\t\t',bcolors.OKBLUE+ "Geometry of the outer isophote (c) Mosenkov A. 2015" + bcolors.ENDC
	print '\t\t\t',bcolors.OKBLUE+ "*****************************************" + bcolors.ENDC
	
	# Output parameters:
	PARS_OUT = {}
	
	PARS_OUT['Disc_iso'] = [float(nan),float(nan),float(nan)]
	PARS_OUT['Galaxy_iso'] = [float(nan),float(nan),float(nan),float(nan),float(nan)]
	PARS_OUT['Bulge_iso'] = [float(nan),float(nan),float(nan)]
	

	hdu = pyfits.open(input_image)
	data = hdu[0].data
	ySize, xSize = data.shape


	if zc==0.:
	  zc=ySize/2.
	if xc==0.:
	  xc=xSize/2.

	int_level = 10**(0.4*(m0-mag_level+5.*log10(pix2sec)))
	print bcolors.OKGREEN+"Intensity level corresponding to the given magnitude" + bcolors.ENDC+": %.3f DN" % (int_level)

	subprocess.call("ds9 %s -scale log -cmap b -contour yes -contour limits %.1f %.1f -contour smooth 7 -contour nlevels 1 -contour convert -regions save con.reg -exit" % (input_image,int_level,int_level), shell=True)
	reg = open('con.reg','r')

	lines = reg.readlines()
	line = max(lines, key=len)

	iso = open('temp_iso.reg','w')
	print >>iso, lines[0],lines[1],'global color=blak\n',lines[3],line
	os.remove('con.reg')
	iso.close()

	reg = open('temp_iso.reg','r')

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
	

	#try:
	if True:
	  # Find geometry of the bulge
	  X = []
	  Y = []
	  for k in range(len(x)):
	      if fabs(x[k]-xc)<=R_bulge_dom:
		X.append(x[k])
		Y.append(y[k])

	  sma_b,smb_b,c0_b,xc_b,zc_b = superell_fit.main(X,Y,xc,zc,out=1)
	  print bcolors.OKGREEN+'Semimajor- and semi-minor axes of the bulge [pix]: ' + bcolors.ENDC+' %.3f\t%.3f' % (sma_b,smb_b) 

	#except:
	#  zz = 1
	  
	try:  
	  # Find geometry of the disc
	  X = []
	  Y = []
	  for k in range(len(x)):
	      if fabs(x[k]-xc)>R_disc_dom:
		X.append(x[k])
		Y.append(y[k])

	  sma_d,smb_d,c0_d,xc_d,zc_d = superell_fit.main(X,Y,xc,zc,out=1)
	  print bcolors.OKGREEN+'Semimajor- and semi-minor axes of the disc [pix]: ' + bcolors.ENDC+' %.3f\t%.3f' % (sma_d,smb_d) 

	except:
	  zz = 1
	  
	try:  
	  sma_g,smb_g,c0_g,xc_g,zc_g = superell_fit.main(x,y,xc,zc,out=1)
	  print bcolors.OKGREEN+'Semimajor- and semi-minor axes of the galaxy [pix]: ' + bcolors.ENDC+' %.3f\t%.3f' % (sma_g,smb_g) 
	  plt.show()
	  exit()
	except:
	  zz = 1	  
  
	try:  
	  PARS_OUT['Bulge_iso'] = [sma_b,smb_b,c0_b]
	except:
	  zz = 1
	  
	try:  
	  PARS_OUT['Disc_iso'] = [sma_d,smb_d,c0_d]
	except:
	  zz = 1  
	  
	try:  
	  PARS_OUT['Galaxy_iso'] = [xc_g,zc_g,sma_g,smb_g,c0_g]
	except:
	  zz = 1
	  
	#plt.plot(R,B/A,'o')
	#plt.plot(R,C0,'o',color='red')
	#plt.show()  
	iso = open(input_image.split('.fits')[0]+'_iso_components.reg','w')
	print >>iso, 'ellipse(%.3f,%.3f,%.3f,%.3f,%.3f) # color=green' % (xc_g,zc_g,sma_d,smb_d,0.)
	print >>iso, 'ellipse(%.3f,%.3f,%.3f,%.3f,%.3f) # color=blue' % (xc_g,zc_g,sma_b,smb_b,0.)
	print >>iso, 'ellipse(%.3f,%.3f,%.3f,%.3f,%.3f) # color=red' % (xc_g,zc_g,sma_g,smb_g,0.)
	iso.close()
	
	zoom = 1.
	max_i = data[int(zc_g),int(xc_g)]
	command1 = "ds9 -geometry %ix%i -zoom %.2f " % (1.01*xSize*zoom,260+ySize*zoom,zoom)
	command1 = command1 + " %s -scale log -scale limits 0 %.1f -invert " % (input_image,max_i)
	reg = "-regions load %s " % (input_image.split('.fits')[0]+'_iso_components.reg')
	
	command1 = command1 + reg + "-regions load %s -colorbar no -saveimage png %s -exit" % ('temp_iso.reg','geom_contour.png')
	subprocess.call(command1, shell=True)
	os.remove('temp_iso.reg')

	return PARS_OUT


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Isophote analysis")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("magLevel", help="Input the outer isophote level in [mag/arcsec^2]",type=float)
    parser.add_argument("radiusDisc", help="Input the radius at which disc dominates over the bulge [pix]")
    parser.add_argument("radiusBulge", help="Input the radius at which bulge still dominates over the disc [pix]")
    parser.add_argument("ZeroPoint", help="Input Zero Point in [mag/arcsec^2]",type=float) 
    parser.add_argument("Scale", help="Input scale in [arcsec/pix]",type=float)
    parser.add_argument("xc", nargs='?', const=0., help="Input x-coordinate of the galaxy center [pix]",type=float,default=0.)
    parser.add_argument("zc", nargs='?', const=0., help="Input z-coordinate of the galaxy center [pix]",type=float,default=0.)    
    
    args = parser.parse_args()
    
    input_image = args.inputImage
    mag_level = float(args.magLevel)
    R_disc_dom = float(args.radiusDisc)
    R_bulge_dom = float(args.radiusBulge)    
    m0 = float(args.ZeroPoint)
    pix2sec = float(args.Scale)
    xc = args.xc
    zc = args.zc

    main(input_image,mag_level,R_disc_dom,R_bulge_dom,m0,pix2sec,xc=xc,zc=zc)      
    #main('cropped_r.fits','mask.fits',26.3,28.240,1.,3.25,0.7)
    
    