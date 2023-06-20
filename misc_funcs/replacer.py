# THE SCRIPT TO RUN SEXTRACTOR
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
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
import pyfits
import argparse
from photutils import datasets
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import EllipticalAperture

tmp_out = sys.stdout


PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
sys.path.append(PATH_TO_PACKAGE+'/FindFluxes')
sys.path.append('/home/amosenko/CurrentWork/ImaPrep')

import mask_indiv
import ds9_contour
import sextr

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

def ds9_show(fileToShow,go_fur='no'):
    ds9Proc = subprocess.Popen(["ds9", fileToShow,
                                "-scale", "histequ","-cmap","Rainbow"])
    if go_fur=='no':
      ds9Proc.wait()


def stats(inframe,iso_region,mask_image='none',layer=0):
		iso_reg = open(iso_region,'r')
		paths=[]
		ySize, xSize = inframe.shape

		for line in iso_reg:
			# Detect polygon regions and add them to the path list
			if 'polygon(' in line:
				param = map(float, ((line.split('polygon(')[1]).split(')')[0]).split(','))
				param2 = [None]*(len(param)/2)
				for i in xrange(0,len(param2)):	param2[i] = (int(param[2*i]),int(param[2*i+1])) 
				param2.append(param2[0])
				codes = []
				codes.append(Path.MOVETO)
				for i in xrange(1,len(param2)-1): codes.append(Path.LINETO)
				codes.append(Path.CLOSEPOLY)
				path = Path(param2, codes)
				paths.append(path)


				nx, ny =inframe.shape[1], inframe.shape[0]
				INT = []
				if mask_image=='none':
				  for i, j, path in product(range(0,nx), range(0,ny), paths):
				    if path.contains_point((i,j)):
					  INT.append(inframe[j][i])
				else:
				  hdulist1 = pyfits.open(mask_image)
				  mask = hdulist1[0].data
				  for i, j, path in product(range(0,nx), range(0,ny), paths):
				    if (path.contains_point((i,j)) and mask[j][i]==0.):
					  INT.append(inframe[j][i])
				iso_reg.close()
				return sum(INT),np.mean(INT),np.median(INT),np.std(INT)
			      
			if 'ellipse(' in line:
				if mask_image!='none':
				    hdulist_infr = pyfits.open(mask_image)
				    mask_infr = hdulist_infr[0].data

				    mask_astropy = np.zeros_like(mask_infr, dtype=bool)
				    nx, ny =mask_infr.shape[1], mask_infr.shape[0]
				    for k in range(ny):
				      for i in range(nx):
					if mask_infr[k,i]>0:
					  mask_astropy[k,i] = True
					else:
					  mask_astropy[k,i] = False


				params = line[8:-2].split(",")
				xc = float(params[0])
				yc = float(params[1])
				ellA = float(params[2])
				ellB = float(params[3])
				try:
				  ellPA = float(params[4])
				except:
				  ellPA = float(params[4].split(')')[0])

				if ellA < ellB:
				    ellA, ellB = ellB, ellA
				    ellPA += 90 

			  	apertures = EllipticalAperture([(xc,yc)], ellA,ellB,radians(ellPA))
			  	if mask_image!='none':
				  phot_table = aperture_photometry(inframe, apertures,mask=mask_astropy)
				else:
				  phot_table = aperture_photometry(inframe, apertures)
				Flux_ell_masked = float(phot_table['aperture_sum'])
				iso_reg.close()
				return Flux_ell_masked,0.,0.,0.


def read_sextr(sextr_file):
  f = open(sextr_file, "r") 
  lines = f.readlines()
  def find_column():
    N_com =0
    for k in range(len(lines)):
      if '#' in lines[k]:
	N_com = N_com + 1
	if 'NUMBER' in lines[k]:
	  number_col = int(lines[k].split()[1])-1
	if 'X_WORLD' in lines[k]:
	  x_world_col = int(lines[k].split()[1])-1
	if 'Y_WORLD' in lines[k]:
	  y_world_col = int(lines[k].split()[1])-1
	if 'FLUX_MODEL' in lines[k]:
	  flux_model_col = int(lines[k].split()[1])-1

    return N_com,number_col,x_world_col,y_world_col,flux_model_col
  N_com,number_col,x_world_col,y_world_col,flux_model_col = find_column()
  #print N_com
  X = [];Y=[];flux_model=[];Number=[]
  for k in range(N_com,len(lines),1):
      star = lines[k].split()
      X.append(float(star[x_world_col]))
      Y.append(float(star[y_world_col]))
      flux_model.append(float(star[flux_model_col]))
      Number.append(int(star[number_col]))
  return X,Y,flux_model,Number


def replace(input_image,m0,GAIN,pix2sec,fwhm):
    #STEP 1: create region file
    print 'Please create region files with the first line for the galaxy and other lines for contaminants'
    ds9_show(input_image)
    
    reg_file = raw_input('Please enter the name of the region file [mask.reg]?') or 'mask.reg'
    
    #STEP 2: mask objects:
    output_image,output_mask = mask_indiv.mask(input_image,reg_file,skip_line=1)
    os.remove(output_image)

    hdulist2 = pyfits.open(output_mask)
    mask = hdulist2[0].data

    hdulist1 = pyfits.open(input_image)
    image_data = hdulist1[0].data
    ny,nx = image_data.shape

    nan_frame = np.copy(image_data)
    for k in range(ny):
      for i in range(nx):
	if mask[k,i]>0.:
	  nan_frame[k,i] = float('nan')
    filled = ds9_contour.replace_nans(nan_frame, 5, 0.5, kernel_size=2, method='localmean')
    
    hdu = pyfits.PrimaryHDU(filled)
    hdu.writeto('interp.fits',clobber=True)
    
    Flux,a,b,c = stats(filled,reg_file,mask_image='none',layer=0) 
    print 'Total flux: %.3f' % (Flux)
    
    check_se_res = 'no'
    while check_se_res!='yes':
      sextr.run_sextr(input_image,m0,GAIN,pix2sec,fwhm,se_file='galex_nuv.sex')
      ds9_show('models.fits',go_fur='yes')
      ds9_show('aper.fits',go_fur='yes')   
      check_se_res = raw_input('Are you happy with Sextractor results [YES/no]?') or 'yes'
    
    coords = ''
    X = []
    Y = []
    while coords!='q':
      coords = raw_input('Input coordinates (x,y) of objects to be subtracted from the flux of the galaxy. Press q for quit? ')
      if coords!='q':
	x,y = coords.split(',')
	X.append(int(x))
	Y.append(int(y))

    hdulist2 = pyfits.open('segm.fits')
    segm = hdulist2[0].data
    
      
    print "Reading Sextractor catalogue to match the objects"
    x_se,y_se,flux_model,Number = read_sextr('field.cat')
    Number = list(Number)
    
    F_sub = 0.
    for k in range(len(X)):
      number = segm[int(Y[k]),int(X[k])]
      print 'Flux for object %i: %.3f'  % (number,flux_model[Number.index(number)])
      F_sub = F_sub + flux_model[Number.index(number)]
    
    print 'Final flux of the galaxy (with the subtracted contaminants: %.3f DN\t%.3f mag' % (Flux-F_sub, m0-2.5*log10(Flux-F_sub))
    return Flux-F_sub,m0-2.5*log10(Flux-F_sub)
    
'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Isophote analysis")
    parser.add_argument("inputImage", help="Input fits image")
    parser.add_argument("ZeroPoint", nargs='?', const=1, help="Input Zero point in [mag arcsec^-2] (optional)",type=float,default=20.0)
    parser.add_argument("Gain", nargs='?', const=1, help="Input gain in [e-/ADU] (optional)",type=float,default=4.0)
    parser.add_argument("Scale", nargs='?', const=1, help="Input scale in [arcsec/pix] (optional)",type=float,default=1.0)
    parser.add_argument("FWHM", nargs='?', const=1, help="Input PSF FWHM (optional)",type=str,default=1.0)
    parser.add_argument("seFile", nargs='?', const=1, help="Input sextractor input file (optional)",type=str,default='default.sex') 

    args = parser.parse_args()

    input_image = args.inputImage
    m0 = args.ZeroPoint
    gain = float(args.Gain)
    pix2sec = float(args.Scale)
    FWHM = float(args.FWHM)
    seFile = args.seFile

    run_sextr(input_image,m0,gain,pix2sec,FWHM,se_file=seFile)   
'''
replace('nuv_crop.fits',20.08,5.,1.5,3.23)
#replace('fuv_crop.fits',18.82,5.,1.5,3.23)