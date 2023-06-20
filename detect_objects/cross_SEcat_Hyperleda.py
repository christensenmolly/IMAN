#!/usr/bin/env python
# DESCRIPTION:
# Script to cross-correlate SE catalog and HyperLeda.


# Import the necessary modules
import pyfits
import numpy as np
import math
import sys
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import os
import subprocess
from astropy import coordinates
import astropy.units as u
from astropy import wcs
from photutils import find_peaks
from astropy.stats import sigma_clip
import warnings

LOCAL_DIR = "/detect_objects"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'ap_photometry'))
sys.path.append(os.path.join(IMAN_DIR, 'sextractor'))

import imp_rebin
import find_good_psf_stars
import read_SEcat


warnings.filterwarnings("ignore")

database= '/home/amosenko/CurrentWork/HYPERLEDA_DATABASE/G_logd25_larger_0.dat' #### WARNING: CHANGE THIS



def read_hyperleda(database, col_names, dtypes=str):
  print 'Reading the database ...'
  f = open(database, 'r')
  #lines = f.readlines()
  #print lines[419130:419133]
  #exit()
  f.readline()
  names = f.readline().split('|')
  names = names[0:len(names)-2]
  Names = []
  for k in range(len(names)):
    Names.append(names[k].strip())

  column_numbers = []
  for col_name in col_names:
    for k in range(len(Names)):
      if col_name == Names[k]:
	column_numbers.append(k)

  
  #print column_numbers
  #exit()
  if isinstance(dtypes, list):
    COLUMNS = []
    columns = np.loadtxt(database, usecols=column_numbers, unpack=True, dtype=str, skiprows=3, delimiter='|',comments='<')
    
    for k in range(len(columns)):
      COLUMNS.append(columns[k])
      COLUMNS[k] = np.array(COLUMNS[k], dtype=dtypes[k])
    
  else:
    columns = np.loadtxt(database, usecols=column_numbers, unpack=True, dtype=dtypes, skiprows=3, delimiter='|',comments='<')
    COLUMNS = columns
  f.close()

  return COLUMNS  



def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx],idx

def read_ellipse(region_file):
   xc = []; yc = []; sma = []; smb = []; PA = []; SEXTR_NUMB = []
   f= open(region_file,'r')
   for line in f:
        if "ellipse" in line:
	  params = line.split(",")
	  xc.append(float(params[0].split('(')[1]))
	  yc.append(float(params[1]))
	  sma.append(float(params[2]))
	  smb.append(float(params[3]))
	  PA.append(float(params[4].split(')')[0]))
	  
	  if 'text={' in line:
		  SEXTR_NUMB.append(str(line.split('text={')[-1].split(',')[0]))
	  else:
		  SEXTR_NUMB.append('nan')
   f.close()
   return xc,yc,sma,smb,PA,SEXTR_NUMB


def find_hyper_leda_objects(input_image, pix2sec, sma, pgc,al2000,de2000,logd25):
    print 'Find objects in the field using HyperLeda catalogue'
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    header = hdulist[0].header
    ny,nx = np.shape(data)
    xc = nx/2.; yc = ny/2.

    
    if 'COMMENT' in header:
	del header['COMMENT']
    if 'HISTORY' in header:
	del header['HISTORY']
    if '' in header:
	del header['']  
    w = wcs.WCS(header)
    pixcrd = np.array([[xc,yc]], np.float_)
    world = w.wcs_pix2world(pixcrd, 1)
    
    # Select object within the cerain Radius
    Radius = math.sqrt(xc**2+yc**2) * pix2sec
    
    HYPERLEDA_OBJECTS = []
    for k in range(len(pgc)):      
      r = math.sqrt( (world[0][0]-al2000[k])**2 + (world[0][1]-de2000[k])**2 ) * 3600.	# arcsec
      if r<Radius and 10**logd25[k]*6.> sma:
	HYPERLEDA_OBJECTS.append(k)
    hdulist.close()
    return HYPERLEDA_OBJECTS
    




def main(input_image, sex_catalog, sma=60., user_inter=False, add_hyper=False, loaded_database=None, border_dist=None):
  hdulist = pyfits.open(input_image)
  data = hdulist[0].data
  header = hdulist[0].header
  ny,nx = np.shape(data)
  xc = nx/2.; yc = ny/2.
  if border_dist!=None:
    dx = nx*border_dist
    dy = ny*border_dist
  else:
    dx = 0.
    dy = 0.

  pix2sec,note = imp_rebin.resolution(input_image)
  pix2sec = float(pix2sec)
    
  if 'COMMENT' in header:
	del header['COMMENT']
  if 'HISTORY' in header:
	del header['HISTORY']
  if '' in header:
	del header['']  
  w = wcs.WCS(header)
  
  hdulist.close()
  
  # Read Sextractor catalogue  
  X_IMAGE,Y_IMAGE,X_WORLD,Y_WORLD,A_IMAGE,B_IMAGE,KRON_RADIUS,THETA_IMAGE,CLASS_STAR = read_SEcat.find_sex_column(sex_catalog, ['X_IMAGE','Y_IMAGE','X_WORLD','Y_WORLD','A_IMAGE','B_IMAGE','KRON_RADIUS','THETA_IMAGE','CLASS_STAR'], float)
  print len(X_IMAGE)
  
  if loaded_database==None:
    # Download HyperLeda database:
    pgc,al2000,de2000,logd25,logr25,pa = read_hyperleda(database, ['pgc','al2000','de2000','logd25','logr25','pa'], dtypes=[int,float,float,float,float,float])
    al2000 = al2000*15.
  else:
    pgc,al2000,de2000,logd25,logr25,pa = loaded_database
  
  if add_hyper==True:
    HYPERLEDA_OBJECTS = find_hyper_leda_objects(input_image, pix2sec, sma, pgc,al2000,de2000,logd25)
    count = 0
    for i in range(len(HYPERLEDA_OBJECTS)):
      KK = HYPERLEDA_OBJECTS[i]
      r = np.sqrt( (X_WORLD-al2000[KK])**2 + (Y_WORLD-de2000[KK])**2 ) * 3600.	# arcsec
      if np.min(r)>10.:
	world = np.array([[al2000[KK],de2000[KK]]], np.float_)
	pixcrd = w.wcs_world2pix(world, 1)
	#print world,pixcrd
	if pixcrd[0][0]>0 and pixcrd[0][0]<nx and pixcrd[0][1]>0 and pixcrd[0][1]<ny:
	  np.append(X_IMAGE,[pixcrd[0][0]])
	  np.append(Y_IMAGE,[pixcrd[0][1]])
	  np.append(X_WORLD,[al2000[KK]])
	  np.append(Y_WORLD,[de2000[KK]])
	  if np.isnan(logd25[KK])==False:
	    np.append(A_IMAGE,[10**logd25[KK]*6./pix2sec])
	  else:
	    np.append(A_IMAGE,[10.])
	  if np.isnan(logd25[KK])==False and logr25[KK]==False:
	    np.append(B_IMAGE,[10**logd25[KK]*6./(pix2sec*10**(logr25[KK]))])
	  else:
	    np.append(B_IMAGE,[10.])
	  np.append(KRON_RADIUS,[1.])
	  if np.isnan(pa[KK])==False:
	    np.append(THETA_IMAGE,[pa[KK]-90.])
	  else:
	    np.append(THETA_IMAGE,[0.])
	  np.append(CLASS_STAR,[-1.])
	  count = count + 1
    print 'Added %i object from HyperLeda' % (count)
    

  fout = open("galaxies.reg", "w")
  fout.truncate(0)
  fout.write("# Region file format: DS9 version 4.1\n")
  fout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 ')
  fout.write('normal roman" select=1 highlite=1 dash=0 fixed=0 ')
  fout.write('edit=1 move=1 delete=1 include=1 source=1\n')
  fout.write('image\n')
  
  count = 0
  X = []; Y = []; RA = []; DEC = []; SMA = []; SMB = []; PA = []; NAME = []; SE_NUMB = []

  for k in range(len(X_IMAGE)):
    ellA = A_IMAGE[k] * KRON_RADIUS[k] 
    ellB = B_IMAGE[k] * KRON_RADIUS[k]
    ellPA = THETA_IMAGE[k]
    if ellA > sma and CLASS_STAR[k]<0.1:
      r = np.sqrt( (X_WORLD[k]-al2000)**2 + (Y_WORLD[k]-de2000)**2 ) * 3600.	# arcsec
      R_near,idx_near = find_nearest(r, 0.)
      borders = find_good_psf_stars.ellipse_borders([X_IMAGE[k], Y_IMAGE[k]],ellA,ellB,ellPA)
      
      if (R_near<10. or CLASS_STAR[k]==-1.) and (borders[0][0]>dx and borders[0][1]>dy and borders[1][0]<nx-dx and borders[1][1]<ny-dy):	# arcsec
	count = count + 1
	fout.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # text={%i,pgc%s}\n" % (X_IMAGE[k], Y_IMAGE[k],
									    ellA, ellB,
									    ellPA, k, pgc[idx_near]))
	print 'Galaxy #%i: pgc%s %i %i %.5f %.5f' % (k+1,pgc[idx_near],X_IMAGE[k],Y_IMAGE[k],X_WORLD[k],Y_WORLD[k])
	
	if user_inter==False:
	    X.append(X_IMAGE[k]); Y.append(Y_IMAGE[k])
	    RA.append(X_WORLD[k]); DEC.append(Y_WORLD[k])
	    SMA.append(ellA); SMB.append(ellB)
	    PA.append(ellPA)
	    NAME.append('pgc'+str(pgc[idx_near]))
	    SE_NUMB.append(str(k))
  fout.close()
  
  print 'TOTAL NUMBER OF GALAXIES IS: %i' % (count)
  
  
  if user_inter==True:
    # Read the image
    hdulist = pyfits.open(input_image)#, mode='update')
    data = hdulist[0].data
    header = hdulist[0].header
    ny,nx = np.shape(data)    
  
    print '\nPlease check the region obejcts (remove or add new ones as ellipses)'
    # View in DS9
    ds9Proc = subprocess.Popen(["ds9", input_image,
				  "-regions", "galaxies.reg",
				  "-scale", "histequ",
				  "-cmap","Rainbow"])
    ds9Proc.wait()        
    
    xc,yc,sma,smb,PAA,SEXTR_NUMB = read_ellipse("galaxies.reg")
    pixcrd = []
    for i in range(len(xc)):
      pixcrd.append([xc[i],yc[i]])

    if 'COMMENT' in header:
	del header['COMMENT']
    if 'HISTORY' in header:
	del header['HISTORY']
    if '' in header:
	del header['']  
    w = wcs.WCS(header)
    pixcrd = np.array(pixcrd, np.float_)
    
    # Convert pixel coordinates :to wcs
    world = w.wcs_pix2world(pixcrd, 1)

    fout = open("galaxies.reg", "w")
    fout.truncate(0)
    fout.write("# Region file format: DS9 version 4.1\n")
    fout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 ')
    fout.write('normal roman" select=1 highlite=1 dash=0 fixed=0 ')
    fout.write('edit=1 move=1 delete=1 include=1 source=1\n')
    fout.write('image\n')

    count = 0
    X = []; Y = []; RA = []; DEC = []; SMA = []; SMB = []; PA = []; NAME = []; SE_NUMB = []
    for k in range(len(xc)):
      ra = world[k][0]
      dec = world[k][1]
      
      r = np.sqrt( (ra-al2000)**2 + (dec-de2000)**2 ) * 3600.	# arcsec
      R_near,idx_near = find_nearest(r, 0.)
      count = count + 1

      if R_near<10.:	# arcsec
	  PGC_NAME = 'pgc'+str(pgc[idx_near])
      else:
	  PGC_NAME = 'nan'
	  
      fout.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # text={%s,%s}\n" % (xc[k], yc[k],
									    sma[k], smb[k],
									    PAA[k],SEXTR_NUMB[k],PGC_NAME))

      if SEXTR_NUMB[k]=='nan':
	print 'Galaxy #%s: %s %i %i %.5f %.5f' % (SEXTR_NUMB[k],PGC_NAME,xc[k],yc[k],ra,dec)
      else:
	print 'Galaxy #%i: %s %i %i %.5f %.5f' % (int(SEXTR_NUMB[k])+1,PGC_NAME,xc[k],yc[k],ra,dec)
	

      X.append(xc[k]); Y.append(yc[k])
      RA.append(ra); DEC.append(dec)
      SMA.append(sma[k]); SMB.append(smb[k])
      PA.append(PAA[k]);NAME.append(PGC_NAME)
      SE_NUMB.append(SEXTR_NUMB[k])
    fout.close()
    hdulist.close()
    print 'TOTAL NUMBER OF GALAXIES IS: %i' % (count)    


  
  
  
  return X,Y,RA,DEC,SMA,SMB,PA,NAME,SE_NUMB
  

#main('new-image.fits', 'field.cat', user_inter=False, add_hyper=True)  
