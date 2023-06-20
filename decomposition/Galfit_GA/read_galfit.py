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



def chi2_nu(number_of_pars,image_file,model_file,sigma_file,mask_file):
  '''
  subprocess.call("imarith " + image_file + " " + model_file + " sub out1.fits", shell=True)
  subprocess.call("imarith out1.fits out1.fits mul out2.fits", shell=True)
  subprocess.call("imarith " + sigma_file + " " + sigma_file + " mul out3.fits", shell=True)
  subprocess.call("imarith out2.fits out3.fits div out.fits", shell=True)
  '''
  
  hdulist1 = pyfits.open(image_file)
  image = hdulist1[0].data

  hdulist2 = pyfits.open(model_file)
  model = hdulist2[0].data

  hdulist3 = pyfits.open(sigma_file)
  sigma = hdulist3[0].data

  hdulist4 = pyfits.open(mask_file)
  mask = hdulist4[0].data
  
  ySize, xSize = image.shape
  SUM = 0.
  Npix = 1
  for k in range(ySize):
    for i in range(xSize):
      if mask[k,i]==0:
	SUM = SUM + (image[k,i]-model[k,i])**2 / (sigma[k,i]**2)
	Npix = Npix + 1
  
  chi2 = SUM/(Npix - number_of_pars)
  #print Npix - number_of_pars
  return chi2


'''

# Component number: 1
 0) edgedisk1 \            #  Component type
 1) 222.3553 47.9388  1 1  # 210 230 40 50 Position x, y
 3) 19.3620     1          # 20 21     Mu(0)   [mag/arcsec^2]
 4) 5.0291      1          #  h_s (disk scale-height)   [pix]
 5) 63.0240     1          #  R_s (disk scale-length)   [pix]
 6) 0.0000      0          #     ----- 
 7) 0.0000      0          #     ----- 
 8) 0.0000      0          #     ----- 
 9) 1.0000      -1         #     ----- 
10) 89.8860     1          #  Position angle (PA) [deg: Up=0, Left=90]
To)  2                     #  Outer truncation by component number(s)
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)
'''


'''

f = open('galfit.02', "r")
file_output = 'model.txt'
ff = open(file_output, "w")
'''
def main(mode,GalfitFile,InputGalfitFile):
  # GalfitFile - input Galfit file with the ranges of free parameters
  # InputGalfitFile - input Galfit file to create a new model (genotype)
  InputGalfitFile = 'model.txt'
  par = {}
  k = 0
  number_of_pars = 0
  with open(GalfitFile) as f:
      with open(InputGalfitFile, "w") as f1:  
	  for line in f:
		k = k + 1
		WRITE='no'
		if line!='\n':
		    if len(line.split())>4 and '#' not in line.split()[0] and 'H)' not in line.split()[0]:
		      #print line.split()[2],line.split()[3]
		      if line.split()[2]=='1' and line.split()[3]!='1':
			try:
			  PAR = float(line.split()[1])
			  PAR_lim_l = float(line.split()[4])
			  PAR_lim_r = float(line.split()[5])
			  par[str(k)] = [PAR,PAR_lim_l,PAR_lim_r]
			  LINE = ''
			  if mode=='ini':
			    PAR = random.uniform(PAR_lim_l,PAR_lim_r)
			    for i in range(len(line.split())):
			      if i==1:
				LINE = LINE + ' ' + str(PAR)
			      else:
				LINE = LINE + ' ' + line.split()[i]
			    LINE = LINE + '\n'
			  f1.write(LINE)
			  WRITE = 'yes'
			  number_of_pars = number_of_pars + 1
			except:
			  print 'There is an ERROR in line %i' % (k)
			  exit()      
		      elif line.split()[3]=='1' and line.split()[4]=='1':
			#Coordinates of the components
			try:
			  x = float(line.split()[1])
			  y = float(line.split()[2])
			  xmin = float(line.split()[6])
			  xmax = float(line.split()[7])
			  ymin = float(line.split()[8])
			  ymax = float(line.split()[9])
			  par[str(k)+'/1'] = [x,xmin,xmax]
			  par[str(k)+'/2'] = [y,ymin,ymax]
			  LINE = ''
			  if mode=='ini':
			    X = random.uniform(xmin,xmax)
			    Y = random.uniform(ymin,ymax)
			    for i in range(len(line.split())):
			      if i==1:
				LINE = LINE + ' ' + str(X)
			      elif i==2:
				LINE = LINE + ' ' + str(Y)
			      else:
				LINE = LINE + ' ' + line.split()[i]
			    LINE = LINE + '\n'
			  f1.write(LINE)
			  WRITE = 'yes'
			  number_of_pars = number_of_pars + 2
			except:
			  print 'There is an ERROR in line %i' % (k)
			  exit()
		      
		else:
		  f1.write(line)
		  WRITE = 'yes'

		if WRITE=='no':
		      #print 'heree'
		      f1.write(line)
  # sort bu the sky:
  '''
  PAr = {}
  for key in sorted(par):
    PAr[str(key)] = par[key]
    print "%s: %s" % (key,par[key]) 

  print PAr
  exit()
  '''
  return par,number_of_pars 


def read_input_files(GalfitFile):
  # This function returns names of input files:
  with open(GalfitFile) as f:
	  for line in f:
	    try:
		if line.split()[0]=='A)':
		  image_file = line.split()[1]
		if line.split()[0]=='B)':
		  model_file = line.split()[1]  
		if line.split()[0]=='C)':
		  sigma_file = line.split()[1]
		if line.split()[0]=='D)':
		  psf_file = line.split()[1]  
		if line.split()[0]=='F)':
		  mask_file = line.split()[1]
	    except:
	      z=1
  return image_file,model_file,sigma_file,psf_file,mask_file


#par,number_of_pars = main('ini',GalfitFile,InputGalfitFile)
#subprocess.call("galfit model.txt", shell=True)
#chi2 = chi2_nu(number_of_pars,'Spitzer_w1_a_r_ext_reb_new_galf.fits','model.fits','Spitzer_sigma_a_r_ext_reb.fits','Spitzer_w1_a_r_ext_mask_reb_new_mask.fits')





def input_new(mode,GalfitFile,InputGalfitFile,pars,free):
  if free=='free':
    NN = '1'
  else:
    NN = '0'

  #InputGalfitFile = 'model.txt'
  # GalfitFile - input Galfit file with the ranges of free parameters
  # InputGalfitFile - input Galfit file to create a new model (genotype)
  par = {}
  k = 0
  number_of_pars = 0
  with open(GalfitFile) as f:
      with open(InputGalfitFile, "w") as f1:  
	  for line in f:
		k = k + 1
		WRITE='no'
		if line!='\n':
		    if len(line.split())>4 and '#' not in line.split()[0] and 'H)' not in line.split()[0]:
		      #print line.split()[2],line.split()[3]
		      if line.split()[2]=='1' and line.split()[3]!='1':
			try:
			  PAR = float(line.split()[1])
			  PAR_lim_l = float(line.split()[4])
			  PAR_lim_r = float(line.split()[5])
			  par[str(k)] = [PAR,PAR_lim_l,PAR_lim_r]
			  LINE = ''

			  if mode=='ini':
			    #print pars
			    for key in pars:
			      #print key, str(k)
			      if key==str(k):
				PAR = pars[key]
			    for i in range(len(line.split())):
			      if i==1:
				LINE = LINE + ' ' + str(PAR)
			      elif i==2:
				LINE = LINE + ' ' + NN
			      else:
				LINE = LINE + ' ' + line.split()[i]
			    LINE = LINE + '\n'
			  f1.write(LINE)
			  WRITE = 'yes'
			  number_of_pars = number_of_pars + 1
			except:
			  print 'There is an ERROR in line %i' % (k)
			  exit()

		      elif line.split()[3]=='1' and line.split()[4]=='1':
			#Coordinates of the components
			try:
			  x = float(line.split()[1])
			  y = float(line.split()[2])
			  xmin = float(line.split()[6])
			  xmax = float(line.split()[7])
			  ymin = float(line.split()[8])
			  ymax = float(line.split()[9])
			  par[str(k)+'/1'] = [x,xmin,xmax]
			  par[str(k)+'/2'] = [y,ymin,ymax]
			  LINE = ''
			  if mode=='ini':
			    for key in pars:
			      if key==str(k)+'/1':
				X = pars[key]
			      if key==str(k)+'/2':
				Y = pars[key]   
			    #X = random.uniform(xmin,xmax)
			    #Y = random.uniform(ymin,ymax)
			    for i in range(len(line.split())):
			      if i==1:
				LINE = LINE + ' ' + str(X)
			      elif i==2:
				LINE = LINE + ' ' + str(Y)
			      elif i==3:
				LINE = LINE + ' ' + NN
			      elif i==4:
				LINE = LINE + ' ' + NN
			      else:
				LINE = LINE + ' ' + line.split()[i]
			    LINE = LINE + '\n'
			  f1.write(LINE)
			  WRITE = 'yes'
			  number_of_pars = number_of_pars + 2
			except:
			  print 'There is an ERROR in line %i' % (k)
			  exit()
		      
		else:
		  f1.write(line)
		  WRITE = 'yes'

		if WRITE=='no':
		      #print 'heree'
		      f1.write(line)
  # sort bu the sky:
  '''
  PAr = {}
  for key in sorted(par):
    PAr[str(key)] = par[key]
    print "%s: %s" % (key,par[key]) 

  print PAr
  exit()
  '''
  return par,number_of_pars 










































'''

def read_header(f,ff):
  lines = f.readlines()
  k = -1
  for line in lines:
    k = k + 1
    if line=='# IMAGE and GALFIT CONTROL PARAMETERS\n':
      for numb_line in range(k-1,k+32,1):
	print >>ff, lines[numb_line].split('\n')[0]
    break
  return k+31,lines
    #break








f = open('galfit.02', "r")
file_output = 'model.txt'
ff = open(file_output, "w")



line_numb,lines = read_header(f,ff)
print line_numb
par = {}
for k in range(line_numb,len(lines),1):
  line = lines[k].split('\n')
  LINE =  line[0]
  if LINE!='':
    if '1)' in LINE and LINE.split()[3]=='1' and LINE.split()[4]=='1':
      #print LINE.split()[1],LINE.split()[6]
      par[str(k)] = [LINE.split()[1],LINE.split()[6], LINE.split()[7],LINE.split()[2],LINE.split()[8],LINE.split()[9]]
    else:
     print LINE
     if LINE.split()[3]=='1':
       par[str(k)] = [LINE.split()[1],LINE.split()[4],LINE.split()[5]]
       print par
       exit()
print par
      
      
  



ff.close()
'''