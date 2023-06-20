#! /usr/bin/env python
import os
import sys

from pylab import *
import pyfits

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import glob
from matplotlib.pyplot import cm
import re
from scipy.odr.odrpack import *
import subprocess
import argparse

warnings.filterwarnings("ignore")

def sort_numbs(arr):
  numbers = []
  for k in range(len(arr)): 
    numb = str(arr[k].split('/')[2].split('_')[3])
    #print numb
    numbers.append(numb)
  a = sorted(numbers)
  new_arr = []
  for k in range(len(a)):
    ind = numbers.index(a[k])
    new_arr.append(arr[ind])
  return new_arr
    
    
    
def sort_bands(arr):
  new_arr = arr[:]
  out_arr = arr[:]
  bands = arr[:]
  bands_out = arr[:]
  #print arr
  for k in range(len(arr)): 
    band = str(arr[k].split('/')[2].split('_')[1])
    telescope = str(arr[k].split('/')[2].split('_')[0])
    #print telescope
    if band=='U':
      numb = 1
    elif band=='u':
      numb = 2
    elif band=='B':
      numb = 3    
    elif band=='g':
      numb = 4
    elif band=='V':
      numb = 5
    elif band=='R':
      numb = 6
    elif band=='r':
      numb = 7    
    elif band=='I':
      numb = 8
    elif band=='i':
      numb = 9
    elif band=='z':
      numb = 10
    elif band=='J':
      #print telescope,arr[k]
      numb = 11
      #print numb
    elif band=='H':
      numb = 12
    elif band=='K':
      if telescope=='2MASS':
	numb = 13
      else:
	numb = 14
    elif band=='w1' or band=='1':
      if telescope=='WISE':
	numb = 15
      else:
	numb = 16
    elif band=='w2' or band=='2':
      if telescope=='WISE':
	numb = 18
      else:
	numb = 17
    elif band=='w3':
      numb = 19
    elif band=='w4':
      numb = 20
    new_arr[k]=numb
    bands[k]=band
  a = sorted(new_arr)
  #print a
  for k in range(len(a)):
    number = new_arr.index(a[k])
    out_arr[k] = arr[number]
    bands_out[k] = out_arr[k].split('/')[2].split('_')[1]
    if (bands_out[k]=='w1' or bands_out[k]=='w2') and out_arr[k].split('/')[2].split('_')[0]=='Spitzer':
      if bands_out[k]=='w1':
	bands_out[k]='I1'
      if bands_out[k]=='w2':
	bands_out[k]='I2'
    if bands_out[k]=='1':
      bands_out[k]='I1'
  return out_arr,bands_out

def line_reg(header1,ima_pix2sec):
  nx = int(header1['NAXIS1'])
  ny = int(header1['NAXIS2'])
  scale = int(round(nx/8.*ima_pix2sec,-1))
  x2 = nx*9.8/10.
  x1 = x2 - scale/ima_pix2sec
  y1 = ny/7.
  y2 = y1
  
  return x1,y1,x2,y2,scale




def main(composed_model_file,min_level,pix2sec,m0):
  hdulist = pyfits.open(composed_model_file)
  primhdr = pyfits.getheader(composed_model_file, 1)
  inframe = hdulist[1].data
    
  hdu = pyfits.PrimaryHDU(-2.5*np.log10(inframe)+m0,primhdr)
  hdu.writeto('model_log.fits',clobber=True)

  primhdr = pyfits.getheader(composed_model_file, 0)
  inframe = hdulist[0].data
    
  noise_level = min_level
    
  max_level = -2.5*log10(np.max(inframe)/8.) + m0
    
  hdu = pyfits.PrimaryHDU(-2.5*np.log10(inframe)+m0,primhdr)
  hdu.writeto('ref_log.fits',clobber=True)
 
  primhdr = pyfits.getheader(composed_model_file, 3)
  inframe = hdulist[3].data
  ySize, xSize = inframe.shape
  header1 = hdulist[0].header
    
  hdu = pyfits.PrimaryHDU(inframe,primhdr)
  hdu.writeto('resid_log.fits',clobber=True)
  '''
  from astropy.io import fits
  data = np.zeros((1, xSize), dtype=np.float64)
  hdu = fits.PrimaryHDU(data=data)
  hdu.writeto('zero.fits',clobber=True)    
  '''
  
  zoom = 1
  #command1 = "ds9 -geometry %ix%i -zoom %.2f " % (1.01*xSize*zoom,260+3.*ySize*zoom,zoom)
  command1 = "ds9 -geometry %ix%i -tile row -zoom %.2f " % (1.01*xSize*zoom,1.5*(ySize*3)*zoom,zoom)
  command1 = command1 + "%s -scale limits %.3f %.3f -cmap A " % ('ref_log.fits',max_level,noise_level)
  x1,y1,x2,y2,scale = line_reg(header1,pix2sec)

  command1 = command1 + "-regions command \"line %.3f %.3f %.3f %.3f # color=black width=2\" " % (x1,y1,x2,y2)      
  box = "-regions command \"box %i %i %i %i 0 # color=blue width=6\" " % (xSize/2.,ySize/2.,xSize,ySize)
  command1 = command1 + box


  command1 = command1 + "%s -scale limits %.3f %.3f -cmap A " % ('model_log.fits',max_level,noise_level)
    
  box = "-regions command \"box %i %i %i %i 0 # color=blue width=6\" " % (xSize/2.,ySize/2.,xSize,ySize)
  command1 = command1 + box


  command1 = command1 + "%s -scale limits 0 1.0 -cmap i8 -invert " % ('resid_log.fits')
    
  box = "-regions command \"box %i %i %i %i 0 # color=blue width=6\" " % (xSize/2.,ySize/2.,xSize,ySize)
  command1 = command1 + box
  
  #command1 = command1 + "%s -scale limits 0 1.0 -cmap i8 -invert " % ('zero.fits')


    
  command1 = command1 + "-colorbar space distance -colorbar fontsize 7 -colorbar size 10 -saveimage png model.png -exit"
  subprocess.call(command1, shell=True)



  try:
      print "Converting to pdf..."
      subprocess.call('sam2p model.png  PDF: model.pdf', shell=True)
  except:
      print "ERROR: The program sam2p is not installed! To install sudo apt-get install sam2p"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create 2D image using the model and the reference image")
    parser.add_argument("model_file", help="Input model image")
    parser.add_argument("min_level", help="Input minimum isophote level to highlight the structure in [mag/arcsec^2]")
    parser.add_argument("Scale", help="Input scale in [arcsec/pix]")
    parser.add_argument("m0", help="Input Zero Point in [mag/arcsec^2]")     
    parser.add_argument("reference_file", nargs='?', const=0., help="Input reference image",type=str,default='none')  
    parser.add_argument("eps", nargs='?', const=0., help="Input yes if you want to convert image to eps",type=str,default='no')      

    args = parser.parse_args()
    
    composed_model_file = args.model_file
    min_level = float(args.min_level)
    pix2sec = float(args.Scale)
    m0 = float(args.m0)
    
    main(composed_model_file,min_level,pix2sec,m0)
