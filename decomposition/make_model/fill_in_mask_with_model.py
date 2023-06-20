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
tmp_out = sys.stdout
from scipy import stats
import re
from scipy.optimize import fsolve
import pyfits
import subprocess
import argparse



def main(input_model,mask,m0,pix2sec,D,output_file):
    hdu = pyfits.open(input_model)
    image_data = hdu[0].data
    primhdr = hdu[0].header
    model_data = hdu[1].data
    
    hdu_mask = pyfits.open(mask)
    mask_data = hdu_mask[0].data
    
    ny,nx = np.shape(image_data)
    
    new_image_data = np.copy(image_data)
    
    total_lum = 0.
    
    for k in range(ny):
      for i in range(nx):
	if mask_data[k,i]!=0:
	  new_image_data[k,i] = model_data[k,i]
	total_lum = total_lum + new_image_data[k,i]

    hdu = pyfits.PrimaryHDU(new_image_data,primhdr)
    if output_file==None:
      output_file = 'interp.fits'
    hdu.writeto(output_file,clobber=True)
    
    
    mag = m0 - 2.5*log10(total_lum)
    MAG = mag - 5.*log10(D) - 25.
    print "Apparent magnitude is:",mag
    print "Absolute magnitude is:",MAG

    return mag,MAG



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plotting profiles")
    parser.add_argument("input_model", help="Input composed model file")
    parser.add_argument("mask", help="Input mask file")
    
    parser.add_argument("ZeroPoint", help="Input Zero Point in [mag/arcsec^2]") 
    parser.add_argument("Scale", help="Input scale in [arcsec/pix]")
    parser.add_argument("D", help="Input distance in [Mpc]")
    parser.add_argument("--o", nargs='?', const=0., help="Name of the output file",type=str,default=None)
    
    args = parser.parse_args()

    input_model = args.input_model
    mask = args.mask
    m0 = float(args.ZeroPoint)
    pix2sec = float(args.Scale)
    D = float(args.D)
    output_file = args.o

    
    main(input_model,mask,m0,pix2sec,D,output_file)