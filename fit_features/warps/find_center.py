#!/usr/bin/python
# -*- coding:  cp1251 -*-

# Import standard modules
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
import re
import glob
import warnings
from scipy import special
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy import optimize

import pyfits
from scipy import interpolate
import argparse

import crea_iso

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
sys.path.append(PATH_TO_PACKAGE+'/FindFluxes')
import ds9_contour
import superell_fit

warnings.filterwarnings("ignore")

#*** Colour fonts ***
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

def main(input_image,m0,pix2sec,xc,yc,inner_level,outer_level): 
  print bcolors.OKBLUE+'\n\n************ Center of the galaxy on the isophotes map ************' + bcolors.ENDC
  isophote_file = crea_iso.main(input_image,m0,pix2sec,xc,yc,inner_level,outer_level)
  
  reg = open(isophote_file,'r')

  X_cen = []; Y_cen = []

  for Line in reg:
	  x = []
	  y = []
	  if 'polygon' in Line:
	    coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]
	    for k in range(0,len(coords)-1,2):
	      x.append(float(coords[k]))
	      y.append(float(coords[k+1]))
	  x=np.array(x)
	  y=np.array(y)

	  # FITTING OF THE ISOPHOTE WITH ELLIPSE
	  try:
		a = superell_fit.fitEllipse(x,y)
		center = ds9_contour.ellipse_center(a)
		axes = ds9_contour.ellipse_axis_length(a)
		phi = degrees(ds9_contour.ellipse_angle_of_rotation(a))
		xc = center[0]
		yc = center[1]
		sma = axes[0]
		smb = axes[1]
	  except:	    
		sma = float(nan)
		smb = float(nan)
		xc = float(nan)
		yc = float(nan)
	  X_cen.append(xc)
	  Y_cen.append(yc)

  print bcolors.OKGREEN+'Estimated center for given isophotes (X,std,Y,std): '+ bcolors.ENDC+ ' %.3f %.3f\t %.3f %.3f ' % (median(X_cen),std(X_cen),median(Y_cen),std(Y_cen))
  return median(X_cen),std(X_cen),median(Y_cen),std(Y_cen)

      
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Creating region file with isophotes")
    parser.add_argument("input_image", help="Input image")
    parser.add_argument("ZeroPoint", help="Input Zero Point in [mag/arcsec^2]") 
    parser.add_argument("Scale", help="Input scale in [arcsec/pix]") 
    parser.add_argument("x0", help="Input x-coordinate of the center")
    parser.add_argument("y0", help="Input y-coordinate of the center")    
    parser.add_argument("inner_level", help="Input inner isophote in [mag arcsec^-2]")
    parser.add_argument("outer_level", help="Input outer isophote in [mag arcsec^-2]")    
    
    
    args = parser.parse_args()
    
    input_image = args.input_image
    m0 = float(args.ZeroPoint)
    pix2sec = float(args.Scale)    
    x0 = float(args.x0)
    y0 = float(args.y0)
    inner_level = float(args.inner_level)
    outer_level = float(args.outer_level)
    
    main(input_image,m0,pix2sec,x0,y0,inner_level,outer_level)      
  
  