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

import astropy.io.fits as pyfits
from scipy import interpolate
import argparse

warnings.filterwarnings("ignore")

#smooth = 5
#nLevels = 30
#print 'HEREEEEEEEEEEEEE'
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

delta = 10.	# Find center of the polygon in this range (in pix) 

def main(input_image,m0,pix2sec,xc=None,yc=None,inner_level=20.,outer_level=26.,smooth=5,nLevels=30):
        print(bcolors.OKBLUE+'\n\n************ Isophotes ************' + bcolors.ENDC)
        print("Plotting...")
        try:
            inner_level = 10**(0.4*(m0-inner_level+5.*log10(pix2sec)))
            outer_level = 10**(0.4*(m0-outer_level+5.*log10(pix2sec)))
            print(inner_level, outer_level)

            subprocess.call("ds9 %s -scale log -cmap b -regions system image -contour yes -contour scale log -contour limits %.6f %.6f -contour smooth %i -contour nlevels %i -contour convert -regions save con.reg -exit" % (input_image,outer_level,inner_level,smooth,nLevels), shell=True)

            reg = open('con.reg','r')

            isophote_file = 'isophotes.reg'
            iso = open(isophote_file,'w')
            iso.write('image\n')
            for Line in reg:
                x = []
                y = []
                if 'polygon(' in Line:
                    coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]
                    for k in range(0,len(coords)-1,2):
                        x.append(float(coords[k]))
                        y.append(float(coords[k+1]))

                x = np.array(x)
                y = np.array(y)
                x0 = mean(x)
                y0 = mean(y)
                if xc==None and yc==None:
                    XC = x0
                    YC = y0
                else:
                    XC = xc
                    YC = yc             

                if fabs(x0-XC)<delta and fabs(y0-YC)<delta:
                    iso.write('polygon(')
                    for k in range(len(x)):
                        if k!=len(x)-1:
                            iso.write('%.3f,%.3f,' % (x[k],y[k]))
                        else:
                            iso.write('%.3f,%.3f' % (x[k],y[k]))
                    iso.write(')\n')
            iso.close()

            os.remove('con.reg')
            subprocess.call("ds9 %s -scale log -cmap b -regions isophotes.reg" % (input_image), shell=True)
            return isophote_file
        except:
            inner_level = 10**(0.4*(m0-inner_level+5.*log10(pix2sec)))
            outer_level = 10**(0.4*(m0-outer_level+5.*log10(pix2sec)))
            #print(inner_level, outer_level)

            subprocess.call("ds9 %s -scale log -cmap b -regions system image -contour yes -contour scale log -contour limits %.1f %.1f -contour smooth %i -contour nlevels %i -contour convert -regions save isophotes.reg" % (input_image,outer_level,inner_level,smooth,nLevels), shell=True)

            return 'isophotes.reg'

#main('combined.fits',28.29,0.396,xc=None,yc=None,inner_level=20.,outer_level=26.,smooth=11,nLevels=30)
#main('combined_someskysub_rot_crop.fits',22.5,0.262,xc=None,yc=None,inner_level=20.,outer_level=26.,smooth=11,nLevels=30) # NGC4452 Legacy
#main('combined_reconv_res_rot_crop.fits',28.098,0.396,xc=None,yc=None,inner_level=20.,outer_level=26.,smooth=11,nLevels=30) # NGC4452 SDSS
#main('combined_someskysub_rot_crop.fits',22.5,0.262,xc=None,yc=None,inner_level=20.,outer_level=26.,smooth=11,nLevels=30) # NGC4469 Legacy
#main('combined_someskysub_rot_crop.fits',22.5,0.262,xc=None,yc=None,inner_level=20.,outer_level=26.,smooth=11,nLevels=30) # NGC4469 Legacy
#main('galaxy_rot_crop.fits',28.66561,1.666,xc=459,yc=275,inner_level=20.,outer_level=27.,smooth=11,nLevels=40) # NGC3034 HERON


#main('combined_rot.fits',28.15,0.396,xc=1516,yc=1518,inner_level=25.,outer_level=25.0,smooth=15,nLevels=1) # NGC509 
#main('interpolated.fits',28.19,0.396,xc=701,yc=379,inner_level=25.,outer_level=25.0,smooth=15,nLevels=1) # NGC4452 
#main('interpolated.fits',28.29,0.396,xc=3031,yc=3033,inner_level=25.,outer_level=25.0,smooth=15,nLevels=1) # NGC4469 
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Creating region file with isophotes")
    parser.add_argument("input_image", help="Input image")
    parser.add_argument("ZeroPoint", help="Input Zero Point in [mag/arcsec^2]") 
    parser.add_argument("Scale", help="Input scale in [arcsec/pix]") 
    parser.add_argument("x0", help="Input x-coordinate of the center")
    parser.add_argument("y0", help="Input y-coordinate of the center")    
    parser.add_argument("inner_level", help="Input inner isophote in [mag arcsec^-2]")
    parser.add_argument("outer_level", help="Input outer isophote in [mag arcsec^-2]")    

    parser.add_argument("--smooth", help="Optional: Smoothness", type=int, default=15) 
    parser.add_argument("--nLevels", help="Optional: Number of isophotes between the given isophotes", type=int, default=1)        
    
    args = parser.parse_args()
    
    input_image = args.input_image
    m0 = float(args.ZeroPoint)
    pix2sec = float(args.Scale)    
    x0 = float(args.x0)
    y0 = float(args.y0)
    inner_level = float(args.inner_level)
    outer_level = float(args.outer_level)

    smooth = int(args.smooth)
    nLevels = int(args.nLevels)
    
    main(input_image, m0, pix2sec, xс, yс, inner_level, outer_level, smooth=smooth, nLevels=nLevels)
