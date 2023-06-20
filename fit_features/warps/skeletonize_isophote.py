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
import heapq
from astropy.io import fits as pyfits
from scipy import interpolate
import argparse

#PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
#PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
#sys.path.append(PATH_TO_PACKAGE+'/FindFluxes')
import polygon_reg
import crea_skeleton_new


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

def main(input_image, iso_reg):
        print(bcolors.OKBLUE+'\n\n************ Warp analysis 2015 ************' + bcolors.ENDC)
        print("Analyzing...")

        
        hdulist = pyfits.open(input_image)
        inframe = hdulist[0].data
        prihdr = hdulist[0].header
        nx, ny =inframe.shape[1], inframe.shape[0]
        
        # 1. Read region file to download the isophotes
        reg = open(iso_reg,'r')

        X_cen = []; Y_cen = []
        std_delta_y = 0.
        STD_DELTA_Y = []
        for Line in reg:
            x = []
            y = []
            if 'polygon' in Line:
                coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]
                for k in range(0,len(coords)-1,2):
                    x.append(float(coords[k]))
                    y.append(float(coords[k+1]))
            else:
                continue
            min_x_ind = x.index(min(x))
            max_x_ind = x.index(max(x))

            x = np.array(x)
            y = np.array(y)
            
            #Radius = max(np.sqrt((x-x0)**2+(y-y0)**2))
            if True:
                if len(STD_DELTA_Y)>2:
                    std_delta_y = mean(STD_DELTA_Y)
                print(Line)
                x_cen,y_cen,Std_delta_y,x_fit,y_fit = crea_skeleton_new.main(input_image,Line,std_delta_y)
                STD_DELTA_Y.append(Std_delta_y)
                print('Skeleton')
                # Picture
                points_reg = input_image.split('.fits')[0]+'_centerlineee.reg'
                
                points_ds9 = open(points_reg,'w')
                for k in range(len(x_fit)):
                    #if y_cen[k]>1500 and y_cen[k]<1540:
                    #print >>points_ds9, "point(%.1f,%.1f) # point=x color=red width=1" % (x_cen[k],y_cen[k])
                    points_ds9.write("circle(%.1f,%.1f,%.1f) # color=red width=1\n" % (x_cen[k],y_cen[k],1.0))
                points_ds9.close()
                exit()

#main('73_stacked.fits', 'stream.reg')
main('76_stacked.fits', 'stream.reg')