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
from scipy.signal import savgol_filter

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

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def main(iso_reg, new_iso_reg):

        
        # 1. Read region file to download the isophotes
        reg = open(iso_reg,'r')


        x = []
        y = []
        for line in reg:
            #print(line)
            if 'circle(' in line:
                                params = line.split("(")[-1].split(",")[0]
                                xc =float(line.split("(")[-1].split(",")[0])
                                yc = float(line.split("(")[-1].split(",")[1])
                                ellA = float(line.split("(")[-1].split(",")[2].split(')')[0])
                                ellB = ellA
                                ellPA = 0.
                                x.append(xc)
                                y.append(yc)

        x = np.array(x)
        y = np.array(y)

        #y = savgol_filter(y, 15, 3)
        y = smooth(y, 5)
            
            

        if True:             
                points_ds9 = open(new_iso_reg,'w')
                for k in range(len(x)):
                    points_ds9.write("point(%.1f,%.1f) # point=box 8 color=blue width=4\n" % (x[k],y[k]))
                points_ds9.close()


main('combined_warp_points.reg', 'smoothed_skeleton_line.reg')