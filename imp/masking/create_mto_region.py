# Import the necessary modules
from astropy.io import fits as pyfits
#import pyfits
import numpy as np
import math
import itertools
from scipy import ndimage
import sys
from itertools import product
from matplotlib.path import Path
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import os
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from astropy import wcs

import read_mto_output



def main(mto_csv, output_region_file = 'mto_mask.reg'):
    ID,X,Y,A,B,theta,total_flux,mu_max,mu_median,mu_mean,R_fwhm,R_e,R10,R90,RA,DEC = read_mto_output.main(mto_csv)

    f = open(output_region_file, 'w')
    f.write('global color=green\n')
    f.write('image\n')
    for k in range(len(ID)):
        if not np.isnan(theta[k]):
            f.write('ellipse(%.4f,%.4f,%.4f,%.4f,%.2f) # text={%s}\n' % (X[k],Y[k],2.*A[k],2.*B[k],np.degrees(theta[k]),str(ID[k])))
    
    f.close()    
