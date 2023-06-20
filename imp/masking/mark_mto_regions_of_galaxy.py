# Import the necessary modules
from astropy.io import fits as pyfits
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
import shapely
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import subprocess

import convert_reg_to_mask
from astropy import wcs

import read_mto_output

def main(mto_segm_image, xc,yc, mto_par_file):
    hdulist = pyfits.open(mto_segm_image)
    segm = hdulist[0].data
    
    ID,X,Y,A,B,theta,total_flux,mu_max,mu_median,mu_mean,R_fwhm,R_e,R10,R90,RA,DEC = read_mto_output.main(mto_par_file)
    circle_cen = Point(xc,yc).buffer(1)
    
    objects_to_mark = []
    
    objects_to_mark.append([xc,yc])
    for k in range(len(ID)):
        if not np.isnan(theta[k]):
                if (xc-X[k])**2 + (yc-Y[k])**2<=(2.*A[k])**2:
                    circle = Point(X[k],Y[k]).buffer(1)
                    ellipse = shapely.affinity.scale(circle,2.*A[k],2.*B[k])
                    rot_ellipse = shapely.affinity.rotate(ellipse, np.degrees(theta[k]), origin='center', use_radians=False)  

                    if rot_ellipse.contains(circle_cen): #rot_ellipse.contains(circle_cen):
                        inds = np.where(segm==ID[k])
                        objects_to_mark.append([inds[1][0],inds[0][0]])
                        #exit()
                        #objects_to_mark.append([X[k],Y[k]])
                        #print(ID[k],X[k],Y[k])

    return objects_to_mark


#objects_to_mark = main(443.,444., 'parameters.csv')

#print(objects_to_mark)
