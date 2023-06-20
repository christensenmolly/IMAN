from astropy.io import fits as pyfits
import numpy as np
import math
import itertools
import matplotlib.pyplot as plt
from scipy import ndimage
import sys
from matplotlib.path import Path
from itertools import product
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import os
import subprocess
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from astropy import wcs
from astroquery import ned
import pickle
import collections
import glob



def check_that_within_iso(x,y, polygon):
    circle = Point(x,y).buffer(1)
    return polygon.contains(circle)


def convert_polygonline_to_shapleyregion(line):
    coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
    pol = []
    for kk in range(0,len(coords)-1,2):
        pol.append((int(float(coords[kk])),int(float(coords[kk+1]))))
    polygon = Polygon(pol)
    return polygon




def main(isophotes_file, input_image_1, colour_image,  mask_image_1, mask_image_2, SB_faint, SB_bright, m0_1, pix2sec, inner=False):
    f_reg = open(isophotes_file, "r")
    for line1 in f_reg:
        if 'polygon(' in line1:
                polygon = convert_polygonline_to_shapleyregion(line1)
    f_reg.close()

    hdulist_1 = pyfits.open(input_image_1)
    data_1 = hdulist_1[0].data
    header_1 = hdulist_1[0].header
    ny,nx = np.shape(data_1)
    
    hdulist_mask_1 = pyfits.open(mask_image_1)
    mask_1 = hdulist_mask_1[0].data    


    hdulist_2 = pyfits.open(colour_image)
    data_2 = hdulist_2[0].data
    header_2 = hdulist_2[0].header
    ny,nx = np.shape(data_2)
    
    hdulist_mask_2 = pyfits.open(mask_image_2)
    mask_2 = hdulist_mask_2[0].data    
    
    color = []
    if inner==False:
        for k in range(ny):
            for i in range(nx):
                if data_1[k,i] >= 10**(0.4*(m0_1-SB_faint)) * (pix2sec)**2 and data_1[k,i] <= 10**(0.4*(m0_1-SB_bright)) * (pix2sec)**2 and mask_1[k,i]==0 and mask_2[k,i]==0:
                    if check_that_within_iso(i, k, polygon):
                        color.append(  data_2[k,i]  )
    else:
        for k in range(ny):
            for i in range(nx):
                if data_1[k,i] > 10**(0.4*(m0_1-SB_bright)) * (pix2sec)**2 and mask_1[k,i]==0 and mask_2[k,i]==0:
                    if check_that_within_iso(i, k, polygon):
                        color.append(  data_2[k,i]  )        
    color = np.array(color)
    
    mean, median, stddev = sigma_clipped_stats(color)
    
    return mean, median, stddev
    
    