#!/usr/bin/python
# DESCRIPTION:
# Script to do automasking of galaxy image using Sextractor.
# Detected sources can be determined by changing snr and min_pix,
# as well as using a different input sextractor file.
# The target galaxy can be unmasked if needed.
# The masked areas can be presented as ellipses or polygons.
# The masked areas can be enlarged (by multiplication or subtraction).
# MINIMAL USAGE: python auto_masking.py [input_image]

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
#from photutils import data_properties, properties_table
import shapely
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from photutils.isophote import EllipseGeometry
from scipy.interpolate import LSQUnivariateSpline
from astropy.io import fits as pyfits

def make_ima(r,I,xc,yc,nx,ny, fill=float('nan')):
    # xc, yc should be given in pixel coordinates!
    # the target grid is spaced in 0.1 pixel intervals so as
    # to ensure no gaps will result on the output array.
    finely_spaced_sma = np.arange(r[0], r[-1], 0.1)

    # interpolate ellipse parameters

    # End points must be discarded, but how many?
    # This seems to work so far
    nodes = r[2:-2]

    intens_array = LSQUnivariateSpline(
        r, I, nodes)(finely_spaced_sma)



    result = np.zeros((ny,nx))
    weight = np.zeros((ny,nx))



    # for each interpolated isophote, generate intensity values on the
    # output image array
    # for index in range(len(finely_spaced_sma)):
    for index in range(1, len(finely_spaced_sma)):
        sma0 = finely_spaced_sma[index]
        print(sma0)
        x0 = xc-1 # Now in numpy
        y0 = yc-1 # Now in numpy
        geometry = EllipseGeometry(x0, y0, sma0, 0., 0.)

        intens = intens_array[index]

        # scan angles. Need to go a bit beyond full circle to ensure
        # full coverage.
        r = sma0
        phi = 0.
        while (phi <= 2*np.pi + geometry._phi_min):
            # get image coordinates of (r, phi) pixel
            x = r * np.cos(phi + 0.) + x0
            y = r * np.sin(phi + 0.) + y0
            i = int(x)
            j = int(y)

            # if outside image boundaries, ignore.
            if (i > 0 and i < nx - 1 and j > 0 and j < ny - 1):
                # get fractional deviations relative to target array
                fx = x - float(i)
                fy = y - float(j)

                # add up the isophote contribution to the overlapping pixels
                result[j, i] += (intens) * (1. - fy) * (1. - fx)
                result[j, i + 1] += (intens) * (1. - fy) * fx
                result[j + 1, i] += (intens) * fy * (1. - fx)
                result[j + 1, i + 1] += (intens) * fy * fx

                # add up the fractional area contribution to the
                # overlapping pixels
                weight[j, i] += (1. - fy) * (1. - fx)
                weight[j, i + 1] += (1. - fy) * fx
                weight[j + 1, i] += fy * (1. - fx)
                weight[j + 1, i + 1] += fy * fx

                # step towards next pixel on ellipse
                phi = max((phi + 0.75 / r), geometry._phi_min)
                r = geometry.radius(phi)
            else:
                break
    # zero weight values must be set to 1.
    weight[np.where(weight <= 0.)] = 1.

    # normalize
    result /= weight

    # fill value
    result[np.where(result == 0.)] = fill

    return result    


def make_ima1(r,I,xc,yc,nx,ny):
    result = np.zeros((ny,nx))

    for index in range(1, len(r)):
        sma0 = r[index]
        print(r[index])
        x0 = xc
        y0 = yc

        intens = I[index]

        # scan angles. Need to go a bit beyond full circle to ensure
        # full coverage.
        rr = sma0
        phi = 0.
        geometry = EllipseGeometry(x0, y0, sma0, 0., 0.)
        while (phi <= 2*np.pi + geometry._phi_min):
            
            
            # get image coordinates of (rr, phi) pixel
            x = rr * np.cos(phi) + x0
            y = rr * np.sin(phi) + y0
            i = int(x)
            j = int(y)

            # if outside image boundaries, ignore.
            if (i > 0 and i < nx - 1 and j > 0 and j < ny - 1):
                result[j, i] = intens

                # step towards next pixel on ellipse
                phi = max((phi + 0.75 / rr), geometry._phi_min)
            else:
                break
    print('here')
    return result    




'''
r = np.arange(0,500)
I = 1./(1.+r)
model = make_ima1(r,I,501,501,1001,1001)
print('Here')
pyfits.PrimaryHDU(data=model).writeto('test.fits', overwrite=True)
'''
