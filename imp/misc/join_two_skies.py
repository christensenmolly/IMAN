#!/usr/bin/env python
# -*- coding: utf8 -*-
# python ~/MEGA/MyPrograms/IMAN/IMP_NEW/join_two_skies.py back_clean_H_in.fits back_clean_H_out.fits H_nan_out.reg H_final.fits
# Import the necessary modules
import pyfits
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
import imp_setup

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def rot_point(p, orig, angle):
    """ performs rotation of a point (x,y) around point (x0,y0)
        by angle radians"""
    x1 = cos(angle)*(p.x-orig.x) - sin(angle)*(p.y-orig.y) + orig.x
    y1 = sin(angle)*(p.x-orig.x) + cos(angle)*(p.y-orig.y) + orig.y
    return Point(x1, y1)

def main(clean_sky_inner, clean_sky_outer, region_file, output_name):
        HDU_inn = pyfits.open(clean_sky_inner)
        data_inn = HDU_inn[0].data
        header = HDU_inn[0].header
        ySize, xSize = data_inn.shape
    
        HDU_out = pyfits.open(clean_sky_outer)
        data_out = HDU_out[0].data
        
        data = np.copy(data_out)
   
        reg = open(region_file,'r')

		
        # Loop over all regions 
        for line in reg:
	  if "ellipse" in line and 'text={outer}' in line:
				      params = line.split(",")
				      cen = Point(float(params[0].split('(')[1]),
						  float(params[1]))
				      ellA = float(params[2])
				      ellB = float(params[3])
				      ellPA = float(params[4].split(')')[0])
				      if ellA < ellB:
					    ellA, ellB = ellB, ellA
					    ellPA += 90
				      
				      break


        cospa = cos(radians(ellPA))
        sinpa = sin(radians(ellPA))
        focusR = (ellA ** 2.0 - ellB ** 2.0) ** 0.5
        focus10 = Point(cen.x + focusR, cen.y)  # Unrotated
        focus20 = Point(cen.x - focusR, cen.y)  #
        focus1 = rot_point(focus10, cen, radians(ellPA))
        focus2 = rot_point(focus20, cen, radians(ellPA))
        # Find pixels inside of the ellipse
        dEll = 2 * ellA


  

        for x in xrange(0, xSize+1):
	    for y in xrange(0, ySize+1):
		    dFocus1 = hypot(x-focus1.x, y-focus1.y)
		    dFocus2 = hypot(x-focus2.x, y-focus2.y)
		    dPoint = dFocus1 + dFocus2
		    if dPoint <= dEll:
		      try:
			    data[y-1,x-1] = data_inn[y-1,x-1]

		      except:
			zz=1

        outHDU = pyfits.PrimaryHDU(data, header=header)
        outHDU.writeto(output_name, clobber=True)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])