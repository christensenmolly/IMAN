#!/usr/bin/env python
# -*- coding: utf8 -*-

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

# Main function, this part actually runs when routine is called
def ellipse_mask(cen, ellA, ellB, ellPA, inframe, xSize, ySize, out):
    if out=='no':
        cospa = cos(radians(ellPA))
        sinpa = sin(radians(ellPA))
        # Check if whole ellipse is inside of the image
        # and obtain size of the ellipse in xy plane
        xMax = 0.0
        yMax = 0.0
        xMin = 1e10
        yMin = 1e10
        for e in linspace(0, 4*pi, 1000):
            cose = cos(e)
            sine = sin(e)
            x = cen.x + ellA * cose * cospa - ellB * sine * sinpa
            y = cen.y + ellB * sine * cospa + ellA * cose * sinpa
            if x > xMax:
                xMax = x
            if y > yMax:
                yMax = y
            if x < xMin:
                xMin = x
            if y < yMin:
                yMin = y
        xMin = max(0, int(round(xMin)))
        xMax = min(xSize, int(round(xMax)))
        yMin = max(0, int(round(yMin)))
        yMax = min(ySize, int(round(yMax)))
        focusR = (ellA ** 2.0 - ellB ** 2.0) ** 0.5
        focus10 = Point(cen.x + focusR, cen.y)  # Unrotated
        focus20 = Point(cen.x - focusR, cen.y)  #
        focus1 = rot_point(focus10, cen, radians(ellPA))
        focus2 = rot_point(focus20, cen, radians(ellPA))
        # Find pixels inside of the ellipse
        dEll = 2 * ellA
        for x in xrange(xMin, xMax+1):
            for y in xrange(yMin, yMax+1):
                dFocus1 = hypot(x-focus1.x, y-focus1.y)
                dFocus2 = hypot(x-focus2.x, y-focus2.y)
                dPoint = dFocus1 + dFocus2
                if dPoint < dEll:
		  try:
			inframe[y,x] = 1.0
		  except:
		    zz=1
    else:
        cospa = cos(radians(ellPA))
        sinpa = sin(radians(ellPA))
        # Check if whole ellipse is inside of the image
        # and obtain size of the ellipse in xy plane
        xMax = 0.0
        yMax = 0.0
        xMin = 1e10
        yMin = 1e10
        for e in linspace(0, 4*pi, 1000):
            cose = cos(e)
            sine = sin(e)
            x = cen.x + ellA * cose * cospa - ellB * sine * sinpa
            y = cen.y + ellB * sine * cospa + ellA * cose * sinpa
            if x > xMax:
                xMax = x
            if y > yMax:
                yMax = y
            if x < xMin:
                xMin = x
            if y < yMin:
                yMin = y
        xMin = max(0, int(round(xMin)))
        xMax = min(xSize, int(round(xMax)))
        yMin = max(0, int(round(yMin)))
        yMax = min(ySize, int(round(yMax)))
        focusR = (ellA ** 2.0 - ellB ** 2.0) ** 0.5
        focus10 = Point(cen.x + focusR, cen.y)  # Unrotated
        focus20 = Point(cen.x - focusR, cen.y)  #
        focus1 = rot_point(focus10, cen, radians(ellPA))
        focus2 = rot_point(focus20, cen, radians(ellPA))
        # Find pixels inside of the ellipse
        dEll = 2 * ellA
        ny,nx = np.shape(inframe)
        for x in xrange(0, nx):
            for y in xrange(0, ny):
                dFocus1 = hypot(x-focus1.x, y-focus1.y)
                dFocus2 = hypot(x-focus2.x, y-focus2.y)
                dPoint = dFocus1 + dFocus2
                if dPoint >= dEll:
		  try:
			inframe[y,x] = 1.0
		  except:
		    zz=1
    #return inframe





def main_mask_objects(inputname,region):
		value = 1.0
		  
		hdulist = pyfits.open(inputname)
		primhdr = pyfits.getheader(inputname, 0)
		inframe = hdulist[0].data
		reg = open(region,'r')
		paths=[]
		ySize, xSize = inframe.shape
		
		# Loop over all regions 
		for line in reg:
		
			# Detect circle regions and fill them with the mask value
			if 'circle(' in line:
				param = ((line.split('circle(')[1]).split(')')[0]).split(',')
				a, b ,r  = int(float(param[0])), int(float(param[1])), int(float(param[2])) 
				y,x = np.ogrid[-b:inframe.shape[0]-b, -a:inframe.shape[1]-a]
				if 'OUTER' in line:
				  mask = x*x + y*y >= r*r
				else:
				  mask = x*x + y*y <= r*r
				inframe[mask] = value
				
			# Detect polygon regions and add them to the path list
			elif 'polygon(' in line:
				param = map(float, ((line.split('polygon(')[1]).split(')')[0]).split(','))
				param2 = [None]*(len(param)/2)
				for i in xrange(0,len(param2)):	param2[i] = (int(param[2*i]),int(param[2*i+1])) 
				param2.append(param2[0])
				codes = []
				codes.append(Path.MOVETO)
				for i in xrange(1,len(param2)-1): codes.append(Path.LINETO)
				codes.append(Path.CLOSEPOLY)
				path = Path(param2, codes)
				paths.append(path)
				
			elif "ellipse" in line:
				params = line.split(",")
				cen = Point(float(params[0].split('(')[1]),
					    float(params[1]))
				ellA = float(params[2])
				ellB = float(params[3])
				ellPA = float(params[4].split(')')[0])
				if ellA < ellB:
				    ellA, ellB = ellB, ellA
				    ellPA += 90
				if 'OUTER' in line:
				  out = 'yes'
				else:
				  out = 'no'
				ellipse_mask(cen, ellA, ellB, ellPA, inframe, xSize, ySize, out)

		
		# Loop over the image and all polygons and fill them with the mask value	
		nx, ny =inframe.shape[1], inframe.shape[0]
		for i, j, path in product(range(0,nx), range(0,ny), paths):
			inframe[j][i]=value if path.contains_point((i,j)) else inframe[j][i]

		return inframe
	      


def add_mask_to_mask(input_mask,input_reg_file):
  new_mask_name = input_reg_file.split('.reg')[0]+'.fits'
  
  new_mask_data = main_mask_objects(input_mask,input_reg_file)

  maskedHDU = pyfits.PrimaryHDU(data=new_mask_data)
  maskedHDU.writeto(new_mask_name,clobber=True)
  return new_mask_name


add_mask_to_mask('mask.fits','add_mask.reg')