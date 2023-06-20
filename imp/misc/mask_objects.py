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

# Main function, this part actually runs when routine is called
def ellipse_mask(cen, ellA, ellB, ellPA, inframe, xSize, ySize):
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
			inframe[y,x] = 0.0
		  except:
		    zz=1

def main(arr,region):
	value = 0.0
	norm = imp_setup.norm

	nframes = arr

	for i in range(len(nframes)):
		inputname = nframes[i]
		if norm:
		  outputname = './mask/'+nframes[i].split('/')[-1].split('.fits')[0]+"_norm.fits"
		else:
		  outputname = './mask/'+nframes[i].split('/')[-1].split('.fits')[0]+"_clean.fits"
		  
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
				params = line[8:-2].split(",")
				cen = Point(float(params[0]),
					    float(params[1]))
				ellA = float(params[2])
				ellB = float(params[3])
				ellPA = float(params[4])
				if ellA < ellB:
				    ellA, ellB = ellB, ellA
				    ellPA += 90
				ellipse_mask(cen, ellA, ellB, ellPA, inframe, xSize, ySize)

		
		# Loop over the image and all polygons and fill them with the mask value	
		nx, ny =inframe.shape[1], inframe.shape[0]
		for i, j, path in product(range(0,nx), range(0,ny), paths):
			inframe[j][i]=value if path.contains_point((i,j)) else inframe[j][i]
		
		# Normalise the frame and write it out
		if norm:
		    total = inframe.sum()
		    print "The total value of the frame is ", str(total)
		    if norm:
			    for i, j in product(range(0,nx), range(0,ny)):
				    inframe[j][i]=inframe[j][i]/total
		hdu = pyfits.PrimaryHDU(inframe,primhdr)
		hdu.writeto(outputname,clobber=True)
		
		if norm:
		    #ADD TOTAL FLUX TO THE END OF THE HEADER
		    hdulist = pyfits.open(outputname, do_not_scale_image_data=True, mode='update')
		    prihdr = hdulist[0].header
		    prihdr.append(('TOTAL_FLUX', total),end=True)
		    hdulist.flush()
		    
	mask_data = np.zeros_like(inframe)
	for i in range(ny):
		for k in range(nx):
			  if inframe[i,k] == value:
				mask_data[i,k] = 1
	hdu = pyfits.PrimaryHDU(mask_data,primhdr)
	hdu.writeto('./mask/mask.fits',clobber=True)
	
def main_mask_objects(inputname,region):
		value = 0.0
		  
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
				params = line[8:-2].split(",")
				cen = Point(float(params[0]),
					    float(params[1]))
				ellA = float(params[2])
				ellB = float(params[3])
				ellPA = float(params[4])
				if ellA < ellB:
				    ellA, ellB = ellB, ellA
				    ellPA += 90
				ellipse_mask(cen, ellA, ellB, ellPA, inframe, xSize, ySize)

		
		# Loop over the image and all polygons and fill them with the mask value	
		nx, ny =inframe.shape[1], inframe.shape[0]
		for i, j, path in product(range(0,nx), range(0,ny), paths):
			inframe[j][i]=value if path.contains_point((i,j)) else inframe[j][i]
		
		    
		mask_data = np.zeros_like(inframe)
		for i in range(ny):
			for k in range(nx):
				  if inframe[i,k] == value:
					mask_data[i,k] = 1
		#hdu = pyfits.PrimaryHDU(mask_data,primhdr)
		#hdu.writeto(outputname,clobber=True)
		return mask_data
	      

 