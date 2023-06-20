#!/usr/bin/python

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
import shutil
import argparse
import os
import subprocess
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats
from photutils import data_properties, properties_table
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon




class PPoint:
    def __init__(self, x, y):
        self.x = x
        self.y = y










def rot_point(p, orig, angle):
    """ performs rotation of a point (x,y) around point (x0,y0)
        by angle radians"""
    x1 = cos(angle)*(p.x-orig.x) - sin(angle)*(p.y-orig.y) + orig.x
    y1 = sin(angle)*(p.x-orig.x) + cos(angle)*(p.y-orig.y) + orig.y
    return PPoint(x1, y1)

def listRightIndex(alist, value):
    return len(alist) - alist[-1::-1].index(value) -1


def ellipse_mask1(cen, ellA, ellB, ellPA, inframe, xSize, ySize, img_inp=None,img_mask=None, outer = False):
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
        focus10 = PPoint(cen.x + focusR, cen.y)  # Unrotated
        focus20 = PPoint(cen.x - focusR, cen.y)  #
        focus1 = rot_point(focus10, cen, radians(ellPA))
        focus2 = rot_point(focus20, cen, radians(ellPA))
        # Find pixels inside of the ellipse
        dEll = 2 * ellA
        
        if outer==True:
	  xMin = 0
	  xMax = xSize -1
	  yMin = 0
	  yMax = ySize -1
	  
	  if img_inp==None:
	    for x in xrange(xMin, xMax+1):
		for y in xrange(yMin, yMax+1):
		    dFocus1 = hypot(x-focus1.x, y-focus1.y)
		    dFocus2 = hypot(x-focus2.x, y-focus2.y)
		    dPoint = dFocus1 + dFocus2
		    if dPoint >= dEll:
		      try:
			    inframe[y-1,x-1] = 0.0
			    img_mask[y-1,x-1] = 1
		      except:
			zz=1
		    else:
		      try:
			    inframe[y-1,x-1] = inframe[y-1,x-1]
			    img_mask[y-1,x-1] = 0
		      except:
			zz=1      
	  else:
	    for x in xrange(xMin, xMax+1):
		for y in xrange(yMin, yMax+1):
		    dFocus1 = hypot(x-focus1.x, y-focus1.y)
		    dFocus2 = hypot(x-focus2.x, y-focus2.y)
		    dPoint = dFocus1 + dFocus2
		    if dPoint >= dEll:
		      try:
			    inframe[y-1,x-1] = img_inp[y-1,x-1]
		      except:
			zz=1

	else:
	  if img_inp==None:
	    for x in xrange(xMin, xMax+1):
		for y in xrange(yMin, yMax+1):
		    dFocus1 = hypot(x-focus1.x, y-focus1.y)
		    dFocus2 = hypot(x-focus2.x, y-focus2.y)
		    dPoint = dFocus1 + dFocus2
		    if dPoint < dEll:
		      try:
			    inframe[y-1,x-1] = 0.0
			    img_mask[y-1,x-1] = 1
		      except:
			zz=1
	  else:
	    for x in xrange(xMin, xMax+1):
		for y in xrange(yMin, yMax+1):
		    dFocus1 = hypot(x-focus1.x, y-focus1.y)
		    dFocus2 = hypot(x-focus2.x, y-focus2.y)
		    dPoint = dFocus1 + dFocus2
		    if dPoint < dEll:
		      try:
			    inframe[y-1,x-1] = img_inp[y-1,x-1]
		      except:
			zz=1
        return inframe

def ellipse_annulus_mask(inner_ellipse, outer_ellipse, mask):
        ySize,xSize = np.shape(mask)
        MASK = np.ones(shape=(ySize,xSize))
	def define_ellipse(cen, ellA, ellB, ellPA):
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
	    return xMin,xMax,yMin,yMax,focus1,focus2

	cen_in, ellA_in, ellB_in, ellPA_in = inner_ellipse
	cen_out, ellA_out, ellB_out, ellPA_out = outer_ellipse
	
	## Inner ellipse
	xMin_in,xMax_in,yMin_in,yMax_in,focus1_in,focus2_in = define_ellipse(cen_in, ellA_in, ellB_in, ellPA_in)

	## Outer ellipse
	xMin_out,xMax_out,yMin_out,yMax_out,focus1_out,focus2_out = define_ellipse(cen_out, ellA_out, ellB_out, ellPA_out)
	
	xMin = min([xMin_in,xMin_out])
	yMin = min([yMin_in,yMin_out])

	xMax = max([xMax_in,xMax_out])
	yMax = max([yMax_in,yMax_out])

        # Find pixels inside of the ellipse annulus
        dEll_in = 2 * ellA_in
        dEll_out = 2 * ellA_out
        

        for x in xrange(xMin, xMax):
            for y in xrange(yMin, yMax):
                dFocus1_in = hypot(x-focus1_in.x, y-focus1_in.y)
                dFocus2_in = hypot(x-focus2_in.x, y-focus2_in.y)
                dPoint_in = dFocus1_in + dFocus2_in

                dFocus1_out = hypot(x-focus1_out.x, y-focus1_out.y)
                dFocus2_out = hypot(x-focus2_out.x, y-focus2_out.y)
                dPoint_out = dFocus1_out + dFocus2_out                
                if dPoint_out < dEll_out and dPoint_in > dEll_in and mask[y,x]==0.:
		  try:
			MASK[y,x] = 0.
		  except:
		    zz=1
	return MASK



def mask(input_image, reg_file, factor=1, output_image=None, output_mask=None, obj_ellipse=None, mask_value=1, show_running=True, mask_DN=None):
		if show_running==True:
		  print 'Masking contaminants...'
		if output_image==None:
		  output_image = input_image.split('.fits')[0]+'_clean.fits'  

		shutil.copy(input_image,output_image)
		
		if output_mask!=None:
		  shutil.copy(input_image,output_mask)
		  hdulist2 = pyfits.open(output_mask, do_not_scale_image_data=True, mode='update')
		  img2 = hdulist2[0].data
  
		hdulist3 = pyfits.open(output_image, do_not_scale_image_data=True, mode='update')
		img3 = hdulist3[0].data


		(dimy,dimx) = img3.shape
		value = 0#float('nan')
		paths=[]
		borders=[]
		f = open(reg_file, "r")
		all_lines = f.readlines()
		Number_of_lines = len(all_lines)
		f.close()
		f = open(reg_file, "r")

		coun = 0
		if factor==1:
		    for line in f:
			    #print line
			    # Detect circle regions and fill them with the mask value
			    if 'circle(' in line:
				      '''
				      param = ((line.split('circle(')[1]).split(')')[0]).split(',')
				      a, b ,r  = int(float(param[0])), int(float(param[1])), int(float(param[2])) 
				      y,x = np.ogrid[-b:img3.shape[0]-b, -a:img3.shape[1]-a]
				      mask = x*x + y*y <= r*r
				      img3[mask] = value
				      '''
				      params = line.split("(")[-1].split(",")[0]
				      cen = PPoint(float(line.split("(")[-1].split(",")[0]),
						  float(line.split("(")[-1].split(",")[1]))
				      ellA = float(line.split("(")[-1].split(",")[2].split(')')[0])
				      ellB = ellA
				      ellPA = 0.
				      ellipse_mask(cen, ellA, ellB, ellPA, img3, dimx, dimy)
    
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
				      coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
				      X = []; Y = []
				      for k in range(0,len(coords)-1,2):
					      X.append(int(float(coords[k])))
					      Y.append(int(float(coords[k+1])))

				      borders.append([min(X),min(Y),max(X),max(Y)])
				    
			    elif "ellipse" in line:
				      params = line.split(",")
				      cen = PPoint(float(params[0].split('(')[1]),
						  float(params[1]))
				      ellA = float(params[2])
				      ellB = float(params[3])
				      ellPA = float(params[4].split(')')[0])
				      if ellA < ellB:
					    ellA, ellB = ellB, ellA
					    ellPA += 90
					    
				      if 'text={outer}' in line:
					outer = True
				      else:
					outer = False
					
				      ellipse_mask(cen, ellA, ellB, ellPA, img3, dimx, dimy, outer=outer)
			    coun = coun + 1
			    #print ' Done %.1f perc.' % (float(coun)*100./float(Number_of_lines)), line
			    if show_running==True:
			      show_complete(coun,Number_of_lines)
		else:
		    for line in f:
			    #print line
			    # Detect circle regions and fill them with the mask value
			    if 'circle(' in line:
				      '''
				      param = ((line.split('circle(')[1]).split(')')[0]).split(',')
				      a, b ,r  = int(float(param[0])/factor), int(float(param[1])/factor), int(float(param[2])/factor) 
				      y,x = np.ogrid[-b:img3.shape[0]-b, -a:img3.shape[1]-a]
				      mask = x*x + y*y <= r*r
				      img3[mask] = value
				      '''
				      params = line.split("(")[-1].split(",")[0]
				      cen = PPoint(float(line.split("(")[-1].split(",")[0])/factor,
						  float(line.split("(")[-1].split(",")[1])/factor)
				      ellA = float(line.split("(")[-1].split(",")[2].split(')')[0])
				      ellB = ellA
				      ellPA = 0.
				      ellipse_mask(cen, ellA, ellB, ellPA, img3, dimx, dimy)
			    # Detect polygon regions and add them to the path list
			    elif 'polygon(' in line:
				      param = map(float, ((line.split('polygon(')[1]).split(')')[0]).split(','))
				      param2 = [None]*(len(param)/2)
				      for i in xrange(0,len(param2)):	param2[i] = (int(param[2*i]),int(param[2*i+1])) 
				      param2.append(param2[0]/factor)
				      codes = []
				      codes.append(Path.MOVETO)
				      for i in xrange(1,len(param2)-1): codes.append(Path.LINETO)
				      codes.append(Path.CLOSEPOLY)
				      path = Path(param2, codes)
				      paths.append(path)
				      coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
				      X = []; Y = []
				      for k in range(0,len(coords)-1,2):
					      X.append(int(float(coords[k])))
					      Y.append(int(float(coords[k+1])))

				      borders.append([min(X),min(Y),max(X),max(Y)])    
			    elif "ellipse" in line:
				      params = line.split(",")
				      cen = PPoint(float(params[0].split('(')[1]),
						  float(params[1]))
				      ellA = float(params[2])
				      ellB = float(params[3])
				      ellPA = float(params[4].split(')')[0])
				      if ellA < ellB:
					    ellA, ellB = ellB, ellA
					    ellPA += 90
				      ellipse_mask(cen, ellA, ellB, ellPA, img3, dimx, dimy)
			    coun = coun + 1
			    #print ' Done %.1f perc.' % (float(coun)*100./float(Number_of_lines)), line
			    if show_running==True:
			      show_complete(coun,Number_of_lines)
		
		if obj_ellipse==None:
		      # Loop over the image and all polygons and fill them with the mask value	
		      nx, ny =img3.shape[1], img3.shape[0]
		      '''
		      for i, j, path in product(range(0,nx), range(0,ny), paths):
			      #print path
			      #exit()
			      img3[j][i]=value if path.contains_point((i,j)) else img3[j][i]
		      hdulist3.flush()
		      '''      
		      for kk in range(len(paths)):
			path = paths[kk]
			border = borders[kk]
			xmin = border[0]
			ymin = border[1]
			xmax = border[2]
			ymax = border[3]
			for i, j in product(range(xmin,xmax), range(ymin,ymax)):
			  try:
				img3[j-1][i-1]=value if path.contains_point((i,j)) else img3[j-1][i-1]
			  except:
			    zz=1
		      hdulist3.flush()
		      
		      if output_mask!=None:
			if mask_DN==None:
			  for i in range(dimy):
				  for k in range(dimx):
				    if img3[i,k] == value:
						  img2[i,k] = mask_value
				    else:
						  img2[i,k] = 0
			else:
			  for i in range(dimy):
				  for k in range(dimx):
				    if img3[i,k] == value:
						  img2[i,k] = mask_value
				    else:
						  img2[i,k] = 0
				    if img3[i,k] == mask_DN and mask_DN!=value:
				      img2[i,k] = mask_value
				      img3[i,k] = value


			hdulist2.flush()
		else:
		      hdulist_inp = pyfits.open(input_image)
		      img_inp = hdulist_inp[0].data
		      
		      [xc,yc,sma,smb,PA] = obj_ellipse
		      cen = PPoint(xc/factor,yc/factor)

		      ellipse_mask(cen, sma, smb, PA, img3, dimx, dimy, img_inp)

		      # Loop over the image and all polygons and fill them with the mask value	
		      nx, ny =img3.shape[1], img3.shape[0]
		      for kk in range(len(paths)):
			path = paths[kk]
			border = borders[kk]
			xmin = border[0]
			ymin = border[1]
			xmax = border[2]
			ymax = border[3]
			for i, j in product(range(xmin,xmax), range(ymin,ymax)):
			  try:
				img3[j-1][i-1]=value if path.contains_point((i,j)) else img3[j-1][i-1]
			  except:
			    zz=1
		      hdulist3.flush()

		      if output_mask!=None:
			if mask_DN==None:
			  for i in range(dimy):
				  for k in range(dimx):
				    if img3[i,k] == value:
						  img2[i,k] = mask_value
				    else:
						  img2[i,k] = 0
			else:
			  for i in range(dimy):
				  for k in range(dimx):
				    if img3[i,k] == value:
						  img2[i,k] = mask_value
				    else:
						  img2[i,k] = 0
				    if img3[i,k] == mask_DN and mask_DN!=value:
				      img2[i,k] = mask_value
				      img3[i,k] = value


			hdulist2.flush()