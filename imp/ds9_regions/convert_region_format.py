#!/usr/bin/python
# -*- coding:  cp1251 -*-
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
import imp_setup
import os
import subprocess

import find_objects
import calibration


def crop_region_mask(image, input_region, output_region, x_l,y_l,x_r,y_r):
      hdulist = pyfits.open(image)
      inframe = hdulist[0].data
      ny,nx = np.shape(inframe)


      if int(math.floor(y_l-0.5))<0:
	  ymin = 0
	  print "ymin is set to 0"
      else:
	  ymin = int(math.floor(y_l-0.5))
      
      if int(math.floor(x_l-0.5))<0:
	  xmin = 0
	  print "xmin is set to 0"
      else:
	  xmin = int(math.floor(x_l-0.5))
      
      if int(math.floor(y_r-0.5))>=ny:
	  ymax = ny-1    
	  print "ymax is set to dim_y"
      else:
	  ymax = int(math.floor(y_r-0.5))
      
      if int(math.floor(x_r-0.5))>=nx:
	  xmax = nx-1
	  print "xmax is set to dim_x"
      else:
	  xmax = int(math.floor(x_r-0.5))



      f = open(input_region, 'r')
      ff = open(output_region, 'w')
      ff.write('image\n')
      for line in f:
			    if 'circle(' in line:
				      params = line.split("(")[-1].split(",")[0]
				      cen = [float(line.split("(")[-1].split(",")[0]),float(line.split("(")[-1].split(",")[1])]
				      ellA = float(line.split("(")[-1].split(",")[2].split(')')[0])
				      ellB = ellA
				      ellPA = 0.
				      [[x_min,y_min],[x_max,y_max]] = calibration.ellipse_borders([cen[0],cen[1]],ellA,ellB,ellPA)
				      if x_max>xmin and x_min<xmax and y_max>ymin and y_min<ymax:
					ff.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f)\n" % (cen[0]-xmin, cen[1]-ymin,
												    ellA, ellB, ellPA))
    
			    # Detect polygon regions and add them to the path list
			    elif 'polygon(' in line:
				      coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
				      X = []; Y = []
				      poly = False
				      for k in range(0,len(coords)-1,2):
					      X.append(float(coords[k]))
					      Y.append(float(coords[k+1]))
					      if float(coords[k])>x_l and float(coords[k])<x_r and float(coords[k+1])>y_l and float(coords[k+1])<y_r:
						poly=True
				      if poly == True:
					for i in range(len(X)):
					  if i==0:
					    ff.write('polygon(%.1f,%.1f,' % (X[i]-xmin,Y[i]-ymin) )
					  elif i>0 and i<len(X)-1:
					    ff.write('%.1f,%.1f,' % (X[i]-xmin,Y[i]-ymin) )
					  else:
					    ff.write('%.1f,%.1f)\n' % (X[i]-xmin,Y[i]-ymin) )

    
			    elif "ellipse" in line:
				      params = line[8:-2].split(",")
				      cen = [float(params[0]),float(params[1])]
				      ellA = float(params[2])
				      ellB = float(params[3])
				      try:
					ellPA = float(params[4])
				      except:
					ellPA = 0.
				      if ellA < ellB:
					  ellA, ellB = ellB, ellA
					  ellPA += 90
				      [[x_min,y_min],[x_max,y_max]] = calibration.ellipse_borders([cen[0],cen[1]],ellA,ellB,ellPA)
				      if x_max>xmin and x_min<xmax and y_max>ymin and y_min<ymax:
					ff.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f)\n" % (cen[0]-xmin, cen[1]-ymin,
												    ellA, ellB, ellPA))  

				    
	    
      f.close()
      ff.close()
      hdulist.close()









def wcs_to_image(image, input_region, output_region):
      p = subprocess.Popen(["ds9",image,"-scale","log","-invert","-regions","load",input_region,"-regions","system","image","-regions","save",output_region,"-exit"])
      p.wait()
      
      hdulist = pyfits.open(image)
      inframe = hdulist[0].data
      ny,nx = np.shape(inframe)
    
      f = open(output_region, 'r')
      ff = open('tmp.reg', 'w')
      ff.write('image\n')
      for line in f:
			    if 'circle(' in line:
				      params = line.split("(")[-1].split(",")[0]
				      cen = [float(line.split("(")[-1].split(",")[0]),float(line.split("(")[-1].split(",")[1])]
				      ellA = float(line.split("(")[-1].split(",")[2].split(')')[0])
				      ellB = ellA
				      ellPA = 0.
				      [[x_min,y_min],[x_max,y_max]] = calibration.ellipse_borders([cen[0],cen[1]],ellA,ellB,ellPA)
    
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
				      x_min = min(X); y_min = min(Y); x_max = max(X); y_max = max(Y)
				    
			    elif "ellipse" in line:
				      params = line[8:-2].split(",")
				      cen = [float(params[0]),float(params[1])]
				      ellA = float(params[2])
				      ellB = float(params[3])
				      try:
					ellPA = float(params[4])
				      except:
					ellPA = 0.
				      if ellA < ellB:
					  ellA, ellB = ellB, ellA
					  ellPA += 90
				      [[x_min,y_min],[x_max,y_max]] = calibration.ellipse_borders([cen[0],cen[1]],ellA,ellB,ellPA)
	  
  
			    if ('circle(' in line or 'polygon(' in line or "ellipse" in line) and x_max>=0 and y_max>=0 and x_min<nx and y_min<ny:
				ff.write(line)
	    
      f.close()
      ff.close()
      hdulist.close()
      shutil.move('tmp.reg', output_region)

      
def image_to_wcs(image, input_region, output_region):
      p = subprocess.Popen(["ds9",image,"-scale","log","-invert","-regions","load",input_region,"-regions","system","wcs","-regions","save",output_region,"-exit"])
      p.wait()  