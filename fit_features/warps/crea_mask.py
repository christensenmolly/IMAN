#!/usr/bin/python
# -*- coding:  cp1251 -*-
# Import the necessary modules
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



def mask(input_image, line):
  if 'polygon(' in line:
        x = []
        y = []
        coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
        for k in range(0,len(coords)-1,2):
          x.append(float(coords[k]))
          y.append(float(coords[k+1]))
        #plt.plot(x,y,'o')
        #plt.show()
        #print line
        hdulist = pyfits.open(input_image)
        img = hdulist[0].data
        
        img2 = np.copy(img)
        img3 = np.copy(img)


        (dimy,dimx) = img3.shape
        value = 0.
        paths=[]


        param = map(float, ((line.split('polygon(')[1]).split(')')[0]).split(','))
        param = list(param)
        param2 = [None]*(int(len(param)/2))
        for i in range(0,len(param2)):    param2[i] = (int(param[2*i]),int(param[2*i+1])) 
        param2.append(param2[0])
        codes = []
        codes.append(Path.MOVETO)
        for i in range(1,len(param2)-1): codes.append(Path.LINETO)
        codes.append(Path.CLOSEPOLY)
        path = Path(param2, codes)
        paths.append(path)

        
        # Loop over the image and all polygons and fill them with the mask value    
        nx, ny =img3.shape[1], img3.shape[0]
        for i, j, path in product(range(0,nx), range(0,ny), paths):
            img3[j][i]=value if path.contains_point((i,j)) else img3[j][i]
        

        for i in range(dimy):
            for k in range(dimx):
              if img3[i,k] == value:
                    img2[i,k] = 1
              else:
                    img2[i,k] = 0


        return img2
    