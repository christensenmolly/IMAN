#!/usr/bin/env python
# -*- coding: utf8 -*-
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
import heapq

from photutils import datasets
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import EllipticalAperture



def stats(image,iso_region,mask_image='none',layer=0):
                inputname = image
                hdulist = pyfits.open(inputname)
                inframe = hdulist[layer].data
                inframe = np.array(inframe,float)
                
                
                iso_reg = open(iso_region,'r')
                paths=[]
                ySize, xSize = inframe.shape

                for line in iso_reg:
                        # Detect polygon regions and add them to the path list
                        if 'polygon(' in line:
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
                                

                                nx, ny =inframe.shape[1], inframe.shape[0]
                                INT = []
                                if mask_image=='none':
                                  for i, j, path in product(range(0,nx), range(0,ny), paths):
                                    if path.contains_point((i,j)):
                                            INT.append(inframe[j][i])
                                else:
                                  hdulist1 = pyfits.open(mask_image)
                                  mask = hdulist1[0].data
                                  for i, j, path in product(range(0,nx), range(0,ny), paths):
                                    if (path.contains_point((i,j)) and mask[j][i]==0.):
                                            INT.append(inframe[j][i])
                                iso_reg.close()
                                return sum(INT),np.mean(INT),np.median(INT),np.std(INT)
                      
                        if 'ellipse(' in line:
                                if mask_image!='none':
                                    hdulist_infr = pyfits.open(mask_image)
                                    mask_infr = hdulist_infr[0].data

                                    mask_astropy = np.zeros_like(mask_infr, dtype=bool)
                                    nx, ny =mask_infr.shape[1], mask_infr.shape[0]
                                    for k in range(ny):
                                      for i in range(nx):
                                        if mask_infr[k,i]>0:
                                            mask_astropy[k,i] = True
                                        else:
                                            mask_astropy[k,i] = False


                                params = line[8:-2].split(",")
                                xc = float(params[0])
                                yc = float(params[1])
                                ellA = float(params[2])
                                ellB = float(params[3])
                                try:
                                        ellPA = float(params[4])
                                except:
                                        ellPA = 0.
                                if ellA < ellB:
                                        ellA, ellB = ellB, ellA
                                        ellPA += 90 
                                
                                apertures = EllipticalAperture([(xc,yc)], ellA,ellB,radians(ellPA))
                                if mask_image!='none':
                                    phot_table = aperture_photometry(inframe, apertures,mask=mask_astropy)
                                else:
                                    phot_table = aperture_photometry(inframe, apertures)
                                Flux_ell_masked = float(phot_table['aperture_sum'])
                                iso_reg.close()
                                return Flux_ell_masked,0.,0.,0.

                
def Pixels_inside(image,iso_region,mask_image=None,layer=0):
                inputname = image
                hdulist = pyfits.open(inputname)
                inframe = hdulist[layer].data
                iso_reg = open(iso_region,'r')
                lines = iso_reg.readlines()
                paths=[]
                ySize, xSize = inframe.shape
                line = heapq.nlargest(2, list(lines), key=len)[0]

                if 'polygon(' in line:
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


                nx, ny =inframe.shape[1], inframe.shape[0]

                coords = []
                INT = []
                if mask_image==None:
                  for i, j, path in product(range(0,nx), range(0,ny), paths):
                    if path.contains_point((i,j)):
                            INT.append(inframe[j][i])
                else:
                  hdulist1 = pyfits.open(mask_image)
                  mask = hdulist1[0].data
                  for i, j, path in product(range(0,nx), range(0,ny), paths):
                    if (path.contains_point((i,j)) and mask[j][i]==0.):
                            coords.append(i)
                    else:

                            inframe[j][i] = 0.

                return inframe,coords
                
