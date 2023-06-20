#!/usr/bin/python

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
import argparse
#import imp_setup
import os
import subprocess
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats
#from photutils import data_properties, properties_table
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
try:
    import find_objects
    import calibration
except:
    z=1

import mask_indiv

class PPoint:
    def __init__(self, x, y):
        self.x = x
        self.y = y    
    
# Main function, this part actually runs when routine is called
def ellipse_mask(cen, ellA, ellB, ellPA, inframe, xSize, ySize, img_inp=None,img_mask=None, outer = False):
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


def unmask_galaxy(segm, xc=None, yc=None, min_radius=10.):
    # Find the masked area which corresponds to the galaxy and unmask it. All objects within at least 2*min_radius will be unmasked.
    
    segm_without_galaxy = np.copy(segm)
    
    
    ny,nx = np.shape(segm_without_galaxy)
    
    segm_only_galaxy = np.zeros(shape=(ny,nx))
    
    if xc is None or yc is None:
        xc = nx/2.
        yc = ny/2.
    
    labels = []
    for y in range(int(yc-2.*min_radius), int(yc+2.*min_radius)):
        for x in range(int(xc-2.*min_radius), int(xc+2.*min_radius)):
            try:
                lable = segm_without_galaxy[y,x]
                if lable not in labels:
                    labels.append(lable)
            except:
                z=1 # Beyond the image borders
    
    for label in labels:
        segm_without_galaxy[segm_without_galaxy == label] = 0.
        segm_only_galaxy[segm == label] = 1.
        
    return segm_without_galaxy, segm_only_galaxy  


def add_galaxy_mask(input_mask_image, output_mask_image, galaxy_ellipse, sma_in_factor=1.250, sma_out_factor=1.601):
    # galaxy_ellipse = [ [ [xCen, yCen], ellA, ellB, PA], [...] ]
    hdulist = pyfits.open(input_mask_image)
    img_mask = hdulist[0].data
    ny,nx = np.shape(img_mask)

    
    for k in range(len(galaxy_ellipse)):
        [[xCen, yCen], ellA, ellB, ellPA] = galaxy_ellipse[k]
        cen = PPoint(xCen, yCen)
        if sma_in_factor==0.:
            # Mask within an ellipse with sma = sma_out_factor * ellA
            mask_indiv.ellipse_mask1(cen, ellA*sma_out_factor, ellB*sma_out_factor, ellPA, img_mask, nx, ny, img_inp=None, img_mask=img_mask, outer = False)
        else:
            # Only area inside an elliptical annulus is non-masked:
            
            # Mask within an ellipse with sma = sma_in_factor * ellA
            #mask_indiv.ellipse_mask(cen, ellA*sma_in_factor, ellB*sma_in_factor, ellPA, img_mask, nx, ny, img_inp=None, img_mask=img_mask, outer = False)            
            # Mask outside an ellipse with sma = sma_out_factor * ellA
            #mask_indiv.ellipse_mask(cen, ellA*sma_out_factor, ellB*sma_out_factor, ellPA, img_mask, nx, ny, img_inp=None, img_mask=img_mask, outer = True)  
            
            #img_mask = mask_indiv.ellipse_annulus_mask([cen, ellA*sma_in_factor, ellB*sma_in_factor, ellPA], [cen, ellA*sma_out_factor, ellB*sma_out_factor, ellPA], img_mask)
            img_mask = mask_indiv.ellipse_annulus_mask([cen, ellA*sma_in_factor, ellB+ellA*(sma_in_factor-1), ellPA], [cen, ellA*sma_out_factor, ellB+ellA*(sma_out_factor-1), ellPA], img_mask)
            
    hdu = pyfits.PrimaryHDU(img_mask)
    hdu.writeto(output_mask_image, clobber=True)       
    
    return img_mask
