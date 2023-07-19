#!/usr/bin/python
# DESCRIPTION:
# Script to convert region mask to fits mask.
# Pixels with counts >0 are a mask in the output mask image.
# MINIMAL USAGE: python convert_reg_to_mask.py [input_image] [region_file]

# Import the necessary modules
from astropy.io import fits as pyfits
#import pyfits
import numpy as np
import math
import itertools
from scipy import ndimage
import sys
from itertools import product
from matplotlib.path import Path
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import os
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


def show_complete(i,N):
  percent = 100. * float(i) / float(N)
  sys.stdout.write("\r%2d%%" % percent)
  sys.stdout.flush()

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


def ellipse_annulus_mask(inner_ellipse, outer_ellipse, mask_value, xSize, ySize):
        #ySize,xSize = np.shape(mask)
        MASK = np.zeros(shape=(ySize,xSize))
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
        

        for x in range(xMin, xMax):
            for y in range(yMin, yMax):
                dFocus1_in = hypot(x-focus1_in.x, y-focus1_in.y)
                dFocus2_in = hypot(x-focus2_in.x, y-focus2_in.y)
                dPoint_in = dFocus1_in + dFocus2_in

                dFocus1_out = hypot(x-focus1_out.x, y-focus1_out.y)
                dFocus2_out = hypot(x-focus2_out.x, y-focus2_out.y)
                dPoint_out = dFocus1_out + dFocus2_out                
                if dPoint_out < dEll_out and dPoint_in > dEll_in:
                  try:
                    MASK[y,x] = mask_value
                  except:
                    zz=1
        return MASK

def read_annulus(annulus_file="annulus.reg"):
   ellipses = []
   for line in open(annulus_file):
        if "ellipse" in line:
          params = line.split(",")
          cen = Point(float(params[0].split('(')[1]),
		      float(params[1]))
          sma_in = float(params[2])
          smb_in = float(params[3])
          sma_out = float(params[4])
          smb_out = float(params[5])
          PA = float(params[6].split(')')[0])
   inner_ellipse = [cen,sma_in,smb_in,PA]
   outer_ellipse = [cen,sma_out,smb_out,PA]
   ann_width = abs(sma_out-sma_in)
   return inner_ellipse,outer_ellipse,ann_width


def mask(input_image, reg_file, output_mask=None, mask_value=1, mask_DN=None, verbosity=True):


                if verbosity: print('Masking contaminants...')

                    
                if output_mask is None:
                    output_mask = reg_file.split('.reg')[0]+'_masked_annulus.fits'  

                
                hdulist = pyfits.open(input_image)
                img = hdulist[0].data
                header = hdulist[0].header
                hdulist.close()

                (dimy,dimx) = img.shape
                img2 = np.zeros(np.shape(img))


                inner_ellipse,outer_ellipse,ann_width = read_annulus(reg_file)
                mask = ellipse_annulus_mask(inner_ellipse, outer_ellipse, mask_value, dimx, dimy)
 


                if verbosity: print('\nCreating mask image...')
                

                outHDU = pyfits.PrimaryHDU(mask, header=header)
                outHDU.writeto(output_mask, overwrite=True)                       
                    
                if verbosity: print('\nDone!')


            
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert region file with a mask to a fits segmentation image")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("reg_file", help="Input region file with an annulus ellipse mask")
    parser.add_argument("--output_mask", nargs='?', const=1, help="Optional: Output mask with masked values >0", type=str, default=None) 
    parser.add_argument("--mask_value", nargs='?', const=1, help="Optional: Masking value in output_mask", type=int, default=1) 
 
    args = parser.parse_args()

    input_image = args.input_image
    reg_file = args.reg_file
    output_mask = args.output_mask   
    mask_value = args.mask_value

    mask(input_image, reg_file, output_mask=output_mask, mask_value=mask_value, verbosity=True)
