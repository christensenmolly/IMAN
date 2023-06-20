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




def main(input_region_file, output_region_file=None, xc=None, yc=None):
      if output_region_file is None:
            output_region_file = input_region_file.split('.reg')[0]+'_object.reg'
    
      circle = Point(xc,yc)#.buffer(1)

      ff = open(output_region_file, 'w')
      ff.write('image\n')
      f = open(input_region_file, "r")
      for line in f:
            if 'polygon(' in line:
                coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
                pol = []
                for kk in range(0,len(coords)-1,2):
                   pol.append((int(float(coords[kk])),int(float(coords[kk+1]))))
                try:
                    polygon = Polygon(pol)
                    if circle.intersects(polygon):
                            ff.write(line)
                except:
                    z=1

      ff.close()
      f.close()    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Masking")
    parser.add_argument("input_region", help="Input region file")
    parser.add_argument("xc", help="Galaxy center x",type=float)
    parser.add_argument("yc", help="Galaxy center y",type=float)  
    
    parser.add_argument("--o", nargs='?', const=1, help="Optional: Output region file",type=str, default=None)  
    

    args = parser.parse_args()

    input_region_file = args.input_region
    output_region_file = args.o
    xc = args.xc
    yc = args.yc
    
    main(input_region_file, output_region_file=output_region_file, xc=xc, yc=yc)