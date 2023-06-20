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


LOCAL_DIR = "/imp/masking"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'sextractor'))
import run_SExtractor

import convert_segm_to_region




def read_region(reg_file):
    f = open(reg_file, "r")    
    for line in f:
        if "ellipse" in line:
            params = line.split(",")
            cen = [float(params[0].split('(')[1]),float(params[1])]
            ellA = float(params[2])
            ellB = float(params[3])
            ellPA = float(params[4].split(')')[0])
            if ellA < ellB:
                ellA, ellB = ellB, ellA
                ellPA += 90    
            break
    f.close()
    return cen[0],cen[1],ellA,ellB,ellPA

def main(input_image, output_region_file='general_mask.reg', snr=2., min_pix=5, region_type='polygon', sextr_setup=None, galaxy_ellipse=None, offset_size=1., offset_pix=0., verbosity=True):
    """
    Function creates a region file with polygons that cover
    all objects (the galaxy can be removed from masking if
    galaxy_ellipse is specified).
    
    galaxy_ellipse = [[xc,yc], sma, smb, PA] (the outer ellipse) OR shapely polygon
    """
    if galaxy_ellipse is not None:
        try:
            if galaxy_ellipse.geom_type == 'Polygon':
                cont = True
                ellipse = galaxy_ellipse
                SMB = None; PA = None
                [(xc,yc)] = galaxy_ellipse.centroid.coords

            else:
                cont = False
                if verbosity: print('ERROR: Some error with the galaxy polygon provided! Exiting.')
                exit()
        except:
            try:
                xc,yc,SMA,SMB,PA = read_region(galaxy_ellipse)
            except:
                [[xc,yc],SMA,SMB,PA] = galaxy_ellipse
            circle = Point(xc,yc).buffer(1)
            if SMA is not None:
                if SMB is not None and PA is not None:
                    ellipse = shapely.affinity.scale(circle,SMA,SMB)
                    rot_ellipse = shapely.affinity.rotate(ellipse, PA, origin='center', use_radians=False)
                else:
                    ellipse = shapely.affinity.scale(circle,SMA,SMA)
    else:
        xc = None
        yc = None
        PA = None
        
    run_SExtractor.call_SE(input_image, snr=snr, min_pix=min_pix, verbosity=verbosity)

    # find SE number of our galaxy (it is in the centre of the image now
    # so we can just get the intensity of the central pixel)

    if galaxy_ellipse is None:
        segmHDU = pyfits.open("segm.fits")
        segmData = segmHDU[0].data
        ySize, xSize = segmData.shape
        xCen = int(xSize / 2.0)
        yCen = int(ySize / 2.0)
        galN = int(segmData[yCen, xCen])
        fout = open(output_region_file, "w")
        fout.truncate(0)
        fout.write("# Region file format: DS9 version 4.1\n")
        fout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 ')
        fout.write('normal roman" select=1 highlite=1 dash=0 fixed=0 ')
        fout.write('edit=1 move=1 delete=1 include=1 source=1\n')
        fout.write('image\n')
        for line in open("field.cat"):
            if line.startswith("#"):
                continue
            params = line.split()
            n = int(params[0])
            if n == galN:
                continue
            xCen = float(params[1])
            yCen = float(params[2])
            kron = float(params[8])
            ellA = kron * float(params[4]) * offset_size + offset_pix
            ellB = ellA * float(params[5])/float(params[4])#kron * float(params[5]) * offset_size
            ellPA = float(params[6])
            fout.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f)\n" % (xCen, yCen,
								    ellA, ellB,
								    ellPA))
        fout.close()
    
    if region_type=='polygon':
      # Convert segm to region
      #convert_segm_to_region.main('segm.fits', 'mask_tmp.reg')
      convert_segm_to_region.main('segm.fits', 'mask_tmp.reg', output_mask_image=None, fits_slice = 0, offset_size=offset_size, offset_pix=offset_pix, xc=None, yc=None, system='image', ignore_value=None)

      ff = open(output_region_file, 'w')
      ff.write('image\n')
      f = open('mask_tmp.reg', "r")
      for line in f:
            if 'polygon(' in line:
                coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
                pol = []
                for kk in range(0,len(coords)-1,2):
                   pol.append((int(float(coords[kk])),int(float(coords[kk+1]))))
                polygon = Polygon(pol)
                if PA is not None:
                    if not rot_ellipse.intersects(polygon):
                        ff.write(line)
                elif PA is None and xc is not None:
                    if not ellipse.intersects(polygon):
                        ff.write(line)
                else:
                    ff.write(line)
      ff.close()
      f.close()
    else:
      fout = open(output_region_file, "w")
      fout.truncate(0)
      fout.write("# Region file format: DS9 version 4.1\n")
      fout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 ')
      fout.write('normal roman" select=1 highlite=1 dash=0 fixed=0 ')
      fout.write('edit=1 move=1 delete=1 include=1 source=1\n')
      fout.write('image\n')
      for line in open("field.cat"):
           if line.startswith("#"):
               continue
           params = line.split()
           xCen = float(params[1])
           yCen = float(params[2])
           kron = float(params[8])
           ellA = kron * float(params[4]) * offset_size + offset_pix
           ellB = ellA * float(params[5])/float(params[4])#kron * float(params[5]) * offset_size
           ellPA = float(params[6])

           intersect = False

           if galaxy_ellipse is not None:
                cont_circle = Point(xCen,yCen).buffer(1)
                cont_ellipse = shapely.affinity.scale(cont_circle,ellA,ellB)
                cont_rot_ellipse = shapely.affinity.rotate(cont_ellipse, ellPA, origin='center', use_radians=False)
                
                if PA is not None:
                    if rot_ellipse.intersects(cont_rot_ellipse):
                        intersect = True
                else:
                    if ellipse.intersects(cont_rot_ellipse):
                        intersect = True
                        
           if intersect==False:
                fout.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f)\n" % (xCen, yCen,
								    ellA, ellB,
								    ellPA))
      fout.close()
      os.remove('field.cat')
      os.remove('segm.fits')
      
      if verbosity: print('Done!')
      return output_region_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky background estimation")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("--outputRegion", help="Optional: Output region file",type=str,default='general_mask.reg') 
    parser.add_argument("--region_type", nargs='?', const=1, help="Optional: Region type",type=str,default='ellipse')
    parser.add_argument("--snr", nargs='?', const=1, help="Optional: Input signal-to-noise ratio of pixels to mask",type=float,default=2.)
    parser.add_argument("--min_pix", nargs='?', const=1, help="Optional: Number of joint pixels",type=int,default=5)
    parser.add_argument("--offset_pix", nargs='?', const=1, help="Optional: Offset to make the regions larger or smaller, in pixels. Default 0..",type=float, default=0.) 
    parser.add_argument("--offset_size", nargs='?', const=1, help="Optional: Offset to make the regions larger or smaller, in units of the region size (multiplication). Default 1.",type=float, default=1.) 
    
    args = parser.parse_args()

    input_image = args.inputImage
    output_region_file = args.outputRegion
    region_type = args.region_type
    snr = args.snr
    min_pix = args.min_pix
    offset_pix = args.offset_pix
    offset_size = args.offset_size

    main(input_image, output_region_file, snr=snr, min_pix=min_pix, region_type=region_type, sextr_setup=None, galaxy_ellipse=None, offset_size=offset_size, offset_pix=offset_pix)
    
