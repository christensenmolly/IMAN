# Import the necessary modules
from astropy.io import fits as pyfits
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
import subprocess

import convert_reg_to_mask
from astropy import wcs



        
    
    


def main(input_image, coordinates=None, output_region='objects_to_unmask.reg', coord_system='pix'):
    hdulist = pyfits.open(input_image)
    img = hdulist[0].data 
    header = hdulist[0].header 
    w = wcs.WCS(header)
    
    if coordinates is None:
        print('No coordinates are given!')
        return 0

    f = open(output_region, 'w')
    f.write('global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write('image\n')

    for k in range(len(coordinates)):
        if coord_system=='wcs':
            RA = float(coordinates[k][0])
            DEC = float(coordinates[k][1])
            
            world = np.array([[RA,DEC]], np.float_)
            pixcrd = w.wcs_world2pix(world, 1)
            xc,yc = pixcrd[0][0],pixcrd[0][1]
        else:
            xc = float(coordinates[k][0])
            yc = float(coordinates[k][1])

        f.write('point(%.3f,%.3f) # point=x\n' % (xc,yc))

    f.close()
    hdulist.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unmasking targets")
    parser.add_argument("inputImage", help="Input fits image with the object")
    parser.add_argument("coordinates", help="Input coordinates of objects to unmask, separated by comma and / for individual galaxies")
    parser.add_argument("--region", nargs='?', const=1, help="Output region file with the objects to unmask", type=str, default='objects_to_unmask.reg')
    parser.add_argument("--coord_system", nargs='?', const=1, help="Optional: System of coordinates (pix or wcs)", type=str, default='pix') 

    args = parser.parse_args()

    input_image = args.inputImage
    coordinates = args.coordinates
    output_region = args.region 
    coord_system = args.coord_system
    
    coord_galaxies = coordinates.split('/')
    right_coordinates = []
    for k in range(len(coord_galaxies)):
        x,y = coord_galaxies[k].split(',')
        right_coordinates.append([float(x),float(y)])

    main(input_image, right_coordinates, output_region, coord_system)
        
        
        
        
