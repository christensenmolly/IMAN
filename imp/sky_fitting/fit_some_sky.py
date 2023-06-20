#!/usr/bin/python
# DESCRIPTION:
# Script to fit sky with polygons specified in a region file.
# MINIMAL USAGE: fit_some_sky.py [input_image] [region_file]

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
from astropy.stats import sigma_clipped_stats

import determine_sky

LOCAL_DIR = "/imp/sky_fitting"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))

import arithm_operations


def main(input_image, region_file, output_image='', sky_image='sky.fits', mask_image=None, polynomial_degree=0):
    print('Fitting some sky...')
    if output_image=='':
        output_image = input_image.split('.fits')[0] + '_someskysub.fits'
    
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    header = hdulist[0].header
    hdulist.close()
    
    Data = np.copy(data)
    
    if mask_image is not None:
        hdulist = pyfits.open(mask_image)
        mask = hdulist[0].data
        hdulist.close()
    else:
        mask = np.zeros(np.shape(data))
    
    Sky = np.zeros(np.shape(data))
    
                
    f = open(region_file, "r")
    N = 0
    for line in f:
        # Detect polygon regions and add them to the path list
        if 'polygon(' in line:
            mask_sky = np.ones(np.shape(data))
            
            print('\tPolygon #%i...' % (N+1))
            param = map(float, ((line.split('polygon(')[1]).split(')')[0]).split(','))
            param = list(param)
            param2 = [None]*(int(len(param)/2))
            for i in range(0,len(param2)):  param2[i] = (int(param[2*i]),int(param[2*i+1])) 
            param2.append(param2[0])
            codes = []
            codes.append(Path.MOVETO)
            for i in range(1,len(param2)-1):  codes.append(Path.LINETO)
            codes.append(Path.CLOSEPOLY)
            path = Path(param2, codes)
            coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
            if 'text=' in line:
                polynomial_degree = int(line.split('text=')[-1].rstrip().split('{')[-1].split('}')[0])
            print('\t\tPolynomial degree %i' % (polynomial_degree))
            X = []; Y = []
            for k in range(0,len(coords)-1,2):
               X.append(int(float(coords[k])))
               Y.append(int(float(coords[k+1])))

            xmin = min(X)
            ymin = min(Y)
            xmax = max(X)
            ymax = max(Y)
            for i, j in product(range(xmin,xmax), range(ymin,ymax)):
                try:
                    if path.contains_point((i,j)):
                        mask_sky[j-1][i-1] = 0
                except:
                    zz=1
            
            total_mask = mask_sky + mask
            mask_astropy = determine_sky.convert_segm_to_boolean(total_mask)

            if polynomial_degree==0:
                mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask_astropy)

                sky = mean#median

                for i, j in product(range(xmin,xmax), range(ymin,ymax)):
                    try:
                        if path.contains_point((i,j)):
                            #Data[j-1][i-1] = data[j-1][i-1] - sky
                            Sky[j-1][i-1] = sky
                    except:
                        zz=1 
            else:
                sky = determine_sky.fit_polynomial(data, int(polynomial_degree), mask=mask_astropy)

                for i, j in product(range(xmin,xmax), range(ymin,ymax)):
                    try:
                        if path.contains_point((i,j)):
                            #Data[j-1][i-1] = data[j-1][i-1] - sky[j-1][i-1]
                            Sky[j-1][i-1] = sky[j-1][i-1]
                    except:
                        zz=1            
            N = N+1
    # Save the output images
    #outHDU = pyfits.PrimaryHDU(Data, header=header)
    #outHDU.writeto(output_image, clobber=True)
      
    #if sky_image is not None:
    outHDU = pyfits.PrimaryHDU(Sky, header=header)
    outHDU.writeto(sky_image, clobber=True)
    
    arithm_operations.main(input_image, sky_image, 'sub', output_image)

    print('Done!')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky fitting within polygons")
    parser.add_argument("inputImage", help="Input fits image ")
    parser.add_argument("inputRegion", help="Optional: Output region file", type=str, default='general_mask.reg') 
    parser.add_argument("--inputMask", help="Input mask image", type=str, default=None) 
    
    parser.add_argument("--n", help="Optional: Polynomial degree (0 is constant level)", type=int, default=0) 
    parser.add_argument("--output_image", help="Optional: Output image", type=str, default='') 
    parser.add_argument("--output_sky", help="Optional: Output sky", type=str, default='sky.fits')
    
    args = parser.parse_args()

    input_image = args.inputImage
    region_file = args.inputRegion
    output_image = args.output_image
    sky_image = args.output_sky
    mask_image = args.inputMask
    polynomial_degree = args.n


    main(input_image, region_file, output_image=output_image, sky_image=sky_image, mask_image=mask_image, polynomial_degree=polynomial_degree)
