#!/usr/bin/python
# DESCRIPTION:
# Script create region file with stars for solving astrometry.
# For the input and reference image, each star should be presented by a circle region with a label (text) which is the unique number of this star. The same stars in the reference and input image should have the same number.
# NOTE: After launching this script, the output region file can be used for add_wcs.py as region_file.
# MINIMAL USAGE: python create_astrom_stars.py [input_image] [reference_image] [region_file_image] [region_file_reference] 
# EXAMPLE: python create_astrom_stars.py HCG041_11Dec2018_ud.fits HCG41_dss.fits ima.reg dss.reg


# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import math
import sys
import shutil
import argparse
import os
import astroquery
from astroquery.simbad import Simbad
from astropy import coordinates
import astropy.units as u
from astropy import wcs
import collections
from astroquery.vizier import Vizier
from photutils import CircularAperture
from astropy.stats import sigma_clipped_stats
from photutils import centroid_2dg
import warnings
from photutils import CircularAnnulus, EllipticalAnnulus
from photutils import CircularAperture, EllipticalAperture
from photutils import aperture_photometry
warnings.filterwarnings("ignore")





def precise_center(data,nx,ny,x,y,R):
  ymin = max([0,int(math.floor(y-R))])
  ymax = min([ny,int(math.floor(y+R))])
  xmin = max([0,int(math.floor(x-R))])
  xmax = min([nx,int(math.floor(x+R))])
  x_cen, y_cen = centroid_2dg(data[ymin:ymax,xmin:xmax])
  x = x_cen+xmin; y = y_cen+ymin
  return x+1.,y+1.


def read_reg_file(region_file, data):
    ny,nx = np.shape(data)
    f = open(region_file, "r")

    lines = f.readlines()
    x = []; y = []; numb = []; rad = []
    N = 0
    for line in lines:
      if 'circle(' in line:
        N = N+1
        xc = float(line.split('(')[-1].split(',')[0])
        yc = float(line.split(',')[1])
        rad_in = float(line.split(',')[2].split(')')[0])

        xc,yc = precise_center(data,nx,ny,xc,yc,rad_in)
        x.append(xc)
        y.append(yc)
        rad.append(rad_in)
        if 'text=' in line:
            numb.append(int(line.split('text={')[-1].split('}')[0]))
        else:
            numb.append(N)
    return x,y,rad,numb
        



def main(input_image, reference_image, region_file_image, region_file_reference):
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    x_image,y_image,rad_image,numb_image = read_reg_file(region_file_image, data)
    
    hdulist = pyfits.open(reference_image)
    data = hdulist[0].data
    header = hdulist[0].header
    x_ref,y_ref,rad_ref,numb_ref = read_reg_file(region_file_reference, data)  

    ff = open('astrom_stars.reg', 'w')
    ff.write('global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    ff.write('image\n')
    w = wcs.WCS(header)
    for k in range(len(x_ref)):
        pix = np.array([[x_ref[k],y_ref[k]]], np.float_)
        world = w.wcs_pix2world(pix, 1)
        RA,DEC = world[0][0],world[0][1]
        ind = numb_image.index(numb_ref[k])
        ff.write('circle(%f,%f,%f) # text={%f,%f}\n' % (x_image[ind],y_image[ind],rad_image[ind],RA,DEC))
    ff.close()
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create region file with at least three good stars for astrometry")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("reference_image", help="Reference fits image")
    
    parser.add_argument("region_file_image", help="Region file with stars marked as circles and subscribed by a number corresponding to the same star in reference image")
    parser.add_argument("region_file_reference", help="Region file with stars marked as circles and subscribed by a number corresponding to the same star in input image") 

    args = parser.parse_args()

    input_image = args.input_image
    reference_image = args.reference_image
    region_file_image = args.region_file_image
    region_file_reference = args.region_file_reference


    main(input_image, reference_image, region_file_image, region_file_reference)



