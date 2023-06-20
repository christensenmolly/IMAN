#!/usr/bin/python
# DESCRIPTION:
# Script to crop the input_image based on the given coordinates
# of the bottom left (xl,yl) and top right (xr, yr) corners. The crop will be
# saved in the output_image. This function keeps the WCS!
# MINIMAL USAGE: python crop_image.py [input_image] [xl yl xr yr]

# Import the necessary modules
import numpy as np
import math
from numpy import *
import sys
from astropy.stats import sigma_clipped_stats
import warnings
import argparse
warnings.filterwarnings("ignore")  

try:
  import astropy.io.fits as pyfits
  from astropy import wcs
except:
  warnings.warn("Astropy is not installed! No WCS has been added to the header!")


def ds9_to_np(x):
    '''
    Function to convert ds9 pixel coordinates to numpy format, i.e
    0.5 pix is 0, 1.5 pix is 1
    '''
    return int(round(x)) - 1


def main(input_image, x_l, y_l, x_r, y_r, output_image=None, hdu=0):
    print('Cropping image...')
    '''
    Function to crop the input_image based on the given coordinates
    of the bottom left and top right corners. The crop will be
    saved in the output_image.
    This function keeps the WCS!
    
    x_l, y_l, x_r, y_r - pixel coordinates in the ds9 format
    (start from 0.5).
    hdu - Operate on the specified FITS header extension (HDU), int.
    '''

    if output_image is None:
        output_image = input_image.split('.fits')[0] + '_crop.fits'
    
    # Convert to numpy coordinates
    x_l = ds9_to_np(x_l)
    y_l = ds9_to_np(y_l)
    x_r = ds9_to_np(x_r)
    y_r = ds9_to_np(y_r)
    
    #print(x_l,y_l,x_r,y_r)
    hdulist = pyfits.open(input_image)
    referenceheader = pyfits.getheader(input_image, hdu)
    inframe = hdulist[hdu].data
    ny,nx = np.shape(inframe)
    header = referenceheader
    
    
    outframe = inframe[y_l:y_r+1,x_l:x_r+1] # NOTE: +1pix - including the right coords.
    
    try:
        header['CRPIX1'] = header['CRPIX1'] - x_l
        header['CRPIX2'] = header['CRPIX2'] - y_l
    except:
        zz=1
        print('No WCS found!')
    hdu = pyfits.PrimaryHDU(outframe, header)
    hdu.writeto(output_image, overwrite=True)
    
    ny_new,nx_new = np.shape(outframe)

    
    print('Done!')
    return nx_new/2., ny_new/2.

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Cropping")
    parser.add_argument("input_image", help="Input image")    
    parser.add_argument("--corners", help="Corners (left and right, all separated by comma): xl,yl,xr,yr", type=str, default=None)
    parser.add_argument("--output_image", help="Optional: Output image", type=str, default=None) 
    parser.add_argument("--hdu", help="Optional: HDU layer", type=int, default=0)

    args = parser.parse_args()

    input_image = args.input_image
    corners = args.corners
    output_image = args.output_image
    hdu = args.hdu
    
    [x_l, y_l, x_r, y_r,] = corners.split(',')

    
    main(input_image, float(x_l), float(y_l), float(x_r), float(y_r), output_image=output_image, hdu=hdu)

