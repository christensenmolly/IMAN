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
import subprocess


def main(input_image, segm_image, out_region_file = None):
    if out_region_file is not None:
        if not os.path.exists(out_region_file):
            open(out_region_file, 'a').close()
    
    hdulist = pyfits.open(input_image)
    img = hdulist[0].data 
    header = hdulist[0].header 
    
    hdulist1 = pyfits.open(segm_image)
    segm = hdulist1[0].data

    
    '''
    outHDU = pyfits.PrimaryHDU(img/np.sum(img), header=header)
    outHDU.writeto('ima_display_red.fits', overwrite=True)     

    outHDU = pyfits.PrimaryHDU(np.log10(segm), header=header)
    outHDU.writeto('segm_display_green.fits', overwrite=True)  

    outHDU = pyfits.PrimaryHDU(segm, header=header)
    outHDU.writeto('segm_display_blue.fits', overwrite=True)  

    

    subprocess.call("ds9 -rgb -red %s -green %s -blue %s -rgb view blue no -rgb red -scale histequ" % ('ima_display_red.fits','segm_display_green.fits','segm_display_blue.fits'), shell=True)
    '''


    outHDU = pyfits.PrimaryHDU(img/np.nansum(img), header=header)
    outHDU.writeto('ima_display_red.fits', overwrite=True)     

    outHDU = pyfits.PrimaryHDU(np.log10(segm), header=header)
    outHDU.writeto('segm_display_blue.fits', overwrite=True)  

    outHDU = pyfits.PrimaryHDU(segm, header=header)
    outHDU.writeto('segm_display_green.fits', overwrite=True)  

    

    #subprocess.call("ds9 -rgb -red %s -green %s -blue %s -rgb view green no -rgb red -scale histequ" % ('ima_display_red.fits','segm_display_green.fits','segm_display_blue.fits'), shell=True)

    subprocess.call("ds9 -rgb -red %s -green %s -blue %s -rgb red -scale histequ -region %s" % ('ima_display_red.fits','segm_display_green.fits','segm_display_blue.fits',out_region_file), shell=True)
    
    os.remove('ima_display_red.fits')
    os.remove('segm_display_green.fits')
    os.remove('segm_display_blue.fits')
    hdulist.close()
    hdulist1.close()
    
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky background estimation")
    parser.add_argument("inputImage", help="Input fits image with the object")
    parser.add_argument("segmImage", help="Input segmenation image")
    parser.add_argument("--region", nargs='?', const=1, help="Optional: Output region file where wrong and new regions are specified.", type=str, default=None) 
    #parser.add_argument("--outputImage", nargs='?', const=1, help="Optional: Output image with subtracted sky", type=str, default=None) 
    #parser.add_argument("--degree", nargs='?', const=1, help="Optional: Input the polinomial order to fit the sky within the annulus. NOT WORKING PROPERLY!", type=int, default=0)
    #parser.add_argument("--manual", nargs='?', const=1, help="Optional: Do you want to manually change the annulus in DS9", type=bool, default=False)
    #parser.add_argument("--annulus_width", nargs='?', const=1, help="The width of the annulus in pix.", type=float, default=30)
    #parser.add_argument("--galaxy_ellipse", nargs='?', const=1, help="Galaxy ellipse, either as a DS9 region file with an ellipse region, or in the format xc,yc,sma,smb,PA", type=str, default=None)
    args = parser.parse_args()

    input_image = args.inputImage
    segm_image = args.segmImage
    region = args.region

    main(input_image, segm_image, out_region_file=region )
