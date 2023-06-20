#!/usr/bin/python
# DESCRIPTION:
# Script to mask pixels which have intensities larger or lower than SBlevel.
# MINIMAL USAGE: python mask_SBlevel.py [input_image] [SBlevel]

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



import convert_segm_to_region



def main(input_image, SB_level, output_mask='SB_mask.fits', mask='>'):
        hdulist = pyfits.open(input_image)
        data = hdulist[0].data
        header = hdulist[0].header
        
        ySize, xSize = data.shape
        data_mask = np.zeros((ySize, xSize))
        
        if mask=='>':
            for k in range(ySize):
                for i in range(xSize):
                    if data[k,i]>SB_level:
                        data_mask[k,i]=1.
        else:
            for k in range(ySize):
                for i in range(xSize):
                    if data[k,i]<SB_level:
                        data_mask[k,i]=1. 


        pyfits.PrimaryHDU(data=data_mask, header=header).writeto(output_mask, overwrite=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky background estimation")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("SB_level", help="SB level in ADU",type=float) 
    parser.add_argument("--output_mask", nargs='?', const=1, help="Optional: Output mask image",type=str,default='SB_mask.fits')
    parser.add_argument("--mask", nargs='?', const=1, help="Optional: Mask pixels which are > (or <) than SBlevel", type=str, default='>')
    
    args = parser.parse_args()

    input_image = args.inputImage
    SB_level = args.SB_level
    output_mask = args.output_mask
    mask = args.mask

    main(input_image, SB_level, output_mask=output_mask, mask=mask)
    
