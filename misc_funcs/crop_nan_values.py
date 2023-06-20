#!/usr/bin/env python
# -*- coding: utf8 -*-

# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import itertools
import matplotlib.pyplot as plt
import sys
import math
import argparse
import os

LOCAL_DIR = "/misc_funcs"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))

import crop_image



def main(input_file, output_file):
    hdulist1 = pyfits.open(input_file)
    data = hdulist1[0].data
    header = hdulist1[0].header
    ny,nx = np.shape(data)
    
    inds = np.argwhere(~np.isnan(data))
    
    inds = np.transpose(inds)
    ymin = np.min(inds[0])
    xmin = np.min(inds[1])
    ymax = np.max(inds[0])
    xmax= np.max(inds[1])
    
    crop_image.main(input_file, xmin, ymin, xmax, ymax, output_image=output_file, hdu=0)
    return xmin,ymin,xmax,ymax

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mask values")
    parser.add_argument("input_file", help="Input  image")
    parser.add_argument("output_file", help="Output image")  
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file


    main(input_file, output_file)