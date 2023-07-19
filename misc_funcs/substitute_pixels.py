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

def main(first_file, second_file, ranges, output_file):
    hdulist1 = pyfits.open(first_file)
    frame1 = hdulist1[0].data
    header = hdulist1[0].header


    hdulist2 = pyfits.open(second_file)
    frame2 = hdulist2[0].data

    sframe = np.copy(frame1)
    
    [xmin,ymin,xmax,ymax] = list(np.array(ranges.split(','), int))

    sframe[ymin-1:ymax-1,xmin-1:xmax-1] = frame2[ymin-1:ymax-1,xmin-1:xmax-1]

    hdu = pyfits.PrimaryHDU(sframe, header)
    hdu.writeto(output_file, overwrite=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Substitution of pixels in an image from another image")
    parser.add_argument("first_file", help="Input first fits image where pixels will be substututed")
    parser.add_argument("second_file", help="Input second fits image by which you will substutute the pixels in the first image")  
    parser.add_argument("ranges", help="Input the range in the format xmin,ymin,xmax,ymax")  
    parser.add_argument("output_file", help="Input the name of the output file")  
    args = parser.parse_args()

    first_file = args.first_file
    second_file = args.second_file
    ranges = args.ranges
    output_file = args.output_file

    main(first_file, second_file, ranges, output_file)
