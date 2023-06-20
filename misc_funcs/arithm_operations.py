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

def main(first_file, second_file, operation, output_file):
    hdulist1 = pyfits.open(first_file)
    frame1 = hdulist1[0].data
    header = hdulist1[0].header

    try:
      frame2 = float(second_file)
    except:
      hdulist2 = pyfits.open(second_file)
      frame2 = hdulist2[0].data

    sframe = np.copy(frame1)

    if operation=='add':
      sframe = frame1 + frame2
    elif operation=='sub':
      sframe = frame1 - frame2  
    elif operation=='mul':
      sframe = frame1 * frame2    
    elif operation=='div':
      sframe = frame1 / frame2  
    elif operation=='norm':
      sframe = sframe/np.sum(sframe)

    hdu = pyfits.PrimaryHDU(sframe,header)
    hdu.writeto(output_file,overwrite=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arithmetical operations")
    parser.add_argument("first_file", help="Input first fits image")
    parser.add_argument("second_file", help="Input second fits image or the value")  
    parser.add_argument("operation", help="Input name of the operation: add, sub, mul, div, norm")  
    parser.add_argument("output_file", help="Input the name of the output file")  
    args = parser.parse_args()

    first_file = args.first_file
    second_file = args.second_file
    operation = args.operation
    output_file = args.output_file

    main(first_file, second_file, operation, output_file)
