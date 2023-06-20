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

def main(input_file, output_file, value):
    hdulist1 = pyfits.open(input_file)
    data = hdulist1[0].data
    header = hdulist1[0].header
    ny,nx = np.shape(data)
    for k in range(ny):
      for i in range(nx):
        if np.isnan(data[k,i]):
            data[k,i]=value


    hdu = pyfits.PrimaryHDU(data, header)
    hdu.writeto(output_file,clobber=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arithmetical operations")
    parser.add_argument("input_file", help="Input  image")
    parser.add_argument("output_file", help="Output image")  
    parser.add_argument("value", help="Input value to replace with")  
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    value = float(args.value)


    main(input_file, output_file, value)
