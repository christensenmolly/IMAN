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

def main(input_file, output_file):
    hdulist1 = pyfits.open(input_file)
    data = hdulist1[0].data
    header = hdulist1[0].header
    new_data = data/np.sum(data)


    hdu = pyfits.PrimaryHDU(new_data, header)
    hdu.writeto(output_file, overwrite=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Norm image")
    parser.add_argument("input_file", help="Input  image")
    parser.add_argument("output_file", help="Output image")  
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file


    main(input_file, output_file)
