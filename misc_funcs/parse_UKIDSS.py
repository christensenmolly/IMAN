#! /usr/bin/env python

import sys
from os.path import exists
from os import remove
from scipy.ndimage.filters import convolve
from numpy import max, sum
import astropy.io.fits as pyfits

def main(input_image, output_image):
    hdulist = pyfits.open(input_image)

   
    image = hdulist[1].data
    header = hdulist[1].header

    hdu = pyfits.PrimaryHDU(data=image, header=header)
    hdu.writeto(output_image, overwrite=True)


if __name__ == "__main__":
    input_image = sys.argv[1]
    output_image = sys.argv[2]

    main(input_image, output_image)

    
