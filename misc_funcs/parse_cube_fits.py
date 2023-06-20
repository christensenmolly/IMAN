#! /usr/bin/env python

import sys
from os.path import exists
from os import remove
from scipy.ndimage.filters import convolve
from numpy import max, sum
import numpy as np
import astropy.io.fits as pyfits



def main(input_image,output_image,layer, header_layer=0,layer_desc='standard'):
    # Usually layer=0, layer_desc='standard'
    
    hdulist = pyfits.open(input_image)

    if layer_desc=='standard':
        image = hdulist[0].data
        header = hdulist[header_layer].header

        hdu = pyfits.PrimaryHDU(data=image[layer], header=header)
        hdu.writeto(output_image,overwrite=True)
    else:
        image = hdulist[layer].data
        header = hdulist[header_layer].header

        hdu = pyfits.PrimaryHDU(data=image, header=header)
        hdu.writeto(output_image,overwrite=True)


if __name__ == "__main__":
    input_image = sys.argv[1]
    output_image = sys.argv[2]
    layer = int(sys.argv[3])
    header_layer = int(sys.argv[4])
    layer_desc = str(sys.argv[5])
    main(input_image,output_image,layer,header_layer,layer_desc)

    
