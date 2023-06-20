#!/usr/bin/env python
# -*- coding: utf8 -*-

# Import the necessary modules
import pyfits
import numpy as np
import itertools
import matplotlib.pyplot as plt
import sys
import math
import argparse
import os
from photutils import CosineBellWindow, create_matching_kernel
from photutils import create_matching_kernel, TopHatWindow
import convolve

def main(input_image1, psf_image1, psf_image2, output_image):
    #window = TopHatWindow(0.35)
    window = CosineBellWindow(alpha=0.35)
    
    imageHDU1 = pyfits.open(psf_image1)
    image1 = imageHDU1[0].data
    header1 = imageHDU1[0].header
    ny1,nx1 = np.shape(image1)  
    
    imageHDU2 = pyfits.open(psf_image2)
    image2 = imageHDU1[0].data
    header2 = imageHDU1[0].header
    ny2,nx2 = np.shape(image2) 
    
    if nx1!=nx2 or ny1!=ny2:
        print('Error! The kernels have different dimensions! Exiting...')
        exit()
    
    
    kernel = create_matching_kernel(image1/np.sum(image1), image2/np.sum(image2), window)
    
    outHDU = pyfits.PrimaryHDU(kernel)
    outHDU.writeto('match_psf.fits', clobber=True)   
    
    convolve.convolution(input_image1, 'match_kernel.fits', output_image)
    #os.remove('match_kernel.fits')
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])