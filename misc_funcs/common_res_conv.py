#! /usr/bin/env python
# Old module to generate a matching kernel between two PSFs
import os
import sys

from pylab import *
import astropy.io.fits as pyfits
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import glob
from matplotlib.pyplot import cm
import re
from scipy.odr.odrpack import *
import pyparsing
import matplotlib.pyplot as pyplot
import matplotlib as mpl
import numpy as np
import subprocess
import argparse
from scipy.ndimage.filters import convolve

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
sys.path.append(PATH_TO_PACKAGE+'/Functions')
import functions

def main(A_image,B_image,psf_A_image,psf_B_image):
    #taken from http://arxiv.org/pdf/1106.5065v2.pdf

    # Read the input image:
    imageHDU = pyfits.open(A_image)[0]
    image = imageHDU.data
    header = imageHDU.header
    
    # Read the input PSF files:
    psf_A = pyfits.open(psf_A_image)[0].data
    psf_A_ny, psf_A_nx = psf_A.shape
    psf_B = pyfits.open(psf_B_image)[0].data
    psf_B_ny, psf_B_nx = psf_B.shape


    # Normalize given kernels:
    psf_A = pyfits.open(psf_A_image)[0].data
    psf_A = psf_A / sum(psf_A)    

    psf_B = pyfits.open(psf_B_image)[0].data
    psf_B = psf_B / sum(psf_B)     
    
    # Find FT:
    FT_A = real(np.fft.fft2(psf_A))
    FT_B = real(np.fft.fft2(psf_B))
    
    # Create kernel for convolution from image A to B:
    K = real(np.fft.ifft2(FT_B/FT_A))

    core = K / sum(K)

    outHDU = pyfits.PrimaryHDU(data=core)
    outHDUList = pyfits.HDUList([outHDU])
    outHDUList.writeto('psf.fits',clobber=True)


    # Convolve the input image with found kernel K:
    convolvedImage = convolve(image, core)
    
    # Delete old file with convolved data if it exists:
    if os.path.exists(B_image):
        os.remove(B_image)
        
    # Create the output file with convolved data
    outHDU = pyfits.PrimaryHDU(data=convolvedImage, header=header)
    outHDUList = pyfits.HDUList([outHDU])
    outHDUList.writeto(B_image,clobber=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Common-Resolution Convolution")
    parser.add_argument("input_image", help="Input image A")
    parser.add_argument("output_image", help="Input name of the output image which will be created to have a common resolution with image B")
    parser.add_argument("psf_A_image", help="Input psf kernel (image) for image A")
    parser.add_argument("psf_B_image", help="Input psf kernel (image) for image B")
  
    args = parser.parse_args()
    
    input_image = args.input_image
    output_image = args.output_image
    psf_A_image = args.psf_A_image
    psf_B_image = args.psf_B_image      

    main(input_image,output_image,psf_A_image,psf_B_image)  