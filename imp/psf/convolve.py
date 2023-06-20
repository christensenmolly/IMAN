#! /usr/bin/env python

import sys
from os.path import exists
from os import remove
from scipy.ndimage.filters import convolve
from scipy import signal
from numpy import max, sum
import astropy.io.fits as pyfits
from astropy.convolution import convolve, convolve_fft

def convolution(imageFileName, coreFileName, convolvedFileName):
    """Function performs convolution of image (fits file) 
    with given core (fits file as well)"""
    print(imageFileName)
    imageHDU = pyfits.open(imageFileName)[0]
    image = imageHDU.data#[0]
    header = imageHDU.header
    core = pyfits.open(coreFileName)[0].data
    core = core / sum(core)
    #convolvedImage = convolve(image, core)
    #convolvedImage = signal.fftconvolve(image, core, mode='same')
    convolvedImage = convolve_fft(image, core, crop=True, allow_huge=True) # Was changed!
    #convolvedImage = signal.convolve2d(image, core, boundary='symm', mode='same')
    # delete old file with convolved data if it exists
    if exists(convolvedFileName):
        remove(convolvedFileName)
    # creation of the file with convolved data
    outHDU = pyfits.PrimaryHDU(data=convolvedImage, header=header)
    outHDUList = pyfits.HDUList([outHDU])
    outHDUList.writeto(convolvedFileName)


if __name__ == "__main__":
    convolution(sys.argv[1],
                sys.argv[2],
                sys.argv[3])

    
