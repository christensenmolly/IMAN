#!/usr/bin/env python

from math import sin
from math import cos
from math import pi
from math import modf
from math import sqrt
import numpy as np
import argparse
from types import SimpleNamespace

from astropy.io import fits as pyfits
from photutils.isophote import IsophoteList
from photutils.isophote import build_ellipse_model
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from scipy.signal import find_peaks, find_peaks_cwt
from scipy.signal import argrelextrema

import argparse
import convert_reg_to_mask

def main(mask_image, output_mask='outliers.fits', sigma=5.):
    hdulist = pyfits.open(mask_image)
    data = hdulist[0].data
    header = hdulist[0].header
    
    if sigma is not None:
        data = ndimage.gaussian_filter(data, sigma=sigma, order=0)
    
    data[data>0.]=1.

    outHDU = pyfits.PrimaryHDU(data, header=header)
    outHDU.writeto(output_mask, overwrite=True)
    print('Done!')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("mask", help="Name of the mask",
                        type=str)
    parser.add_argument("--o", help="Name of the output mask",
                        type=str, default='smoothed_mask.fits')
    parser.add_argument("--sigma", help="Sigma for smoothing",
                        type=float, default=5.)

    
    args = parser.parse_args()
    mask_image = args.mask
    output_mask = args.o
    sigma = args.sigma
    
    main(mask_image=mask_image, output_mask=output_mask, sigma=sigma)
            