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

def main(composed_model_file, mask_image=None, output_mask='outliers.fits', ouliers_level=0.7, sigma=None, reg_file=None):
    hdulist = pyfits.open(composed_model_file)
    data = hdulist[3].data
    header = hdulist[0].header
    ny,nx = np.shape(data)
    
    if reg_file is not None:
        convert_reg_to_mask.mask(composed_model_file, reg_file, output_image=None, output_mask=reg_file.split('.reg')[0]+'_mask.reg', mask_value=1, mask_DN=None, verbosity=True)
        hdulist1 = pyfits.open(reg_file.split('.reg')[0]+'_mask.reg')
        reg_mask = hdulist1[0].data
    else:
        reg_mask = np.ones((ny,nx))
    
    
    if sigma is not None:
        data = ndimage.gaussian_filter(data, sigma=sigma, order=0)
    
    
    
    outliers_mask = np.zeros((ny,nx))

    if mask_image is not None:
        hdulist_mask = pyfits.open(mask_image)
        data_mask = hdulist_mask[0].data   
    else:
        data_mask = np.zeros((ny,nx))
    
    for k in range(ny):
        for i in range(nx):
            if data[k,i]>=ouliers_level and data_mask[k,i]==0 and data[k,i]<0.9:
                outliers_mask[k,i] = 1.*reg_mask[k,i]

    outHDU = pyfits.PrimaryHDU(outliers_mask, header=header)
    outHDU.writeto(output_mask, overwrite=True)
    print('Done!')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("composed_model_file", help="Name of the composed image")
    parser.add_argument("--mask", help="Name of the mask",
                        type=str, default=None)
    parser.add_argument("--o", help="Name of the output mask",
                        type=str, default='outliers.fits')
    parser.add_argument("--level", help="Outliers level (1=100%)",
                        type=float, default=0.6)
    parser.add_argument("--reg_file", help="Name of the region file with where the outliers should be found",
                        type=str, default=None)
    
    args = parser.parse_args()
    composed_model_file = args.composed_model_file
    mask_image = args.mask
    output_mask = args.o
    ouliers_level = args.level
    reg_file = args.reg_file
    
    main(composed_model_file, mask_image=mask_image, output_mask=output_mask, ouliers_level=ouliers_level, reg_file = reg_file)
            