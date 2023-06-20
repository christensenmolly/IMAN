#!/usr/bin/python
# DESCRIPTION:
# Script to check the background. A mask is required.
# MINIMAL USAGE: python check_background.py [input_image] [input_mask]

import warnings

import numpy as np
from astropy.io import fits as pyfits

warnings.filterwarnings("ignore")
from astropy.stats import sigma_clipped_stats
import scipy.ndimage as ndimage
import subprocess
import argparse

def convert_segm_to_boolean(mask):
    return np.ma.make_mask(mask, shrink=False)


def main(input_image, mask_image, hdu_inp=0, sigma_smooth=None):
    print('Checking the background...')
    hdulist = pyfits.open(input_image)
    data = hdulist[hdu_inp].data
    header = hdulist[hdu_inp].header
    ny, nx = np.shape(data)

    hdulist_mask = pyfits.open(mask_image)
    mask = hdulist_mask[0].data

    # Create astropy mask (boolean)
    mask_astropy = convert_segm_to_boolean(mask)

    # Find average sky
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask_astropy)
    # mmode = float(stats.mode(data, axis=None)[0][0]) # Too long!

    print('Sky for %s: %f (%f) +/- %f' % (input_image, mean, median, std))
    print('\t (in rms): %.2f (%.2f) +/- %.2f rms' % (mean / std, median / std, std / std))
    print('\t (in %% of rms): %.2f (%.2f) +/- %.2f %%' % (mean * 100. / std, median * 100. / std, std * 100. / std))

    # Apply gaussian smoothing for the unmasked data

    if sigma_smooth is not None:
        try:
            smooth_data = ndimage.gaussian_filter(data, sigma=sigma_smooth, order=0)
        except:
            smooth_data = ndimage.gaussian_filter(data, sigma=(sigma_smooth, sigma_smooth, 0), order=0)

    smooth_data = (smooth_data * ~mask_astropy) / std

    outHDU = pyfits.PrimaryHDU(smooth_data, header=header)
    outHDU.writeto('tmp.fits', clobber=True)
    subprocess.call("ds9 %s -scale histequ %s -scale histequ -cmap i8" % (input_image, 'tmp.fits'), shell=True)
    return mean, median, std



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fitting sky background with a polynomial")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("inputMask", help="Input mask image",type=str) 
    parser.add_argument("--hdu", help="Optional: HDU layer", type=int, default=0)
    parser.add_argument("--sigma", help="Optional: Sigma for Gaussian smoothing", type=float, default=None) 
    
    args = parser.parse_args()

    input_image = args.inputImage
    mask_image = args.inputMask
    hdu = args.hdu
    sigma_smooth = args.sigma

    main(input_image, mask_image, hdu_inp=hdu, sigma_smooth=sigma_smooth)

